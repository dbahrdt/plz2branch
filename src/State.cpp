#include "State.h"
#include <osmtools/AreaExtractor.h>
#include <osmtools/AreaExtractorFilters.h>
#include <sserialize/algorithm/utilmath.h>
#include <sserialize/strings/stringfunctions.h>
#include <queue>
#include <map>
#include <random>
#include "FileFormat.h"

#include <boost/property_map/property_map.hpp>
#include <boost/graph/directed_graph.hpp>
#include <boost/graph/successive_shortest_path_nonnegative_weights.hpp>
#include <boost/graph/find_flow_cost.hpp>

#include <pcg_random.hpp>


namespace plz2branch {

void registerTypesWithQt() {
	qRegisterMetaType<plz2branch::BranchId>();
	qRegisterMetaType<plz2branch::DistanceWeightConfig>();
	qRegisterMetaType<plz2branch::DistanceWeightConfig::NS>();
	qRegisterMetaType<plz2branch::DistanceWeightConfig::NWM>();
	qRegisterMetaType<plz2branch::DistanceWeightConfig::RWM>();
	qRegisterMetaType<plz2branch::DistanceWeightConfig::BWM>();
	qRegisterMetaType<plz2branch::DistanceWeightConfig::BAM>();
}

void
State::importFiles() {
	osmtools::OsmGridRegionTree<RegionInfo> grt;
	osmpbf::RCFilterPtr filter(
		new osmpbf::AndTagFilter({
			new osmpbf::KeyValueTagFilter("boundary", "postal_code"),
			new osmpbf::KeyOnlyTagFilter("postal_code")
		})
	);
	osmtools::AreaExtractor ae;
	auto processor = [&](const std::shared_ptr<sserialize::spatial::GeoRegion> & region, osmpbf::IPrimitive & primitive) {
		RegionInfo ri = { .bbox = region->boundary(), .plz = (uint32_t) std::atoi( primitive.valueByKey("postal_code").c_str() ) };
		grt.push_back(*region, ri);
	};
	std::cout << "Fetching regions" << std::endl;
	ae.extract(cfg.osmFile, processor, ae.ET_ALL_SPECIAL_BUT_BUILDINGS, filter);
	std::cout << "Found " << grt.regions().size() << " regions" << std::endl;
	
	std::cout << "Computing grt" << std::endl;
	grt.setRefinerOptions(2, 2, 250);
	grt.addPolygonsToRaster(100, 100);
	
	std::cout << "Fetching graph" << std::endl;
	graph = memgraph::Graph::fromPBF(cfg.osmFile, memgraph::Graph::Edge::AT_CAR);
	
	struct It {
		RegionId * dest;
		It & operator++() { return *this; }
		It & operator*() { return *this; }
		It & operator=(uint32_t rId) { *dest = RegionId(rId); return *this;}
	};
	
	node2Region.resize(graph.nodeCount());
	
	std::cout << "Computing Node <-> Region mapping" << std::endl;
	#pragma omp parallel for schedule(dynamic)
	for(std::size_t i=0; i < graph.nodeCount(); ++i) {
		auto const & ni = graph.nodeInfo(i);
		It out{.dest=&node2Region.at(i)};
		grt.find(ni.lat, ni.lon, out);
	}
	
	regionInfo = std::move(grt.values());
	
	if (cfg.plzFile.size()) {
		std::cout << "Importing plz inhabitants file" << std::endl;
		std::unordered_map<uint32_t, uint32_t> plz2in;
		try {
			plz2in = plz2inhabitants(cfg.plzFile);
		}
		catch (std::runtime_error const & e) {
			std::cerr << "Could not import plz file: " << e.what() << std::endl;
		}
		for(auto & x : regionInfo) {
			try {
				x.inhabitants = plz2in.at(x.plz);
			}
			catch (std::out_of_range const &) {
				std::cerr << "Did not find plz " << x.plz << " in plz2inhabitants file. Calculations will be flawed" << std::endl;
			}
		}
	}
	
	std::cout << "Preprocessing complete" << std::endl;
}

memgraph::Graph::NodeId
State::closestNode(sserialize::spatial::GeoPoint const & gp) const {
	sserialize::spatial::DistanceCalculator dc(sserialize::spatial::DistanceCalculator::DCT_GEODESIC_FAST);
	
	memgraph::Graph::NodeId bestId;
	double bestDist = std::numeric_limits<double>::max();
	for(uint32_t i(0), s(graph.nodeCount()); i < s; ++i) {
		auto const & ni = graph.nodeInfo(i);
		double dist = std::abs(dc.calc(ni.lat, ni.lon, gp.lat(), gp.lon()));
		if (dist < bestDist) {
			bestDist = dist;
			bestId = memgraph::Graph::NodeId(i);
		}
	}
	return bestId;
}

std::vector<double>
State::nodeDistances(memgraph::Graph::NodeId const & nId, std::function<double(memgraph::Graph::Edge const &)> weight) const {
	using namespace memgraph;
	constexpr double NoDistance = std::numeric_limits<double>::max();
	std::vector<double> distance(graph.nodeCount(), NoDistance);
	struct QueueElement {
		Graph::NodeId nId;
		double distance{std::numeric_limits<double>::max()};
		QueueElement(Graph::NodeId const & nid, double distance) : nId(nid), distance(distance) {}
		QueueElement(QueueElement const &) = default;
		QueueElement & operator=(QueueElement const&) = default;
		bool operator<(QueueElement const & other) const {
			return !(distance < other.distance);
		}
	};
	std::priority_queue<QueueElement> pq;
	pq.emplace(nId, 0);
	while(pq.size()) {
		auto qe = pq.top();
		pq.pop();
		if (qe.distance < distance.at(qe.nId)) { //current distance is smaller than recorded distance
			distance.at(qe.nId) = qe.distance;
			//check all edges and insert the ones with smaller updated distance in our queue
			for(auto eIt(graph.edgesBegin(qe.nId)), eEnd(graph.edgesEnd(qe.nId)); eIt != eEnd; ++eIt) {
				auto const & edge = *eIt;
				assert(edge.source == qe.nId);
				if (auto nd = qe.distance + weight(edge); nd < distance.at(edge.target)) {
					pq.emplace(Graph::NodeId(edge.target), nd);
				}
			}
		}
	}
	return distance;
}

std::vector<Distance>
State::branchDistance(BranchId const & branchId, DistanceWeightConfig const & dwc) const {
	Branch const & branch = branches.at(branchId.value);
	using NodeId = memgraph::Graph::NodeId;
	//First get the node distances
	std::vector<double> timeNodeDist = nodeDistances(branch.nodeId, [](memgraph::Graph::Edge const & e) -> double { return double(e.distance)/e.speed; });
	std::vector<double> spatialNodeDist = nodeDistances(branch.nodeId, [](memgraph::Graph::Edge const & e) -> double { return e.distance; });
	double reachable = std::count_if(timeNodeDist.begin(), timeNodeDist.end(), [](auto && x) { return x != std::numeric_limits<double>::max(); });
	std::cout << "Branch " << branchId.value << " reaches " << reachable << "/" << timeNodeDist.size() << "="
			<< reachable/timeNodeDist.size()*100 << '%' << " nodes" << std::endl;
	std::vector<Distance> branchDistance(regionInfo.size());
	
	auto nodeSelector = [&](NodeId const & nid) {
		switch (dwc.nodeSelection) {
		case DistanceWeightConfig::NS::All:
			return true;
		case DistanceWeightConfig::NS::OnlyCrossRoads:
			return graph.node(nid).edgeCount() > 1;
		};
		return true;
	};
	
	switch (dwc.nodeWeightModel) {
	case DistanceWeightConfig::NWM::Min:
	{
		std::vector<sserialize::AtomicMin<double>> timeD(regionInfo.size());
		std::vector<sserialize::AtomicMin<double>> spatialD(regionInfo.size());
		#pragma omp parallel for schedule(dynamic)
		for(std::size_t i=0; i < timeNodeDist.size(); ++i) {
			if (!nodeSelector(NodeId(i))) {
				continue;
			}
			RegionId nodeRegion = node2Region.at(i);
			if (!nodeRegion.valid()) {
				continue;
			}
			if (i == branch.nodeId) {
				std::cout << "Node distance for branch " << branchId.value << ": " << timeNodeDist[i] << std::endl;
			}
			timeD.at(nodeRegion.value()).update(timeNodeDist[i]);
			spatialD.at(nodeRegion.value()).update(spatialNodeDist[i]);
		}
		for(std::size_t i(0), s(branchDistance.size()); i < s; ++i) {
			branchDistance[i].time = timeD[i].load();
			branchDistance[i].spatial = spatialD[i].load();
		}
	}
	break;
	case DistanceWeightConfig::NWM::Max:
	{
		std::vector<sserialize::AtomicMax<double>> timeD(regionInfo.size());
		std::vector<sserialize::AtomicMax<double>> spatialD(regionInfo.size());
		#pragma omp parallel for schedule(dynamic)
		for(std::size_t i=0; i < timeNodeDist.size(); ++i) {
			if (!nodeSelector(NodeId(i))) {
				continue;
			}
			RegionId nodeRegion = node2Region.at(i);
			if (!nodeRegion.valid()) {
				continue;
			}
			timeD.at(nodeRegion.value()).update(timeNodeDist[i]);
			spatialD.at(nodeRegion.value()).update(spatialNodeDist[i]);
		}
		for(std::size_t i(0), s(branchDistance.size()); i < s; ++i) {
			branchDistance[i].time = timeD[i].load();
			branchDistance[i].spatial = spatialD[i].load();
		}
	}
	break;
	case DistanceWeightConfig::NWM::Mean:
	{
		std::vector< std::atomic<double> > timeD(regionInfo.size());
		std::vector< std::atomic<double> > spatialD(regionInfo.size());
		std::vector< std::atomic<uint32_t> > numNodes(regionInfo.size());
		
		#pragma omp parallel for schedule(dynamic)
		for(std::size_t i=0; i < timeNodeDist.size(); ++i) {
			if (!nodeSelector(NodeId(i))) {
				continue;
			}

			RegionId nodeRegion = node2Region.at(i);
			if (!nodeRegion.valid()) {
				continue;
			}
			if (timeNodeDist.at(i) != std::numeric_limits<double>::max()) {
				timeD.at(nodeRegion.value()).fetch_add(timeNodeDist.at(i), std::memory_order_relaxed);
				spatialD.at(nodeRegion.value()).fetch_add(spatialNodeDist.at(i), std::memory_order_relaxed);
				numNodes.at(nodeRegion.value()).fetch_add(1, std::memory_order_relaxed);
			}
		}
		for(std::size_t i(0), s(branchDistance.size()); i < s; ++i) {
			branchDistance[i].time = timeD[i].load()/numNodes.at(i).load();
			branchDistance[i].spatial = spatialD[i].load()/numNodes.at(i).load();
		}
	}
	break;
	case DistanceWeightConfig::NWM::Median:
	{
		std::vector<sserialize::GuardedVariable<std::vector<double>>> timeD(regionInfo.size());
		std::vector<sserialize::GuardedVariable<std::vector<double>>> spatialD(regionInfo.size());
		
		#pragma omp parallel for schedule(dynamic)
		for(std::size_t i=0; i < timeNodeDist.size(); ++i) {
			if (!nodeSelector(NodeId(i))) {
				continue;
			}
			
			RegionId nodeRegion = node2Region.at(i);
			if (!nodeRegion.valid()) {
				continue;
			}
			auto nd = timeNodeDist.at(i);
			auto sd = spatialNodeDist.at(i);
			
			if (nd != std::numeric_limits<double>::max()) {
				timeD.at(nodeRegion.value()).syncedWithoutNotify([&](auto & v) { v.push_back(nd); });
				spatialD.at(nodeRegion.value()).syncedWithoutNotify([&](auto & v) { v.push_back(sd); });
			}
		}
		for(std::size_t i(0), s(branchDistance.size()); i < s; ++i) {
			branchDistance[i].time = sserialize::statistics::median(timeD.at(i).unsyncedValue().begin(), timeD.at(i).unsyncedValue().end(), double(0));
			branchDistance[i].spatial = sserialize::statistics::median(spatialD.at(i).unsyncedValue().begin(), spatialD.at(i).unsyncedValue().end(), double(0));
		}
	}
	break;
	}//end switch
	
	{
		RegionId brRId = node2Region.at(branch.nodeId);
		if (brRId.valid()) {
			std::cout << "Branch " << branchId.value << " is in plz " << regionInfo.at(brRId.value()).plz << " with distance " << branchDistance.at(brRId.value()).time  << std::endl;
		}
		else {
			std::cout << "Branch " << branchId.value << " is outside of all known plz" << std::endl;
		}
	}
	//reweight with inhabitants model
	switch (dwc.regionWeightModel) {
		case DistanceWeightConfig::RWM::Inhabitants:
			for(std::size_t i(0), s(branchDistance.size()); i < s; ++i) {
				branchDistance.at(i).time *= regionInfo.at(i).inhabitants;
			}
			break;
		case DistanceWeightConfig::RWM::Equal:
			break;
	};
	//reweight with number of employees
	switch (dwc.branchWeightModel) {
		case DistanceWeightConfig::BWM::Employees:
			for(auto & x : branchDistance) {
				x.time /= branch.employees;
			}
			break;
		case DistanceWeightConfig::BWM::Equal:
			break;
	};
	
	return branchDistance;
}

void
State::computeBranchAssignments(DistanceWeightConfig const & dwc) {
	emit_textInfo("Computing branch distances");
	std::vector<std::vector<Distance>> branchDistances(branches.size());
	{
		std::shared_lock<std::shared_mutex> lck(dataMtx);
		#pragma omp parallel for schedule(dynamic)
		for(std::size_t i = 0; i < branches.size(); ++i) {
			Branch const & branch = branches.at(i);
			branchDistances.at(i) = branchDistance(BranchId{.value=uint32_t(i)}, dwc);
			double reachable = std::count_if(branchDistances.at(i).begin(), branchDistances.at(i).end(), [](auto && x) { return x.time != std::numeric_limits<double>::max();});
			std::cout << "Branch " << i << " reaches " << reachable << "/" << branchDistances.at(i).size() << "="
			<< reachable/branchDistances.at(i).size()*100 << '%' << " plz" << std::endl;
			emit_textInfo(QString("Finished computing distances for branch %1").arg(i));
		}
	}
	emit_textInfo("Finished computing branch distances");
	std::unique_lock<std::shared_mutex> lck(dataMtx);
	
	for(Branch & branch : branches) {
		branch.assignedRegions.clear();
	}
	
	switch (dwc.branchAssignmentModel) {
	case DistanceWeightConfig::BAM::ByDistance:
	{
		emit_textInfo("Computing closest branch for each PLZ");
		//For each plz, compute the branch that is closest
		sserialize::ProgressInfo pinfo;
		pinfo.begin(regionInfo.size());
		for(uint32_t plzId(0), s(regionInfo.size()); plzId < s; ++plzId) {
			double bestDist = std::numeric_limits<double>::max();
			uint32_t bestId = std::numeric_limits<uint32_t>::max();
			for(uint32_t brId(0); brId < branches.size(); ++brId) {
				double branchPlzDist = branchDistances.at(brId).at(plzId).time;
				if (branchPlzDist < bestDist) {
					bestId = brId;
					bestDist = branchPlzDist;
				}
			}
			if (bestId != std::numeric_limits<uint32_t>::max()) {
				assert(branchDistances.at(bestId).at(plzId).time == bestDist);
				branches.at(bestId).assignedRegions.emplace_back(RegionId{plzId}, branchDistances.at(bestId).at(plzId));
			}
			else {
				std::cout << "Could not find branch for plz " << regionInfo.at(plzId).plz << std::endl;
				std::cout << "Branch 0 distance to plz: " << branchDistances.at(0).at(plzId).time << std::endl;
			}
	// 		pinfo(plzId);
		}
		pinfo.end();
	}
		break;
	case DistanceWeightConfig::BAM::Evolutionary:
	{
		using Resident = std::vector<uint16_t>;
		using Population = std::vector<Resident>;
		pcg_extras::seed_seq_from<std::random_device> seed_source;
		pcg32 rg(seed_source);
		auto u01d = std::uniform_real_distribution<double>(0,1);
		auto bd = std::uniform_int_distribution<uint16_t>(0, branches.size()-1);
		auto sd = std::geometric_distribution<uint32_t>(0.5); //take half of all upper
		Population pop;
		sserialize::MinMax<double> mm_residentCosts;
		std::vector<double> residentCosts;
		auto less = [&](Distance const & a, Distance const & b) {
			return a.time < b.time;
		};
		auto mutate = [&](Resident & r, double mutationRate) -> Resident {
			Resident result;
			for(auto & x : r) {
				if (u01d(rg) < mutationRate) {
					result.push_back( bd(rg) );
				}
				else {
					result.push_back(x);
				}
			}
			return result;
		};
		auto sex = [&](Resident & r1, Resident & r2) -> Resident {
			Resident child;
			for(std::size_t i(0), s(r1.size()); i < s; ++i) {
				Distance const & cost1 = branchDistances.at(r1.at(i)).at(i);
				Distance const & cost2 = branchDistances.at(r2.at(i)).at(i);
				if (less(cost1, cost2)) {
					child.push_back(r1.at(i));
				}
				else {
					child.push_back(r2.at(i));
				}
			}
			return child;
		};
		auto rcost = [&](Resident const & r) {
			double cost = 0;
			for(std::size_t rId(0); rId < regionInfo.size(); ++rId) {
				cost += branchDistances.at(r.at(rId)).at(rId).time;
			}
			return cost;
		};
		auto pcost = [&](Population const & p) -> double {
			double cost = 0;
			for(auto const & r : p) {
				cost += rcost(r);
			}
			return cost;
		};
		auto pbest = [&](Population const & p) -> std::pair<uint32_t, double> {
			uint32_t best_pos = std::numeric_limits<uint32_t>::max();
			double best_cost = std::numeric_limits<double>::max();
			for(std::size_t i(0), s(p.size()); i < s; ++i) {
				double cost = rcost(p.at(i));
				if (cost < best_cost) {
					best_cost = cost;
					best_pos = i;
				}
			}
			return std::make_pair(best_pos, best_cost);
		};
		auto updateResidentCosts = [&]() {
			mm_residentCosts.reset();
			residentCosts.resize(0);
			residentCosts.reserve(pop.size());
			for(auto const & r : pop) {
				residentCosts.push_back( rcost(r) );
				mm_residentCosts.update(residentCosts.back());
			}
		};
		bool allValid = std::count_if(branches.begin(), branches.end(), [&](auto && x) { return x.maxRegions < regionInfo.size(); }) > 0;
		auto valid = [&](Resident const & r) {
			if (allValid) {
				return true;
			}
			std::vector<uint16_t> counts(branches.size(), 0);
			for(auto const & x : r) {
				counts.at(x) += 1;
			}
			for(std::size_t i(0); i < counts.size(); ++i) {
				if (branches.at(i).maxRegions < counts.at(i)) {
					return false;
				}
			}
			return true;
		};
		//Select residents based on their cost relative to the best
		//The best is selected with propability 1, the worst with probability 0
		//The others are selected with probability (cost - worst)/(best - worst)
		auto select = [&](Population const & p) -> Population {
			Population result;
			for(std::size_t i(0), s(residentCosts.size()); i < s; ++i) {
				if (!valid(p.at(i))) { //remove invalid residents
					continue;
				}
				double prob = 1-(residentCosts.at(i)-mm_residentCosts.min())/(mm_residentCosts.max() - mm_residentCosts.min());
				if (u01d(rg) <= prob) {
					result.push_back( p.at(i) );
				}
			}
			return result;
		};
		//Fill our population

		std::size_t popSize = 1024;
		for(std::size_t i(0); i < popSize; ++i) {
			Resident r;
			for(std::size_t rId(0); rId < regionInfo.size(); ++rId) {
				r.push_back( bd(rg) );
			}
			pop.push_back(std::move(r));
		}
		//Now play evolution
		//We mutate a resident with probability 1-(cost - worst)/(best - worst)
		//Hence worst is mutated with probability 1 whereas the best is no mutated
		//
		//Additionally the residents may have sex with each other with probability?
		//The better ones should have a higher probability to have sex with each other than with a worse one?
		std::size_t numGenerations = 1000;
		double mutationRate = 0.01;
		updateResidentCosts();
		for(std::size_t generation(0); generation < numGenerations; ++generation) {
			emit_textInfo(QString("Generation %1: min=%2 max=%3").arg(generation).arg(mm_residentCosts.min()).arg(mm_residentCosts.max()));
			//Mutate
			for(std::size_t resId(0), s(pop.size()); resId < s; ++resId) {
				double prob = (residentCosts.at(resId)-mm_residentCosts.min())/(mm_residentCosts.max()-mm_residentCosts.min());
				if (u01d(rg) < prob) {
					pop.push_back( mutate(pop.at(resId), mutationRate) );
				}
			}
			emit_textInfo(QString("Generation %1: mutation: pop.size()=%2").arg(generation).arg(pop.size()));
			//Each should have sex with a random partner (and probability corresponding to their cost distance?)
			std::uniform_int_distribution<std::size_t> resRand(0, pop.size());
			for(std::size_t resId(0), s(pop.size()); resId < s; ++resId) {
				std::size_t partnerId = resRand(rg);
				pop.push_back( sex(pop.at(resId), pop.at(partnerId)) );
			}
			emit_textInfo(QString("Generation %1: sex: pop.size()=%2").arg(generation).arg(pop.size()));
			//Select the population down to the original population size
			updateResidentCosts();
			while(pop.size() > popSize && mm_residentCosts.min() != mm_residentCosts.max()) {
				emit_textInfo(QString("Generation %1: select: pop.size=%2 cost(min/max)=%3/%4)=").arg(generation).arg(pop.size()).arg(mm_residentCosts.min()).arg(mm_residentCosts.max()));
				pop = select(pop);
				updateResidentCosts();
			}
			emit_textInfo(QString("Generation %1: select: pop.size()=%2").arg(generation).arg(pop.size()));
		}
		//now get the best one and use it as a result
		auto best = pbest(pop);
		for(std::size_t rId(0); rId < regionInfo.size(); ++rId) {
			auto brId = pop.at(best.first).at(rId);
			branches.at(brId).assignedRegions.emplace_back(RegionId(rId), branchDistances.at(brId).at(rId));
		}
		break;
	}
	case DistanceWeightConfig::BAM::Optimal:
	{
		//We map all weights to uint64_t with totalWeights mapping to std::numeric_limits<uint64_t>/2
		double totalWeight = 0;
		for(auto const & x : branchDistances) {
			for(auto const & y : x) {
				totalWeight += y.time;
			}
		}
		emit_textInfo(QString("Total weight: %1").arg(totalWeight));
		auto intWeight = [totalWeight](double v) -> int64_t {
			int64_t result = (v/totalWeight) * double(uint64_t(1) << 58);
			if (v >= 1 && result == 0) {
				throw std::runtime_error("Precision of uint32_t not high enough");
			}
			return result;
		};
		using Graph = boost::adjacency_list<>;
		using vertex_descriptor = boost::graph_traits<Graph>::vertex_descriptor;
		using edge_descriptor = boost::graph_traits<Graph>::edge_descriptor;
		emit_textInfo("Constructing graph...");
		Graph g;
		vertex_descriptor source = boost::add_vertex(g);
		vertex_descriptor target = boost::add_vertex(g);
		std::vector<vertex_descriptor> branch_nodes;
		std::vector<edge_descriptor> s_branch_edges;
		std::map<edge_descriptor, int64_t> weights;
		std::map<edge_descriptor, int64_t> capacities;
		std::map<edge_descriptor, edge_descriptor> reverse_edges;
		for(auto const & branch : branches) {
			vertex_descriptor v = boost::add_vertex(g);
			branch_nodes.push_back(v);
			auto e = boost::add_edge(source, v, g).first;
			//real edge
			s_branch_edges.push_back(e);
			weights[e] = 0;
			capacities[e] = branch.maxRegions;
			//reverse edge
			auto re = boost::add_edge(v, source, g).first;
			weights[re] = -weights[e];
			capacities[re] = 0;
			reverse_edges[e] = re;
			reverse_edges[re] = e;
		}
		std::vector<vertex_descriptor> regionNodes;
		std::vector<edge_descriptor> branch_region_edges;
		std::vector<edge_descriptor> region_t_edges;
		for(std::size_t rId(0); rId < regionInfo.size(); ++rId) {
			vertex_descriptor region_v = boost::add_vertex(g);
			regionNodes.push_back( region_v );
			for(std::size_t brId(0); brId < branches.size(); ++brId) {
				auto const & branch_v = branch_nodes.at(brId);
				auto e = boost::add_edge(branch_v, region_v, g).first;
				branch_region_edges.push_back(e);
				weights[e] = intWeight( branchDistances.at(brId).at(rId).time );
				capacities[e] = 1;
				//reverse edge
				auto re = boost::add_edge(region_v, branch_v, g).first;
				weights[re] = -weights[e];
				capacities[re] = 0;
				reverse_edges[e] = re;
				reverse_edges[re] = e;
			}
			auto e = boost::add_edge(region_v, target, g).first;
			region_t_edges.push_back(e);
			weights[e] = 0;
			capacities[e] = 1;
			//reverse edge
			auto re = boost::add_edge(target, region_v, g).first;
			weights[re] = -weights[e];
			capacities[re] = 0;
			reverse_edges[e] = re;
			reverse_edges[re] = e;
		}
		emit_textInfo( QString("Graph has %1 vertices and %2 edges").arg(boost::num_vertices(g)).arg(boost::num_edges(g)) );
		emit_textInfo("Executing boost::successive_shortest_path_nonnegative_weights");
		std::map<edge_descriptor, int64_t> residual_capacities;
		boost::successive_shortest_path_nonnegative_weights(
			g,
			source,
			target,
			boost::weight_map(boost::associative_property_map(weights))
			.capacity_map(boost::associative_property_map(capacities))
			.residual_capacity_map(boost::associative_property_map(residual_capacities))
			.reverse_edge_map(boost::associative_property_map(reverse_edges))
		);
		emit_textInfo(QString("residual_capacities.size()=%1").arg(residual_capacities.size()));
		emit_textInfo("Assembling mapping");
		
		//Assemble the mapping:
		//The edges between the branches an the regions are the ones that we picked
		//An edge (branch, region) with flow=1 (residual capacity=0) means
		//that the branch handles the region
		uint32_t plzWithBranch = 0;
		auto edgeIt = branch_region_edges.begin();
		for(uint32_t rId(0); rId < regionInfo.size(); ++rId) {
			for(uint32_t brId(0); brId < branches.size(); ++brId, ++edgeIt) {
				assert(residual_capacities.at(*edgeIt) == 0 || residual_capacities.at(*edgeIt) == 1);
				if (residual_capacities.count(*edgeIt) && residual_capacities.at(*edgeIt) == 0) {
					branches.at(brId).assignedRegions.emplace_back(RegionId{rId}, branchDistances.at(brId).at(rId));
					++plzWithBranch;
				}
			}
		}
		assert(plzWithBranch == regionInfo.size());
		if (plzWithBranch != regionInfo.size()) {
			std::cerr << "Failed to compute assignment for " << regionInfo.size() - plzWithBranch << " regions" << std::endl;
			for(uint32_t brId(0); brId < branches.size(); ++brId) {
				auto const & br = branches.at(brId);
				auto flow = br.maxRegions - residual_capacities.at(s_branch_edges.at(brId));
				if (br.assignedRegions.size() !=  flow) {
					std::cerr << "Branch " << brId << " has " << br.assignedRegions.size() << " regions and flow of " << flow << std::endl;
				}
			}
		}
		emit_textInfo("Assembled mapping");
	}
		break;
	};
	emit_textInfo("Finished computing closest branch for each PLZ");

	emit_textInfo("Computing costs");
	totalCost = 0;
	for(auto & br : branches) {
		br.cost = 0;
		for(auto const & [_, dist] : br.assignedRegions) {
			br.cost += dist.time;
		}
		totalCost += br.cost;
	}
	emit_textInfo("Finished computing branch assignments");
}

void
State::writeBranchAssignments(std::ostream & out) {
	std::shared_lock<std::shared_mutex> lck(dataMtx);
	out << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for(std::size_t i(0), s(branches.size()); i < s; ++i) {
		Branch const & branch = branches[i];
		out << branch.name.toStdString() << '\t';
		out << branch.coord.lat() << '\t' << branch.coord.lon() << '\t';
		out << branch.employees << '\t';
		for(std::size_t j(0), js(branch.assignedRegions.size()); j < js;) {
			out << regionInfo.at(branch.assignedRegions.at(j).first.value()).plz;
			++j;
			if (j < js) {
				out << '\t';
			}
		}
		out << '\n';
	}
}

void
State::exportBranchAssignments(std::ostream & out) {
	std::shared_lock<std::shared_mutex> lck(dataMtx);
	auto prec = out.precision(std::numeric_limits<double>::digits10);
	out << "PLZ\tBranch\tDistance [m]\n";
	for(Branch const & branch : branches) {
		std::string bn = branch.name.toStdString();
		for(auto const & [rId, dist] : branch.assignedRegions) {
			auto plz = regionInfo.at(rId.value()).plz;
			if (plz < 10000) {
				out << '0';
			}
			out << plz << '\t' << bn << '\t' << dist.spatial/1000 << '\n';
		}
	}
	out.precision(prec);
}

void
State::createBranch(double lat, double lon, QString name, uint32_t employees) {
	std::array<Qt::GlobalColor, 8> colors{
		Qt::red,
		Qt::blue,
		Qt::darkRed,
		Qt::darkGreen,
		Qt::darkBlue,
		Qt::darkCyan,
		Qt::darkMagenta,
		Qt::darkYellow
	};
	Branch branch;
	{
		std::shared_lock<std::shared_mutex> lck;
		branch.nodeId = closestNode(sserialize::spatial::GeoPoint(lat, lon));
	}
	branch.name = name;
	branch.coord.lat() = graph.nodeInfo(branch.nodeId).lat;
	branch.coord.lon() = graph.nodeInfo(branch.nodeId).lon;
	{
		std::unique_lock<std::shared_mutex> lck;
		branch.color = QColor(colors.at(branches.size()%colors.size()));
		branches.push_back(branch);
	}
	emit dataChanged();
}

void
State::createBranch(double lat, double lon) {
	createBranch(lat, lon, "", 1);
}

void
State::createBranches(std::istream & data) {
	std::vector< std::tuple<double, double, QString, uint32_t> > tmp;
	std::size_t lineNum = 0;
	while (!data.eof() && data.good()) {
		++lineNum;
		std::string line;
		std::getline(data, line);
		if (!line.size()) {
			continue;
		}
		auto sl = sserialize::split< std::vector<std::string> >(line, '\t');
		if (sl.size() < 3) {
			std::cerr << "Ignoring line " << lineNum << std::endl;
			continue;
		}
		std::string name;
		double lat, lon;
		uint32_t employees = 1;
		name = sl.at(0);
		try {
			lat = sserialize::toDouble(sl.at(1));
			lon = sserialize::toDouble(sl.at(2));
			if (sl.size() > 3) {
				employees = std::stoi(sl.at(3));
			}
		}
		catch (std::exception const & e) {
			std::cerr << "Ignoring line " << lineNum << " due to exception: " << e.what() << std::endl;
		}
		//Ignore plz
		tmp.emplace_back(lat, lon, QString::fromStdString(name), employees);
	}
	#pragma omp parallel for schedule(dynamic)
	for(std::size_t i=0; i < tmp.size(); ++i) {
		std::apply([this](auto &&... params) { createBranch(params...); }, tmp[i]);
	}
}

void State::clear() {
	std::unique_lock<std::shared_mutex> lck(dataMtx);
	enabledBranches.clear();
	branches.clear();
	totalCost = 0;
}

void
State::emit_textInfo(QString const & str) {
	std::cout << str.toStdString() << std::endl;
// 	emit textInfo(str);
}

} //end namespace plz2branch
