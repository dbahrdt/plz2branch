#include "State.h"
#include <osmtools/AreaExtractor.h>
#include <osmtools/AreaExtractorFilters.h>
#include <sserialize/algorithm/utilmath.h>
#include <queue>
#include <map>
#include "FileFormat.h"

#include <boost/property_map/property_map.hpp>
#include <boost/graph/directed_graph.hpp>
#include <boost/graph/successive_shortest_path_nonnegative_weights.hpp>
#include <boost/graph/find_flow_cost.hpp>


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
State::nodeDistances(memgraph::Graph::NodeId const & nId) const {
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
	auto weight = [](Graph::Edge const & e) { return double(e.distance)/e.speed; };
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
				if (qe.distance + weight(edge) < distance.at(edge.target)) {
					pq.emplace(Graph::NodeId(edge.target), qe.distance + weight(edge));
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
	std::vector<double> nodeDist = nodeDistances(branch.nodeId);
	double reachable = std::count_if(nodeDist.begin(), nodeDist.end(), [](auto && x) { return x != std::numeric_limits<double>::max(); });
	std::cout << "Branch " << branchId.value << " reaches " << reachable << "/" << nodeDist.size() << "="
			<< reachable/nodeDist.size()*100 << '%' << " nodes" << std::endl;
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
		std::vector<sserialize::AtomicMin<double>> distances(regionInfo.size());
		#pragma omp parallel for schedule(dynamic)
		for(std::size_t i=0; i < nodeDist.size(); ++i) {
			if (!nodeSelector(NodeId(i))) {
				continue;
			}
			RegionId nodeRegion = node2Region.at(i);
			if (!nodeRegion.valid()) {
				continue;
			}
			if (i == branch.nodeId) {
				std::cout << "Node distance for branch " << branchId.value << ": " << nodeDist[i] << std::endl;
			}
			distances.at(nodeRegion.value()).update(nodeDist[i]);
		}
		for(std::size_t i(0), s(branchDistance.size()); i < s; ++i) {
			branchDistance[i].value = distances[i].load();
		}
	}
	break;
	case DistanceWeightConfig::NWM::Max:
	{
		std::vector<sserialize::AtomicMax<double>> distances(regionInfo.size());
		#pragma omp parallel for schedule(dynamic)
		for(std::size_t i=0; i < nodeDist.size(); ++i) {
			if (!nodeSelector(NodeId(i))) {
				continue;
			}
			RegionId nodeRegion = node2Region.at(i);
			if (!nodeRegion.valid()) {
				continue;
			}
			distances.at(nodeRegion.value()).update(nodeDist[i]);
		}
		for(std::size_t i(0), s(branchDistance.size()); i < s; ++i) {
			branchDistance[i].value = distances[i].load();
		}
	}
	break;
	case DistanceWeightConfig::NWM::Mean:
	{
		std::vector< std::atomic<double> > distances(regionInfo.size());
		std::vector< std::atomic<uint32_t> > numNodes(regionInfo.size());
		
		#pragma omp parallel for schedule(dynamic)
		for(std::size_t i=0; i < nodeDist.size(); ++i) {
			if (!nodeSelector(NodeId(i))) {
				continue;
			}

			RegionId nodeRegion = node2Region.at(i);
			if (!nodeRegion.valid()) {
				continue;
			}
			if (nodeDist.at(i) != std::numeric_limits<double>::max()) {
				distances.at(nodeRegion.value()).fetch_add(nodeDist.at(i), std::memory_order_relaxed);
				numNodes.at(nodeRegion.value()).fetch_add(1, std::memory_order_relaxed);
			}
		}
		for(std::size_t i(0), s(branchDistance.size()); i < s; ++i) {
			branchDistance[i].value = distances[i].load()/numNodes.at(i).load();
		}
	}
	break;
	case DistanceWeightConfig::NWM::Median:
	{
		std::vector<sserialize::GuardedVariable<std::vector<double>>> distances(regionInfo.size());
		
		#pragma omp parallel for schedule(dynamic)
		for(std::size_t i=0; i < nodeDist.size(); ++i) {
			if (!nodeSelector(NodeId(i))) {
				continue;
			}
			
			RegionId nodeRegion = node2Region.at(i);
			if (!nodeRegion.valid()) {
				continue;
			}
			auto nd = nodeDist.at(i);
			
			if (nd != std::numeric_limits<double>::max()) {
				distances.at(nodeRegion.value()).syncedWithoutNotify([&](auto & v) { v.push_back(nd); });
			}
		}
		for(std::size_t i(0), s(branchDistance.size()); i < s; ++i) {
			branchDistance[i].value = sserialize::statistics::median(distances.at(i).unsyncedValue().begin(), distances.at(i).unsyncedValue().end(), double(0));
		}
	}
	break;
	}//end switch
	
	{
		RegionId brRId = node2Region.at(branch.nodeId);
		if (brRId.valid()) {
			std::cout << "Branch " << branchId.value << " is in plz " << regionInfo.at(brRId.value()).plz << " with distance " << branchDistance.at(brRId.value()).value  << std::endl;
		}
		else {
			std::cout << "Branch " << branchId.value << " is outside of all known plz" << std::endl;
		}
	}
	//reweight with inhabitants model
	switch (dwc.regionWeightModel) {
		case DistanceWeightConfig::RWM::Inhabitants:
			for(std::size_t i(0), s(branchDistance.size()); i < s; ++i) {
				branchDistance.at(i).value *= regionInfo.at(i).inhabitants;
			}
			break;
		case DistanceWeightConfig::RWM::Equal:
			break;
	};
	//reweight with number of employees
	switch (dwc.branchWeightModel) {
		case DistanceWeightConfig::BWM::Employees:
			for(auto & x : branchDistance) {
				x.value /= branch.employees;
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
			double reachable = std::count_if(branchDistances.at(i).begin(), branchDistances.at(i).end(), [](auto && x) { return x.value != std::numeric_limits<double>::max();});
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
				double branchPlzDist = branchDistances.at(brId).at(plzId).value;
				if (branchPlzDist < bestDist) {
					bestId = brId;
					bestDist = branchPlzDist;
				}
			}
			if (bestId != std::numeric_limits<uint32_t>::max()) {
				branches.at(bestId).assignedRegions.emplace_back(RegionId{plzId}, Distance{.value=bestDist});
			}
			else {
				std::cout << "Could not find branch for plz " << regionInfo.at(plzId).plz << std::endl;
				std::cout << "Branch 0 distance to plz: " << branchDistances.at(0).at(plzId).value << std::endl;
			}
	// 		pinfo(plzId);
		}
		pinfo.end();
	}
		break;
	case DistanceWeightConfig::BAM::Evolutionary:
	{
		emit_textInfo("Evolutionary algo not implemented");
		break;
	}
	case DistanceWeightConfig::BAM::Optimal:
	{
		//We map all weights to uint64_t with totalWeights mapping to std::numeric_limits<uint64_t>/2
		double maxWeight = 0;
		for(auto const & x : branchDistances) {
			for(auto const & y : x) {
				maxWeight = std::max(maxWeight, y.value);
			}
		}
		emit_textInfo(QString("Max weight: %1").arg(maxWeight));
		auto intWeight = [maxWeight](double v) -> int64_t {
			int64_t result = (v/maxWeight) * double(std::numeric_limits<uint32_t>::max());
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
				weights[e] = intWeight( branchDistances.at(brId).at(rId).value );
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
			br.cost += dist.value;
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
	std::stringstream ls;
	while (!data.eof() && data.good()) {
		std::string line;
		std::getline(data, line);
		if (!line.size()) {
			continue;
		}
		ls << line;
		std::string name;
		double lat, lon;
		uint32_t employees;
		ls >> name >> lat >> lon >> employees;
		//Ignore plz
	}
	#pragma omp parallel for schedule(dynamic)
	for(std::size_t i=0; i < tmp.size(); ++i) {
		std::apply([this](auto &&... params) { createBranch(params...); }, tmp[i]);
	}
}

void
State::emit_textInfo(QString const & str) {
	std::cout << str.toStdString() << std::endl;
// 	emit textInfo(str);
}

} //end namespace plz2branch
