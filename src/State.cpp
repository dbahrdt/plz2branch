#include "State.h"
#include <osmtools/AreaExtractor.h>
#include <osmtools/AreaExtractorFilters.h>
#include <sserialize/algorithm/utilmath.h>
#include <queue>

namespace plz2branch {

void registerTypesWithQt() {
	qRegisterMetaType<plz2branch::BranchId>();
	qRegisterMetaType<plz2branch::DistanceWeightConfig>();
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
		case DistanceWeightConfig::All:
			return true;
		case DistanceWeightConfig::OnlyCrossRoads:
			return graph.node(nid).edgeCount() > 1;
		};
		return true;
	};
		
	switch (dwc.weightModel) {
	case DistanceWeightConfig::Min:
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
	case DistanceWeightConfig::Max:
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
	case DistanceWeightConfig::Mean:
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
	case DistanceWeightConfig::Median:
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
	emit_textInfo("Finished computing closest branch for each PLZ");
	//Color the branches
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
	for(std::size_t i(0), s(branches.size()); i < s; ++i) {
		branches.at(i).color = QColor(colors.at(i%colors.size()));
	}
	emit_textInfo("Finished computing branch assignments");
}

void
State::writeBranchAssignments(std::ostream & out) {
	std::shared_lock<std::shared_mutex> lck(dataMtx);
	out << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for(std::size_t i(0), s(branches.size()); i < s; ++i) {
		Branch const & branch = branches[i];
		out << "{\n\tid:" << i << ",\n";
		out << "\nname:" << branch.name.toStdString() << ",\n";
		out << "\tcoord: (" << branch.coord.lat() << ", " << branch.coord.lon() << "),\n";
		out << "\tplz: [";
		for(auto const & [rId, rDist] : branch.assignedRegions) {
			RegionInfo const & ri = regionInfo.at(rId.value());
			out << ri.plz << ", ";
		}
		out << "]\n}";
	}
}

void
State::createBranch(double lat, double lon, QString name) {
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
		branches.push_back(branch);
	}
	emit dataChanged();
}

void
State::createBranch(double lat, double lon) {
	createBranch(lat, lon, "");
}

void
State::createBranches(std::istream & data) {
	while (!data.eof() && data.good()) {
		std::string name;
		double lat, lon;
		data >> name >> lat >> lon;
		if (data.good()) {
			createBranch(lat, lon, QString::fromStdString(name));
		}
	}
}

void
State::emit_textInfo(QString const & str) {
	std::cout << str.toStdString() << std::endl;
// 	emit textInfo(str);
}

} //end namespace plz2branch
