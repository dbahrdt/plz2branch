#pragma once

#include <shared_mutex>
#include <unordered_set>
#include <memgraph/Graph.h>
#include <osmtools/OsmGridRegionTree.h>
#include <QColor>
#include <QString>
#include <QObject>

namespace plz2branch {

struct Config {
	Config() {}
	std::string osmFile;
};

struct RegionInfo {
	std::shared_ptr<sserialize::spatial::GeoRegion> shape;
	uint32_t plz{std::numeric_limits<uint32_t>::max()};
	inline bool valid() const { return plz != std::numeric_limits<uint32_t>::max(); }
};

struct RegionId {
	uint32_t value{std::numeric_limits<uint32_t>::max()};
	inline bool valid() const { return value != std::numeric_limits<uint32_t>::max(); }
};

struct BranchId {
	uint32_t value{std::numeric_limits<uint32_t>::max()};
	inline bool valid() const { return value != std::numeric_limits<uint32_t>::max(); }
};

struct DistanceWeightConfig {
	enum NodeSelection {
		All,
		OnlyCrossRoads
	};
	enum WeightModel {
		Min,
		Max,
		Median,
		Mean
	};
	DistanceWeightConfig() {}
	DistanceWeightConfig(NodeSelection ns, WeightModel wm) : nodeSelection(ns), weightModel(wm) {}
	DistanceWeightConfig(DistanceWeightConfig const&) = default;
	DistanceWeightConfig & operator=(DistanceWeightConfig const&) = default;
	NodeSelection nodeSelection{All};
	WeightModel weightModel{Min};
};

std::string to_string(DistanceWeightConfig::NodeSelection v);
std::string to_string(DistanceWeightConfig::WeightModel v);

void from_string(std::string const & src, DistanceWeightConfig::WeightModel & dest);
void from_string(std::string const & src, DistanceWeightConfig::NodeSelection & dest);

struct Distance {
	double value{std::numeric_limits<double>::max()};
};

struct Branch {
	QString name;
	QColor color{Qt::green};
	sserialize::spatial::GeoPoint coord;
	memgraph::Graph::NodeId nodeId;
	std::vector<std::pair<RegionId, Distance>> assignedRegions;
	std::vector<Distance> dist; //regionId->distance
};

class State: public QObject {
Q_OBJECT
public:
	State(QObject * parent) : QObject(parent) {}
	~State() override {}
public:
	Config cfg;
public:
	memgraph::Graph graph;
	std::vector<RegionInfo> regionInfo;
	std::vector<RegionId> node2Region;
	
	std::vector<RegionInfo> regions;
	std::vector<Branch> branches;
	std::unordered_set<uint32_t> enabledBranches;
	
	mutable std::shared_mutex dataMtx;
public:
	State(Config const & cfg) : cfg(cfg) {}
public:
	void importFiles();
public:
	memgraph::Graph::NodeId closestNode(sserialize::spatial::GeoPoint const & gp) const;
	std::vector<double> nodeDistances(memgraph::Graph::NodeId const & nId) const;
	std::vector<Distance> branchDistance(BranchId const & branchId, DistanceWeightConfig const & dwc) const;
	void computeBranchAssignments(DistanceWeightConfig const & dwc);
public:
	void writeBranchAssignments(std::ostream & out);
public slots:
	void createBranch(double lat, double lon);
signals:
	void dataChanged();
	void textInfo(QString const & str);
private:
	void emit_textInfo(QString const & str);
};

void registerTypesWithQt();

}//end namespace plz2branch

Q_DECLARE_METATYPE(plz2branch::BranchId)
Q_DECLARE_METATYPE(plz2branch::DistanceWeightConfig);
