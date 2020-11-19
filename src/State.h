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
	std::string plzFile;
};

struct RegionInfo {
	sserialize::spatial::GeoRect bbox;
	uint32_t plz{std::numeric_limits<uint32_t>::max()};
	uint32_t inhabitants{1};
	inline bool valid() const { return plz != std::numeric_limits<uint32_t>::max(); }
};

class RegionId {
public:
	RegionId() {}
	RegionId(uint32_t v) : m_value(v) {}
	RegionId(RegionId const&) = default;
	RegionId & operator=(RegionId const &) = default;
	~RegionId() {}
public:
	uint32_t value() const {
		assert(valid());
		return m_value;
	}
	inline bool valid() const { return m_value != std::numeric_limits<uint32_t>::max(); }
private:
	uint32_t m_value{std::numeric_limits<uint32_t>::max()};
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
	enum NodeWeightModel {
		Min,
		Max,
		Median,
		Mean
	};
	enum class RWM {
		Equal,
		Inhabitants //weight by number of inhabitants in the given plz
	};
	enum class BWM {
		Equal,
		Employees //weight by number of employees
	};
	enum class BAM { //Branch assignment
		ByDistance,
		Greedy,
	};
	DistanceWeightConfig() {}
	DistanceWeightConfig(NodeSelection ns, NodeWeightModel nwm, RWM rwm, BWM bwm, BAM bam) :
	nodeSelection(ns),
	nodeWeightModel(nwm),
	regionWeightModel(rwm),
	branchWeightModel(bwm),
	branchAssignmentModel(bam)
	{}
	DistanceWeightConfig(DistanceWeightConfig const&) = default;
	DistanceWeightConfig & operator=(DistanceWeightConfig const&) = default;
	NodeSelection nodeSelection{All};
	NodeWeightModel nodeWeightModel{Min};
	RWM regionWeightModel{RWM::Inhabitants};
	BWM branchWeightModel{BWM::Employees};
	BAM branchAssignmentModel{BAM::ByDistance};
};

std::string to_string(DistanceWeightConfig::NodeSelection v);
std::string to_string(DistanceWeightConfig::NodeWeightModel v);

void from_string(std::string const & src, DistanceWeightConfig::NodeWeightModel & dest);
void from_string(std::string const & src, DistanceWeightConfig::NodeSelection & dest);

struct Distance {
	double value{std::numeric_limits<double>::max()};
};

struct Branch {
	QString name;
	sserialize::spatial::GeoPoint coord;
	uint32_t employees; //number of employees
	std::vector<std::pair<RegionId, Distance>> assignedRegions;
	
	//Stuff calculated
	memgraph::Graph::NodeId nodeId;
	QColor color{Qt::green};
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
	void createBranch(double lat, double lon, QString name);
	void createBranch(double lat, double lon);
	void createBranches(std::istream & data);
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
