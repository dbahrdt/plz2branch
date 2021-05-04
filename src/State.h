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
	enum class NS {
		All,
		OnlyCrossRoads
	};
	enum class NWM {
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
		Evolutionary,
		Optimal
	};
	DistanceWeightConfig() {}
	DistanceWeightConfig(NS ns, NWM nwm, RWM rwm = RWM::Equal, BWM bwm = BWM::Equal, BAM bam = BAM::ByDistance) :
	nodeSelection(ns),
	nodeWeightModel(nwm),
	regionWeightModel(rwm),
	branchWeightModel(bwm),
	branchAssignmentModel(bam)
	{}
	DistanceWeightConfig(DistanceWeightConfig const&) = default;
	DistanceWeightConfig & operator=(DistanceWeightConfig const&) = default;
	NS nodeSelection{NS::All};
	NWM nodeWeightModel{NWM::Min};
	RWM regionWeightModel{RWM::Inhabitants};
	BWM branchWeightModel{BWM::Employees};
	BAM branchAssignmentModel{BAM::ByDistance};
};

struct Distance {
	double time{std::numeric_limits<double>::max()};
	double spatial{std::numeric_limits<double>::max()};
};

struct Branch {
	QString name;
	sserialize::spatial::GeoPoint coord;
	uint32_t employees{1}; //number of employees
	uint32_t maxRegions{std::numeric_limits<uint32_t>::max()}; //maximum number of regions
	std::vector<std::pair<RegionId, Distance>> assignedRegions;
	
	//Stuff calculated
	memgraph::Graph::NodeId nodeId;
	QColor color{Qt::green};
	double cost{0};
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
	
	double totalCost{0};
	
	mutable std::shared_mutex dataMtx;
public:
	State(Config const & cfg) : cfg(cfg) {}
public:
	void importFiles();
public:
	memgraph::Graph::NodeId closestNode(sserialize::spatial::GeoPoint const & gp) const;
	std::vector<double> nodeDistances(memgraph::Graph::NodeId const & nId, std::function<double(memgraph::Graph::Edge const &)>) const;
	std::vector<Distance> branchDistance(BranchId const & branchId, DistanceWeightConfig const & dwc) const;
	void computeBranchAssignments(DistanceWeightConfig const & dwc);
public:
	void writeBranchAssignments(std::ostream & out);
	void exportBranchAssignments(std::ostream & out);
public slots:
	void createBranch(double lat, double lon, QString name, uint32_t employees);
	void createBranch(double lat, double lon);
	void createBranches(std::istream & data);
	void clear();
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
Q_DECLARE_METATYPE(plz2branch::DistanceWeightConfig::NS);
Q_DECLARE_METATYPE(plz2branch::DistanceWeightConfig::NWM);
Q_DECLARE_METATYPE(plz2branch::DistanceWeightConfig::RWM);
Q_DECLARE_METATYPE(plz2branch::DistanceWeightConfig::BWM);
Q_DECLARE_METATYPE(plz2branch::DistanceWeightConfig::BAM);
