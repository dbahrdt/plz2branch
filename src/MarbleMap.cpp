#include "MarbleMap.h"
#include <marble/GeoPainter.h>
#include <marble/MarbleWidgetPopupMenu.h>
#include <marble/MarbleWidgetInputHandler.h>
#include <marble/GeoDataLatLonBox.h>
#include <QAction>

namespace plz2branch {

MarbleMap::MyBaseLayer::MyBaseLayer(const QStringList& renderPos, qreal zVal, std::shared_ptr<State> state) :
m_zValue(zVal),
m_renderPosition(renderPos),
m_state(std::move(state))
{}

QStringList MarbleMap::MyBaseLayer::renderPosition() const {
	return m_renderPosition;
}

qreal MarbleMap::MyBaseLayer::zValue() const {
    return m_zValue;
}

MarbleMap::MyNodesLayer::MyNodesLayer(const QStringList& renderPos, qreal zVal, std::shared_ptr<State> state):
MyLockableBaseLayer(renderPos, zVal, std::move(state))
{}

bool MarbleMap::MyNodesLayer::render(Marble::GeoPainter* painter, Marble::ViewportParams*, const QString&, Marble::GeoSceneLayer*) {	

	
	std::shared_lock<std::shared_mutex> lck(state().dataMtx);
	for(uint32_t branchId : state().enabledBranches) {
		Branch const & bri = state().branches.at(branchId);
		painter->setPen(bri.color);
		painter->setBrush(QBrush(bri.color, Qt::BrushStyle::SolidPattern));
		
		Marble::GeoDataCoordinates gp(bri.coord.lon(), bri.coord.lat(), 0.0, Marble::GeoDataCoordinates::Degree);
		painter->drawEllipse(gp, 20, 20);
		painter->drawText(gp, QString::number(branchId));
	}
	return true;
}

//BEGIN MyGeometryLayer

MarbleMap::MyGeometryLayer::MyGeometryLayer(const QStringList & renderPos, qreal zVal, std::shared_ptr<State> state) :
MyLockableBaseLayer(renderPos, zVal, state)
{}

MarbleMap::MyGeometryLayer::~MyGeometryLayer() {}

bool MarbleMap::MyGeometryLayer::render(Marble::GeoPainter* painter, Marble::ViewportParams* viewport, const QString&
 renderPos, Marble::GeoSceneLayer* layer) {
	std::shared_lock<std::shared_mutex> lck(state().dataMtx);
	for(uint32_t branchId : state().enabledBranches) {
		Branch const & bri = state().branches.at(branchId);
		painter->setPen(bri.color);
		painter->setBrush(QBrush(bri.color, Qt::Dense7Pattern));
		for(const auto & [rId, rDist] : bri.assignedRegions) {
			sserialize::spatial::GeoRect const & bbox = state().regionInfo.at(rId.value()).bbox;
			Marble::GeoDataLatLonBox lb(bbox.maxLat(), bbox.minLat(), bbox.maxLon(), bbox.minLon(), Marble::GeoDataCoordinates::Degree);
			painter->drawRect(lb.center(), lb.width(Marble::GeoDataCoordinates::Degree), lb.height(Marble::GeoDataCoordinates::Degree), true);
			
			auto msg = QString("plz=%1, br=%2, dt=%3, ds=%4").arg(QString::number(state().regionInfo.at(rId.value()).plz), QString::number(branchId), QString::number(rDist.time), QString::number(rDist.spatial));
			painter->drawText(lb.center(), msg);
		}
	}
	return true;
}

//END MyGeometryLayer

MarbleMap::MarbleMap(QWidget * parent, std::shared_ptr<State> state) :
MarbleWidget(parent),
m_state(state)
{
	m_nodesLayer = new MyNodesLayer({"HOVERS_ABOVE_SURFACE"}, 0.0, state);
	m_plzLayer = new MyGeometryLayer({"HOVERS_ABOVE_SURFACE"}, 0.0, state);
	
	QAction * branchCreateME = new QAction("Create branch", this);
	
	addLayer(m_plzLayer);
	addLayer(m_nodesLayer);
	
	popupMenu()->addAction(Qt::MouseButton::RightButton, branchCreateME);
	
	//get mouse clicks
	connect(this->inputHandler(), SIGNAL(rmbRequest(int,int)), this, SLOT(rmbRequested(int,int)));
	
	connect(branchCreateME, SIGNAL(triggered(bool)), this, SLOT(slot_branchCreated()));
	
	centerOn(Marble::GeoDataCoordinates(9.17, 48.78, 100.0, Marble::GeoDataCoordinates::Degree));
}

MarbleMap::~MarbleMap() {
	removeLayer(m_nodesLayer);
	removeLayer(m_plzLayer);
	delete m_nodesLayer;
	delete m_plzLayer;
}

void MarbleMap::dataChanged() {
	this->update();
}

void MarbleMap::zoomToBranch(BranchId brId) {
	Marble::GeoDataCoordinates geo;
	{
		std::shared_lock<std::shared_mutex> lck(m_state->dataMtx);
		auto const & bri = m_state->branches.at(brId.value);
		geo = Marble::GeoDataCoordinates(bri.coord.lon(), bri.coord.lat(), 100.0, Marble::GeoDataCoordinates::Degree);
	}
	centerOn(geo);
	setZoom(100000);
}

void MarbleMap::rmbRequested(int x, int y) {
	this->geoCoordinates(x, y, m_lastMouseClickLon, m_lastMouseClickLat, Marble::GeoDataCoordinates::Degree);
}

void MarbleMap::slot_branchCreated() {
	branchCreated(m_lastMouseClickLat, m_lastMouseClickLon);
}

}//end namespace simpleroute
