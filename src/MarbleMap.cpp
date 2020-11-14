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
	painter->setBrush(Qt::BrushStyle::SolidPattern);
	
	std::shared_lock<std::shared_mutex> lck(state().dataMtx);
	for(uint32_t branchId : state().enabledBranches) {
	
		Marble::GeoDataCoordinates gp;
		{
			std::shared_lock<std::shared_mutex> lck(state().dataMtx);
			Branch const & bri = state().branches.at(branchId);
			painter->setPen(bri.color);
			gp = Marble::GeoDataCoordinates(bri.coord.lon(), bri.coord.lat(), 0.0, Marble::GeoDataCoordinates::Degree);
		}
		painter->drawEllipse(gp, 10, 10);
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
	
	QPen pen(( QColor(Qt::blue) ));
	QColor fillColor(0, 0, 255, 50);
	QBrush brush(fillColor);
	
	painter->setBrush(Qt::BrushStyle::SolidPattern);
	
	std::shared_lock<std::shared_mutex> lck(state().dataMtx);
	for(uint32_t branchId : state().enabledBranches) {
		std::shared_lock<std::shared_mutex> lck(state().dataMtx);
		Branch const & bri = state().branches.at(branchId);
		painter->setPen(bri.color);
		for(auto rId : bri.assignedRegions) {
			sserialize::spatial::GeoRect bbox = state().regionInfo.at(rId.value).shape->boundary();
			Marble::GeoDataLatLonBox lb(bbox.maxLat(), bbox.minLat(), bbox.maxLon(), bbox.minLon(), Marble::GeoDataCoordinates::Degree);
			painter->drawRect(lb.center(), lb.width(Marble::GeoDataCoordinates::Degree), lb.height(Marble::GeoDataCoordinates::Degree), true);
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
	
	QAction * branchCreateME = new QAction("Create branch", this);
	
	addLayer(m_nodesLayer);
	
	popupMenu()->addAction(Qt::MouseButton::RightButton, branchCreateME);
	
	//get mouse clicks
	connect(this->inputHandler(), SIGNAL(rmbRequest(int,int)), this, SLOT(rmbRequested(int,int)));
	
	connect(branchCreateME, SIGNAL(triggered(bool)), this, SLOT(slot_branchCreated()));
	
	centerOn(Marble::GeoDataCoordinates(9.17, 48.78, 100.0, Marble::GeoDataCoordinates::Degree));
}

MarbleMap::~MarbleMap() {
	removeLayer(m_nodesLayer);
	delete m_nodesLayer;
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
