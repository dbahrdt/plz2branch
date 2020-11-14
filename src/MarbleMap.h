#pragma once

#include <marble/MarbleWidget.h>
#include <marble/LayerInterface.h>
#include <marble/GeoDataLineString.h>

#include <unordered_set>
#include <shared_mutex>

#include "State.h"

namespace plz2branch {

class MarbleMap : public Marble::MarbleWidget {
	Q_OBJECT
public:
	MarbleMap(QWidget * parent, std::shared_ptr<State> state);
	virtual ~MarbleMap();
public slots:
	void dataChanged();
	void zoomToBranch(BranchId brId);
signals:
	void branchCreated(double lat, double lon);
private:
	class MyBaseLayer: public Marble::LayerInterface {
	private:
		qreal m_zValue;
		QStringList m_renderPosition;
		std::shared_ptr<State> m_state;
	protected:
		State const & state() const { return *m_state; }
		State & state() { return *m_state; }
	public:
		MyBaseLayer(const QStringList & renderPos, qreal zVal, std::shared_ptr<State> state);
		~MyBaseLayer() override {}
		QStringList renderPosition() const override;
		qreal zValue() const override;
	};
	
	class MyLockableBaseLayer: public MyBaseLayer {
	private:
		std::shared_mutex m_lock;
	protected:
		inline std::shared_mutex & lock() { return m_lock;}
	public:
		MyLockableBaseLayer(const QStringList & renderPos, qreal zVal, std::shared_ptr<State> state) :
		MyBaseLayer(renderPos, zVal, state) {}
		~MyLockableBaseLayer() override {}
	};
	
	class MyNodesLayer: public MyLockableBaseLayer {
	public:
		MyNodesLayer(const QStringList & renderPos, qreal zVal, std::shared_ptr<State> state);
		~MyNodesLayer() override {}
		bool render(Marble::GeoPainter *painter, Marble::ViewportParams * viewport, const QString & renderPos, Marble::GeoSceneLayer * layer)  override;
	};
	
	class MyGeometryLayer: public MyLockableBaseLayer {
	public:
		MyGeometryLayer(const QStringList & renderPos, qreal zVal, std::shared_ptr<State> state);
		virtual ~MyGeometryLayer();
	public:
		bool render(Marble::GeoPainter *painter, Marble::ViewportParams * viewport, const QString & renderPos, Marble::GeoSceneLayer * layer) override;
	};
	
private slots:
	void rmbRequested(int x, int y);
	void slot_branchCreated();
private:
	std::shared_ptr<State> m_state;
	MyNodesLayer * m_nodesLayer;
	MyGeometryLayer * m_plzLayer;
	double m_lastMouseClickLat;
	double m_lastMouseClickLon;
};


}//end namespace plz2branch
