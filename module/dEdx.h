#ifndef DEDX_H
#define DEDX_H

#include <fun4all/SubsysReco.h>

#include <tpc/TpcClusterMover.h>
#include <tpc/TpcGlobalPositionWrapper.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertex.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterContainer.h>

#include <string>
#include <vector>

class PHCompositeNode;
class PHNode;
class TH1;
class TH2;
class TFile;
class TTree;
class TrkrCluster;
class PHCompositeNode;
class ActsGeometry;
class SvtxTrack;
class TrackSeed;
class TrkrClusterContainer;
class TrkrHitSetContainer;
class PHG4TpcCylinderGeomContainer;
class PHG4CylinderGeomContainer;
class TpcDistortionCorrectionContainer;

class dEdx : public SubsysReco
{
 public:

  dEdx(const std::string &name = "dEdx", const std::string &file = "output.root");

  ~dEdx() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  float eta_to_theta(float eta);
  float calc_dedx(TrackSeed* tpcseed, TrkrClusterContainer* clustermap, PHG4TpcCylinderGeomContainer* tpcgeom);
  void ResetTreeVectors();
  void createBranches();

  std::string GetTrackMapName() {return m_trackMapName;}
  void SetTrackMapName(std::string name) {m_trackMapName = name;}

 private:
  int m_event = 0;
  std::string _outfilename;
  TFile *_outfile = nullptr;
  TTree *_tree = nullptr;

  TpcClusterMover m_clusterMover;
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  int _runNumber = 0;
  int _eventNumber = 0;
  int _ntracks = 0;
  std::vector<int> _nhits;
  std::vector<int> _nmaps;
  std::vector<int> _nintt;
  std::vector<int> _ntpc;
  std::vector<int> _nmms;

  std::vector<float> _p;
  std::vector<float> _pt;
  std::vector<float> _theta;
  std::vector<float> _eta;
  std::vector<float> _phi;
  std::vector<int> _charge;
  std::vector<float> _quality;
  std::vector<int> _crossing;
  std::vector<int> _vtxvalid;
  std::vector<float> _vx;
  std::vector<float> _vy;
  std::vector<float> _vz;
  std::vector<float> _dedx;
  std::vector<float> _ClusAdcPerLayerThickness_allz;
  std::vector<float> _ClusAdcPerLayerThickness_z0;
  std::vector<float> _ClusAdcPerLayerThickness_z1;
  std::vector<float> _ClusAdcPerLayerThickness_z2;
  std::vector<float> _ClusAdcPerLayerThickness_z3;
  std::vector<float> _ClusAdcPerLayerThickness_z4;
  std::vector<float> _ClusAdcPerLayerThickness_z5;
  std::vector<float> _ClusAdcPerLayerThickness_z6;
  std::vector<float> _ClusAdcPerLayerThickness_z7;
  std::vector<float> _ClusAdcPerLayerThickness_z8;
  std::vector<float> _ClusAdcPerLayerThickness_z9;

  SvtxTrackMap* trackMap = nullptr;
  SvtxVertexMap* vertexMap = nullptr;
  ActsGeometry* acts_Geometry = nullptr;
  TrkrClusterContainer* clusterMap = nullptr;
  PHG4TpcCylinderGeomContainer* tpcGeom = nullptr;

  std::string m_trackMapName = "SvtxTrackMap";
};

#endif
// DEDX_H
