#include "dEdx.h"
#include <ffaobjects/EventHeaderv1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertex.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssocv1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackAnalysisUtils.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/recoConsts.h>

#include <CLHEP/Vector/ThreeVector.h>
#include <math.h>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

namespace
{
  template <class T>
  inline T square(const T& t)
  {
    return t * t;
  }
  template <class T>
  inline T r(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track)
  {
    std::vector<TrkrDefs::cluskey> out;
    for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (seed)
      {
        std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
      }
    }

    return out;
  }
}  // namespace

//____________________________________________________________________________..
dEdx::dEdx(const std::string &name, const std::string &file):
 SubsysReco(name),
 _outfilename(file),
 _outfile(nullptr),
 _tree(nullptr)
{
  std::cout << "dEdx::dEdx(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
dEdx::~dEdx()
{
  std::cout << "dEdx::~dEdx() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int dEdx::Init(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  std::cout << "dEdx::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  delete _outfile;
  _outfile = new TFile(_outfilename.c_str(), "RECREATE");

  createBranches();

  m_event=0;

  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsTrkFitter::createNodes");
  }

  PHNodeIterator dstIter(topNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void dEdx::createBranches()
{
  delete _tree;
  _tree = new TTree("tree", "dEdx");
  _tree->Branch("runNumber", &_runNumber);
  _tree->Branch("eventNumber", &_eventNumber);
  _tree->Branch("ntracks", &_ntracks);
  _tree->Branch("nhits", &_nhits);
  _tree->Branch("nmaps", &_nmaps);
  _tree->Branch("nintt", &_nintt);
  _tree->Branch("ntpc", &_ntpc);
  _tree->Branch("nmms", &_nmms);
  _tree->Branch("p", &_p);
  _tree->Branch("pt", &_pt);
  _tree->Branch("theta", &_theta);
  _tree->Branch("eta", &_eta);
  _tree->Branch("phi", &_phi);
  _tree->Branch("charge", &_charge);
  _tree->Branch("quality", &_quality);
  _tree->Branch("crossing", &_crossing);
  _tree->Branch("dedx", &_dedx);
  _tree->Branch("ClusAdcPerLayerThickness_allz", &_ClusAdcPerLayerThickness_allz);
  _tree->Branch("ClusAdcPerLayerThickness_z0", &_ClusAdcPerLayerThickness_z0);
  _tree->Branch("ClusAdcPerLayerThickness_z1", &_ClusAdcPerLayerThickness_z1);
  _tree->Branch("ClusAdcPerLayerThickness_z2", &_ClusAdcPerLayerThickness_z2);
  _tree->Branch("ClusAdcPerLayerThickness_z3", &_ClusAdcPerLayerThickness_z3);
  _tree->Branch("ClusAdcPerLayerThickness_z4", &_ClusAdcPerLayerThickness_z4);
  _tree->Branch("ClusAdcPerLayerThickness_z5", &_ClusAdcPerLayerThickness_z5);
  _tree->Branch("ClusAdcPerLayerThickness_z6", &_ClusAdcPerLayerThickness_z6);
  _tree->Branch("ClusAdcPerLayerThickness_z7", &_ClusAdcPerLayerThickness_z7);
  _tree->Branch("ClusAdcPerLayerThickness_z8", &_ClusAdcPerLayerThickness_z8);
  _tree->Branch("ClusAdcPerLayerThickness_z9", &_ClusAdcPerLayerThickness_z9);
}

void dEdx::ResetTreeVectors()
{
  _runNumber = 0;
  _eventNumber = 0;
  _ntracks = 0;
  _nhits.clear();
  _nmaps.clear();
  _nintt.clear();
  _ntpc.clear();
  _nmms.clear();
  _p.clear();
  _pt.clear();
  _theta.clear();
  _eta.clear();
  _phi.clear();
  _charge.clear();
  _quality.clear();
  _crossing.clear();
  _dedx.clear();
  _ClusAdcPerLayerThickness_allz.clear();
  _ClusAdcPerLayerThickness_z0.clear();
  _ClusAdcPerLayerThickness_z1.clear();
  _ClusAdcPerLayerThickness_z2.clear();
  _ClusAdcPerLayerThickness_z3.clear();
  _ClusAdcPerLayerThickness_z4.clear();
  _ClusAdcPerLayerThickness_z5.clear();
  _ClusAdcPerLayerThickness_z6.clear();
  _ClusAdcPerLayerThickness_z7.clear();
  _ClusAdcPerLayerThickness_z8.clear();
  _ClusAdcPerLayerThickness_z9.clear();
}

//____________________________________________________________________________..
int dEdx::process_event(PHCompositeNode* topNode)
{
  std::cout << "dEdx::process_event event " << m_event << std::endl;

  ResetTreeVectors();

  if(!trackMap)
  {
    trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
    if(!trackMap)
    {
      std::cout << m_trackMapName << " not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (!acts_Geometry)
  {
    acts_Geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    if (!acts_Geometry)
    {
      std::cout << "ActsGeometry not on node tree. Exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(!clusterMap)
  {
    clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if(!clusterMap)
    {
      std::cout << "TRKR_CLUSTER not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(!tpcGeom)
  {
    tpcGeom = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    if(!tpcGeom)
    {
      std::cout << "CYLINDERCELLGEOM_SVTX not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  // global position wrapper
  m_globalPositionWrapper.loadNodes(topNode);

  m_clusterMover.initialize_geometry(tpcGeom);
  m_clusterMover.set_verbosity(0);

  PHNodeIterator nodeIter(topNode);
  PHNode* evtNode = dynamic_cast<PHNode*>(nodeIter.findFirst("EventHeader"));
  if (evtNode)
  {
    EventHeaderv1* evtHeader = findNode::getClass<EventHeaderv1>(topNode, "EventHeader");
    std::cout<<"runNumber = "<<evtHeader->get_RunNumber()<<" , m_evtNumber = "<<evtHeader->get_EvtSequence()<<std::endl;
    _runNumber = evtHeader->get_RunNumber();
    _eventNumber = evtHeader->get_EvtSequence();
  }
  else
  {
    _runNumber = 0;
    _eventNumber = -1;
  }

  SvtxTrack *track = nullptr;

  _ntracks = trackMap->size();
  for (auto &iter : *trackMap)
  {
    track = iter.second;

    if(!track) continue;

    _p.push_back(track->get_p());
    _pt.push_back(track->get_pt());
    _theta.push_back(eta_to_theta(track->get_eta()));
    _eta.push_back(track->get_eta());
    _phi.push_back(track->get_phi());
    _charge.push_back(track->get_charge());
    _quality.push_back(track->get_quality());
    _crossing.push_back(track->get_crossing());

    float fcorr = fabs(sin(eta_to_theta(track->get_eta())));

    int m_nmaps = 0;
    int m_nintt = 0;
    int m_ntpc = 0;
    int m_nmms = 0;

    auto tpcseed = track->get_tpc_seed();
    auto silseed = track->get_silicon_seed();

    if (tpcseed)
    {
      _dedx.push_back(calc_dedx(tpcseed, clusterMap, tpcGeom));
    }

    // get the fully corrected cluster global positions
    std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;
    std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_moved;
    for (const auto& ckey : get_cluster_keys(track))
    {
      auto cluster = clusterMap->findCluster(ckey);

      // Fully correct the cluster positions for the crossing and all distortions
      Acts::Vector3 global = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cluster, track->get_crossing() );

      // add the global positions to a vector to give to the cluster mover
      global_raw.emplace_back(std::make_pair(ckey, global));
    }

    // move the corrected cluster positions back to the original readout surface
    global_moved = m_clusterMover.processTrack(global_raw);

    Acts::Vector3 clusglob_moved(0, 0, 0);
    float adc_allz=0; int nclus_allz=0;
    float adc_z0=0; int nclus_z0=0;
    float adc_z1=0; int nclus_z1=0;
    float adc_z2=0; int nclus_z2=0;
    float adc_z3=0; int nclus_z3=0;
    float adc_z4=0; int nclus_z4=0;
    float adc_z5=0; int nclus_z5=0;
    float adc_z6=0; int nclus_z6=0;
    float adc_z7=0; int nclus_z7=0;
    float adc_z8=0; int nclus_z8=0;
    float adc_z9=0; int nclus_z9=0;
    for (const auto& pair : global_moved)
    {
      auto ckey = pair.first;
      auto cluster = clusterMap->findCluster(ckey);
      clusglob_moved = pair.second;

      auto detid = TrkrDefs::getTrkrId(ckey);
      if (detid == TrkrDefs::TrkrId::mvtxId)
      {
        m_nmaps++;
      }
      else if (detid == TrkrDefs::TrkrId::inttId)
      {
        m_nintt++;
      }
      else if (detid == TrkrDefs::TrkrId::tpcId)
      {
        m_ntpc++;
      }
      else if (detid == TrkrDefs::TrkrId::micromegasId)
      {
        m_nmms++;
      }

      // only counts TPC R2 and R3
      auto layer = TrkrDefs::getLayer(ckey);
      if (layer<23)
      {
        continue;
      }

      float clusgz = clusglob_moved.z();
      adc_allz+=cluster->getAdc() * fcorr;
      nclus_allz++;
      if (clusgz>-100 && clusgz<=-80)
      {
        adc_z0+=cluster->getAdc() * fcorr;
        nclus_z0++;
      }
      else if (clusgz>-80 && clusgz<=-60)
      {
        adc_z1+=cluster->getAdc() * fcorr;
        nclus_z1++;
      }
      else if (clusgz>-60 && clusgz<=-40)
      {
        adc_z2+=cluster->getAdc() * fcorr;
        nclus_z2++;
      }
      else if (clusgz>-40 && clusgz<=-20)
      {
        adc_z3+=cluster->getAdc() * fcorr;
        nclus_z3++;
      }
      else if (clusgz>-20 && clusgz<=0)
      {
        adc_z4+=cluster->getAdc() * fcorr;
        nclus_z4++;
      }
      else if (clusgz>0 && clusgz<=20)
      {
        adc_z5+=cluster->getAdc() * fcorr;
        nclus_z5++;
      }
      else if (clusgz>20 && clusgz<=40)
      {
        adc_z6+=cluster->getAdc() * fcorr;
        nclus_z6++;
      }
      else if (clusgz>40 && clusgz<=60)
      {
        adc_z7+=cluster->getAdc() * fcorr;
        nclus_z7++;
      }
      else if (clusgz>60 && clusgz<=80)
      {
        adc_z8+=cluster->getAdc() * fcorr;
        nclus_z8++;
      }
      else if (clusgz>80 && clusgz<=100)
      {
        adc_z9+=cluster->getAdc() * fcorr;
        nclus_z9++;
      }
    }

    int m_nhits = m_nmaps + m_nintt + m_ntpc + m_nmms;

    _nhits.push_back(m_nhits);
    _nmaps.push_back(m_nmaps);
    _nintt.push_back(m_nintt);
    _ntpc.push_back(m_ntpc);
    _nmms.push_back(m_nmms);

    adc_allz /= nclus_allz;
    adc_z0 /= nclus_z0;
    adc_z1 /= nclus_z1;
    adc_z2 /= nclus_z2;
    adc_z3 /= nclus_z3;
    adc_z4 /= nclus_z4;
    adc_z5 /= nclus_z5;
    adc_z6 /= nclus_z6;
    adc_z7 /= nclus_z7;
    adc_z8 /= nclus_z8;
    adc_z9 /= nclus_z9;

    _ClusAdcPerLayerThickness_allz.push_back(adc_allz);
    _ClusAdcPerLayerThickness_z0.push_back(adc_z0);
    _ClusAdcPerLayerThickness_z1.push_back(adc_z1);
    _ClusAdcPerLayerThickness_z2.push_back(adc_z2);
    _ClusAdcPerLayerThickness_z3.push_back(adc_z3);
    _ClusAdcPerLayerThickness_z4.push_back(adc_z4);
    _ClusAdcPerLayerThickness_z5.push_back(adc_z5);
    _ClusAdcPerLayerThickness_z6.push_back(adc_z6);
    _ClusAdcPerLayerThickness_z7.push_back(adc_z7);
    _ClusAdcPerLayerThickness_z8.push_back(adc_z8);
    _ClusAdcPerLayerThickness_z9.push_back(adc_z9);
  }

  _tree->Fill();

  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

float dEdx::calc_dedx(TrackSeed* tpcseed, TrkrClusterContainer* clustermap, PHG4TpcCylinderGeomContainer* tpcgeom)
{
  std::vector<TrkrDefs::cluskey> clusterKeys;
  clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
                     tpcseed->end_cluster_keys());

  std::vector<float> dedxlist;
  for (unsigned long cluster_key : clusterKeys)
  {
    auto detid = TrkrDefs::getTrkrId(cluster_key);
    if (detid != TrkrDefs::TrkrId::tpcId)
    {
      continue;  // the micromegas clusters are added to the TPC seeds
    }
    unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
    TrkrCluster* cluster = clustermap->findCluster(cluster_key);
    float adc = cluster->getAdc();
    PHG4TpcCylinderGeom* GeoLayer_local = tpcgeom->GetLayerCellGeom(layer_local);
    float thick = GeoLayer_local->get_thickness();
    float r = GeoLayer_local->get_radius();
    float alpha = (r * r) / (2 * r * TMath::Abs(1.0 / tpcseed->get_qOverR()));
    float beta = atan(tpcseed->get_slope());
    float alphacorr = cos(alpha);
    if (alphacorr < 0 || alphacorr > 4)
    {
      alphacorr = 4;
    }
    float betacorr = cos(beta);
    if (betacorr < 0 || betacorr > 4)
    {
      betacorr = 4;
    }
    adc /= thick;
    adc *= alphacorr;
    adc *= betacorr;
    dedxlist.push_back(adc);
    sort(dedxlist.begin(), dedxlist.end());
  }
  int trunc_min = 0;
  if (dedxlist.size() < 1)
  {
    return std::numeric_limits<float>::quiet_NaN();
  }
  int trunc_max = (int) dedxlist.size() * 0.7;
  float sumdedx = 0;
  int ndedx = 0;
  for (int j = trunc_min; j <= trunc_max; j++)
  {
    sumdedx += dedxlist.at(j);
    ndedx++;
  }
  sumdedx /= ndedx;
  return sumdedx;
}

float dEdx::eta_to_theta(float eta)
{
  return 2*atan(exp(-eta));
}

//____________________________________________________________________________..
int dEdx::End(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  _outfile->cd();
  _outfile->Write();
  _outfile->Close();
  std::cout << "dEdx::End(PHCompositeNode *topNode) Endding" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
