////////////////////////////////////////////////////////////////////////
// Class:       Prototyping
// Plugin Type: analyzer (art v2_05_01)
// File:        Prototyping_module.cc
//
// Generated at Wed Feb 14 10:17:20 2018 by Wouter Van De Pontseele using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <fstream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"

#include "TTree.h"

class Prototyping;

class Prototyping : public art::EDAnalyzer
{
public:
  explicit Prototyping(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Prototyping(Prototyping const &) = delete;
  Prototyping(Prototyping &&) = delete;
  Prototyping &operator=(Prototyping const &) = delete;
  Prototyping &operator=(Prototyping &&) = delete;

  // Required functions.
  void clear();
  void clear_PFParticle();
  void analyze(art::Event const &e) override;
  void endSubRun(const art::SubRun &sr);
  void reconfigure(fhicl::ParameterSet const &p) override;

private:
  std::string m_pfp_producer;
  bool m_is_lite;
  bool m_is_data;

  //handle to tpc energy calibration provider
  const lariov::TPCEnergyCalibProvider &energyCalibProvider = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();

  TTree *fPOTTree;
  uint fRun_sr, fSubrun_sr;
  float fPot;

  TTree *fPFParticlesTree;
  uint fRun, fSubrun, fEvent;
  uint fNumPfp;
  int fPdgCode;
  uint fNumDaughters;
  bool fIsPrimary;
  float fVx, fVy, fVz;
  uint fNhits, fNclusters;
  uint fNhitsU, fNhitsV, fNhitsY;
  // track
  bool fTrack_Valid;
  bool fTrack_HasMomentum;
  float fTrack_StartX;
  float fTrack_StartY;
  float fTrack_StartZ;
  float fTrack_EndX;
  float fTrack_EndY;
  float fTrack_EndZ;
  float fTrack_Length;
  float fTrack_StartMomentumX;
  float fTrack_StartMomentumY;
  float fTrack_StartMomentumZ;
  float fTrack_EndMomentumX;
  float fTrack_EndMomentumY;
  float fTrack_EndMomentumZ;
  float fTrack_Theta;
  float fTrack_Phi;
  float fTrack_ZenithAngle;
  float fTrack_AzimuthAngle;
  //shower

  TTree *fClustersTree;
  float fClusterCharge, fClusterWidth, fClusterPosition;
  uint fClusterNhits, fClusterPlane;

  TTree *fHitsTree;
  uint fPlane;
  uint fWire;
  float fCharge;

  TTree *fSpacePointsTree;
  float fx, fy, fz, fChargeU, fChargeV, fChargeY;
};

Prototyping::Prototyping(fhicl::ParameterSet const &p)
    : EDAnalyzer(p) // ,
                    // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;

  this->reconfigure(p);

  fPOTTree = tfs->make<TTree>("pot", "POT Tree");
  fPOTTree->Branch("run", &fRun_sr, "run/i");
  fPOTTree->Branch("subrun", &fSubrun_sr, "subrun/i");
  fPOTTree->Branch("pot", &fPot, "pot/d");

  fPFParticlesTree = tfs->make<TTree>("PFParticles", "PFParticles Tree");
  fPFParticlesTree->Branch("event", &fEvent, "event/i");
  fPFParticlesTree->Branch("run", &fRun, "run/i");
  fPFParticlesTree->Branch("subrun", &fSubrun, "subrun/i");
  fPFParticlesTree->Branch("num_pfp", &fNumPfp, "num_pfp/i");
  fPFParticlesTree->Branch("pdg_code", &fPdgCode, "pdg_code/I");
  fPFParticlesTree->Branch("subrun", &fSubrun, "subrun/i");
  fPFParticlesTree->Branch("num_daughters", &fNumDaughters, "num_daughters/i");
  fPFParticlesTree->Branch("is_primary", &fIsPrimary, "is_primary/i");
  fPFParticlesTree->Branch("n_hits", &fNhits, "n_hits/i");
  fPFParticlesTree->Branch("n_clusters", &fNclusters, "n_clusters/i");
  fPFParticlesTree->Branch("pfp_vx", &fVx, "pfp_vx/F");
  fPFParticlesTree->Branch("pfp_vy", &fVy, "pfp_vy/F");
  fPFParticlesTree->Branch("pfp_vz", &fVz, "pfp_vz/F");
  //track
  fPFParticlesTree->Branch("track_valid", &fTrack_Valid, "track_valid/O");
  fPFParticlesTree->Branch("track_startx", &fTrack_StartX, "track_startx/F");
  fPFParticlesTree->Branch("track_starty", &fTrack_StartY, "track_starty/F");
  fPFParticlesTree->Branch("track_startz", &fTrack_StartZ, "track_startz/F");
  fPFParticlesTree->Branch("track_endx", &fTrack_EndX, "track_endx/F");
  fPFParticlesTree->Branch("track_endy", &fTrack_EndY, "track_endy/F");
  fPFParticlesTree->Branch("track_endz", &fTrack_EndZ, "track_endz/F");
  fPFParticlesTree->Branch("track_length", &fTrack_Length, "track_length/F");
  fPFParticlesTree->Branch("track_hasmomentum", &fTrack_HasMomentum, "track_hasmomentum/O");
  fPFParticlesTree->Branch("track_startmomentumx", &fTrack_StartMomentumX, "track_startmomentumx/F");
  fPFParticlesTree->Branch("track_startmomentumy", &fTrack_StartMomentumY, "track_startmomentumy/F");
  fPFParticlesTree->Branch("track_startmomentumz", &fTrack_StartMomentumZ, "track_startmomentumz/F");
  fPFParticlesTree->Branch("track_endmomentumx", &fTrack_EndMomentumX, "track_endmomentumx/F");
  fPFParticlesTree->Branch("track_endmomentumy", &fTrack_EndMomentumY, "track_endmomentumy/F");
  fPFParticlesTree->Branch("track_endmomentumz", &fTrack_EndMomentumZ, "track_endmomentumz/F");
  fPFParticlesTree->Branch("track_theta", &fTrack_Theta, "track_theta/F");
  fPFParticlesTree->Branch("track_phi", &fTrack_Phi, "track_phi/F");
  fPFParticlesTree->Branch("track_zenith", &fTrack_ZenithAngle, "track_zeninth/F");
  fPFParticlesTree->Branch("track_azimuth", &fTrack_AzimuthAngle, "track_azimuth/F");

  fClustersTree = tfs->make<TTree>("Clusters", "Clusters Tree");
  fClustersTree->Branch("event", &fEvent, "event/i");
  fClustersTree->Branch("run", &fRun, "run/i");
  fClustersTree->Branch("subrun", &fSubrun, "subrun/i");
  fClustersTree->Branch("pdg_code", &fPdgCode, "pdg_code/I");
  fClustersTree->Branch("charge", &fClusterCharge, "charge/F");
  fClustersTree->Branch("width", &fClusterWidth, "width/F");
  fClustersTree->Branch("n_hits", &fClusterNhits, "n_hits/i");
  fClustersTree->Branch("plane", &fClusterPlane, "plane/i");
  fClustersTree->Branch("position", &fClusterPosition, "position/F");

  if (m_is_lite == false)
  {
    fHitsTree = tfs->make<TTree>("Hits", "Hits Tree");
    fHitsTree->Branch("event", &fEvent, "event/i");
    fHitsTree->Branch("run", &fRun, "run/i");
    fHitsTree->Branch("subrun", &fSubrun, "subrun/i");
    fHitsTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
    fHitsTree->Branch("plane", &fPlane, "plane/i");
    fHitsTree->Branch("wire", &fWire, "wire/i");
    fHitsTree->Branch("charge", &fCharge, "charge/d");

    fSpacePointsTree = tfs->make<TTree>("SpacePoints", "Space Points Tree");
    fSpacePointsTree->Branch("event", &fEvent, "event/i");
    fSpacePointsTree->Branch("run", &fRun, "run/i");
    fSpacePointsTree->Branch("subrun", &fSubrun, "subrun/i");
    fSpacePointsTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
    fSpacePointsTree->Branch("x", &fx, "x/F");
    fSpacePointsTree->Branch("y", &fy, "y/F");
    fSpacePointsTree->Branch("z", &fz, "z/F");
    fSpacePointsTree->Branch("chargeU", &fChargeU, "chargeU/F");
    fSpacePointsTree->Branch("chargeV", &fChargeV, "chargeV/F");
    fSpacePointsTree->Branch("chargeY", &fChargeY, "chargeY/F");
  }
}

void Prototyping::reconfigure(fhicl::ParameterSet const &p)
{
  m_pfp_producer = p.get<std::string>("pfp_producer", "pandoraNu");
  m_is_lite = p.get<bool>("is_lite", true);
  m_is_data = p.get<bool>("is_data", false);
}

// Clear once per event
void Prototyping::clear()
{
  fRun_sr = 0;
  fSubrun_sr = 0;
  fPot = 0;

  fRun = 0;
  fSubrun = 0;
  fEvent = 0;

  fNumPfp = 0;
  /*
  // fClustersTree;
  fClusterCharge = -9999;
  fClusterWidth = -9999;
  fClusterPosition = -9999;
  fClusterNhits = 0;
  fClusterPlane = 0;

  // fHitsTree;
  fPlane = 0;
  fWire = 0;
  fCharge = -9999;

  // fSpacePointsTree;
  fx = -9999;
  fy = -9999;
  fz = -9999;
  fChargeU = -9999;
  fChargeV = -9999;
  fChargeY = -9999;
  */
}

// Clear once per PFP
void Prototyping::clear_PFParticle()
{
  // fPFParticlesTree;
  fPdgCode = 0;
  fNumDaughters = 0;
  fIsPrimary = false;
  fVx = -9999;
  fVy = -9999;
  fVz = -9999;
  fNhits = 0;
  fNclusters = 0;
  fNhitsU = 0;
  fNhitsV = 0;
  fNhitsY = 0;
  //// track
  fTrack_Valid = false;
  fTrack_HasMomentum = false;
  fTrack_StartX = -9999;
  fTrack_StartY = -9999;
  fTrack_StartZ = -9999;
  fTrack_EndX = -9999;
  fTrack_EndY = -9999;
  fTrack_EndZ = -9999;
  fTrack_Length = -9999;
  fTrack_StartMomentumX = -9999;
  fTrack_StartMomentumY = -9999;
  fTrack_StartMomentumZ = -9999;
  fTrack_EndMomentumX = -9999;
  fTrack_EndMomentumY = -9999;
  fTrack_EndMomentumZ = -9999;
  fTrack_Theta = -9999;
  fTrack_Phi = -9999;
  fTrack_ZenithAngle = -9999;
  fTrack_AzimuthAngle = -9999;
  //shower
}
