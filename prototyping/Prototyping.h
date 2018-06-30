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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"

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
  void analyze(art::Event const &e) override;
  void reconfigure(fhicl::ParameterSet const &p) override;

private:
  std::string m_pfp_producer;
  bool m_is_lite;

  //handle to tpc energy calibration provider
  const lariov::TPCEnergyCalibProvider &energyCalibProvider = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();

  TTree *fCalibrationTree;
  float xx, yy, zz;
  float corr_u, corr_v, corr_y;

  TTree *fPFParticlesTree;
  int fPdgCode;
  double fvx, fvy, fvz;
  int fNhits, fNclusters;
  int fNhitsU, fNhitsV, fNhitsY;
  int fRun, fSubrun, fEvent;

  TTree *fClustersTree;
  double fClusterCharge, fClusterWidth, fClusterPosition;
  int fClusterNhits, fClusterPlane;

  TTree *fHitsTree;
  int fPlane;
  int fWire;
  double fCharge;

  TTree *fSpacePointsTree;
  double fx, fy, fz, fChargeU, fChargeV, fChargeY;
};

Prototyping::Prototyping(fhicl::ParameterSet const &p)
    : EDAnalyzer(p) // ,
                    // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;

  this->reconfigure(p);

  fCalibrationTree = tfs->make<TTree>("Calibration", "Calibration Tree");
  fCalibrationTree->Branch("x", &xx, "x/F");
  fCalibrationTree->Branch("y", &yy, "y/F");
  fCalibrationTree->Branch("z", &zz, "z/F");
  fCalibrationTree->Branch("Correction U", &corr_u, "Correction U/F");
  fCalibrationTree->Branch("Correction V", &corr_v, "Correction V/F");
  fCalibrationTree->Branch("Correction Y", &corr_y, "Correction Y/F");
  
  fPFParticlesTree = tfs->make<TTree>("PFParticles", "PFParticles Tree");
  fPFParticlesTree->Branch("event", &fEvent, "event/i");
  fPFParticlesTree->Branch("run", &fRun, "run/i");
  fPFParticlesTree->Branch("subrun", &fSubrun, "subrun/i");

  fClustersTree = tfs->make<TTree>("Clusters", "Clusters Tree");
  fClustersTree->Branch("event", &fEvent, "event/i");
  fClustersTree->Branch("run", &fRun, "run/i");
  fClustersTree->Branch("subrun", &fSubrun, "subrun/i");
  fClustersTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fClustersTree->Branch("charge", &fClusterCharge, "charge/d");
  fClustersTree->Branch("width", &fClusterWidth, "width/d");
  fClustersTree->Branch("n_hits", &fClusterNhits, "n_hits/i");
  fClustersTree->Branch("plane", &fClusterPlane, "plane/i");
  fClustersTree->Branch("position", &fClusterPosition, "position/d");
  fClustersTree->Branch("event", &fEvent, "event/i");
  fClustersTree->Branch("run", &fRun, "run/i");
  fClustersTree->Branch("subrun", &fSubrun, "subrun/i");

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
    fSpacePointsTree->Branch("x", &fx, "x/d");
    fSpacePointsTree->Branch("y", &fy, "y/d");
    fSpacePointsTree->Branch("z", &fz, "z/d");
    fSpacePointsTree->Branch("chargeU", &fChargeU, "chargeU/d");
    fSpacePointsTree->Branch("chargeV", &fChargeV, "chargeV/d");
    fSpacePointsTree->Branch("chargeY", &fChargeY, "chargeY/d");
  }
}

void Prototyping::reconfigure(fhicl::ParameterSet const &p)
{
  m_pfp_producer = p.get<std::string>("pfp_producer", "pandoraNu");
  m_is_lite = p.get<bool>("is_lite", false);
}
