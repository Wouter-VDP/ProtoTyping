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
#include <set>

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
#include "lardataobj/RecoBase/OpFlash.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

//#include "larsim/MCCheater/BackTracker.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "TTree.h"
#include "TVector3.h"

#include "GeometryHelper.h"

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
  void clear_MCParticle();
  void clear_PFParticle();
  void clear_Cluster();
  void clear_Flashes();
  // see cc file
  void analyze(art::Event const &e) override;
  void endSubRun(const art::SubRun &sr);
  void reconfigure(fhicl::ParameterSet const &p) override;

  void fill_flash(art::ValidHandle<std::vector<recob::OpFlash>> const &simple_cosmic_handle, uint number, TTree *tree);

  // Larpandora way of backtracking
  art::Ptr<simb::MCTruth> TrackIDToMCTruth(art::Event const &e, std::string _geant_producer, int geant_track_id);

private:
  // FCL parameters
  std::string m_pfp_producer;
  std::string m_cosmic_opflash_producer;
  std::string m_cosmic_simpleflash_producer;
  std::string m_cosmic_ophit_producer;
  std::string m_beam_opflash_producer;
  std::string m_beam_simpleflash_producer;
  std::string m_beam_ophit_producer;

  bool m_is_lite;
  bool m_is_data;
  bool m_verb;

  // Other firvate fields
  GeometryHelper geoHelper;
  //art::ServiceHandle<cheat::BackTracker> bt;

  std::set<std::string> string_process; // This variable counts the different processes invloved.

  const lariov::TPCEnergyCalibProvider &energyCalibProvider = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();
  art::ServiceHandle<geo::Geometry> geo;
  std::map<std::string, uint> map_process =
      {
          {"0", 0},
          {"CoulombScat", 1},
          {"Decay", 2},
          {"annihil", 3},
          {"compt", 4},
          {"conv", 5},
          {"dInelastic", 6},
          {"eBrem", 7},
          {"eIoni", 8},
          {"hBertiniCaptureAtRest", 9},
          {"hIoni", 10},
          {"hadElastic", 11},
          {"muBrems", 12},
          {"muIoni", 13},
          {"muMinusCaptureAtRest", 14},
          {"muPairProd", 15},
          {"muonNuclear", 16},
          {"nCapture", 17},
          {"neutronInelastic", 18},
          {"phot", 19},
          {"photonNuclear", 20},
          {"pi+Inelastic", 21},
          {"pi-Inelastic", 22},
          {"primary", 23},
          {"protonInelastic", 24},
      };

  // Fields for in the tree!
  TTree *fPOTTree;
  uint fRun_sr, fSubrun_sr;
  float fPot;
  float fDatasetPrescaleFactor; // Prescaling factor for data.
  uint fNevents = 0;

  TTree *fEventTree;
  uint fRun, fSubrun, fEvent;
  uint fNumSimpleBeamFlashes;
  uint fNumOpBeamFlashes;
  uint fNumSimpleCosmicFlashes;
  uint fNumOpCosmicFlashes;
  uint fNumMcp;
  uint fNumMcp_saved; // criteria to keep only relevant particles.
  uint fNumPfp;

  TTree *fMCParticlesTree;
  float fMc_E;
  int fMc_PdgCode;
  bool fMc_StartInside;
  bool fMc_EndInside;
  bool fMc_PartInside; // This means that the track is crossing, starts/ends inside, is completely inside.
  bool fMc_kBeamNeutrino;  // Uses the backtracker
  float fMc_StartX;
  float fMc_StartY;
  float fMc_StartZ;
  float fMc_StartX_tpc; // if the track starts outside but crosses, this will be stored here.
  float fMc_StartY_tpc;
  float fMc_StartZ_tpc;
  float fMc_StartX_sce; // spacecharge corrected version of the tpc edge or the inside start end point.
  float fMc_StartY_sce;
  float fMc_StartZ_sce;
  float fMc_EndX;
  float fMc_EndY;
  float fMc_EndZ;
  float fMc_EndX_tpc;
  float fMc_EndY_tpc;
  float fMc_EndZ_tpc;
  float fMc_EndX_sce;
  float fMc_EndY_sce;
  float fMc_EndZ_sce;
  float fMc_Length;
  float fMc_LengthTPC;
  float fMc_StartMomentumX;
  float fMc_StartMomentumY;
  float fMc_StartMomentumZ;
  uint fMc_Process; // std::string
  int fMc_StatusCode;
  float fMc_Time;

  TTree *fSimpleCosmicFlashesTree;
  TTree *fOpCosmicFlashesTree;
  TTree *fSimpleBeamFlashesTree;
  TTree *fOpBeamFlashesTree;
  float fFlash_Time;
  uint fFlash_TotalPE;
  uint fFlash_num10percentPMT; // The number of PMT that is responsible for more than 10% of the total flash.
  float fFlash_Z;
  float fFlash_sigmaZ;
  float fFlash_Y;
  float fFlash_sigmaY;
  float fFlash_Width;
  float fFlash_AbsTime;

  TTree *fPFParticlesTree;
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
};

Prototyping::Prototyping(fhicl::ParameterSet const &p)
    : EDAnalyzer(p) // ,
                    // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;

  this->reconfigure(p);

  //// Check if things are set up properly:
  std::cout << std::endl;
  std::cout << "[Prototyping constructor] Checking set-up" << std::endl;
  //// Check if spacecharge correction is working
  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  std::vector<double> sce_start = SCE->GetPosOffsets(25, 110, 250); // Seemingly this results in 0
  auto scecorr = SCE->GetPosOffsets(25, 110, 250);
  //double g4Ticks = detClocks->TPCG4Time2Tick(mct.GetNeutrino().Lepton().T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
  //double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr[0];
  double yOffset = scecorr[1];
  double zOffset = scecorr[2];
  std::cout << "[Prototyping constructor] Spacecharge correction at test point " << yOffset << ", " << yOffset << ", " << zOffset << std::endl;
  std::cout << std::endl;
  //// End set up tests

  //// Tree for every subrun
  fPOTTree = tfs->make<TTree>("pot", "POT Tree");
  fPOTTree->Branch("run", &fRun_sr, "run/i");
  fPOTTree->Branch("subrun", &fSubrun_sr, "subrun/i");
  fPOTTree->Branch("pot", &fPot, "pot/d");
  fPOTTree->Branch("n_events", &fNevents, "n_events/i");

  //// Tree for every event
  fEventTree = tfs->make<TTree>("Event", "Event Tree");
  fEventTree->Branch("event", &fEvent, "event/i");
  fEventTree->Branch("run", &fRun_sr, "run/i");
  fEventTree->Branch("subrun", &fSubrun_sr, "subrun/i");
  fEventTree->Branch("pot", &fPot, "pot/d");
  fEventTree->Branch("n_events", &fNevents, "n_events/i");
  fEventTree->Branch("dataset_prescale_factor", &fDatasetPrescaleFactor, "dataset_prescale_factor/F");

  fEventTree->Branch("num_simplebeamflashes", &fNumSimpleBeamFlashes, "num_simplebeamflashes/i");
  fEventTree->Branch("num_opbeamflashes", &fNumOpBeamFlashes, "num_opbeamflashes/i");
  fEventTree->Branch("num_simplecosmicflashes", &fNumSimpleCosmicFlashes, "num_simplecosmicflashes/i");
  fEventTree->Branch("num_opcosmicflashes", &fNumOpCosmicFlashes, "num_opcosmicflashes/i");

  fEventTree->Branch("num_pfp", &fNumPfp, "num_pfp/i");
  if (!m_is_data)
  {
    fEventTree->Branch("num_mcp", &fNumMcp, "num_mcp/i");
    fEventTree->Branch("num_mcp_saved", &fNumMcp_saved, "num_mcp_saved/i");
  }

  //// Tree for MC particles
  if (!m_is_data)
  {
    fMCParticlesTree = tfs->make<TTree>("MCParticles", "MCParticles Tree");
    fMCParticlesTree->Branch("event", &fEvent, "event/i");
    fMCParticlesTree->Branch("run", &fRun, "run/i");
    fMCParticlesTree->Branch("subrun", &fSubrun, "subrun/i");
    fMCParticlesTree->Branch("num_mcp", &fNumMcp, "num_mcp/i");
    fMCParticlesTree->Branch("num_mcp_saved", &fNumMcp_saved, "num_mcp_saved/i");
    fMCParticlesTree->Branch("mc_energy", &fMc_E, "mc_energy/F");
    fMCParticlesTree->Branch("mc_pdg_code", &fMc_PdgCode, "mc_pdg_code/I");
    fMCParticlesTree->Branch("mc_status_code", &fMc_StatusCode, "mc_status_code/I");
    fMCParticlesTree->Branch("mc_process", &fMc_Process, "mc_pocess/i");
    fMCParticlesTree->Branch("mc_start_inside", &fMc_StartInside, "mc_start_inside/O");
    fMCParticlesTree->Branch("mc_end_inside", &fMc_EndInside, "mc_end_inside/O");
    fMCParticlesTree->Branch("fMc_part_inside", &fMc_PartInside, "mc_part_inside/O");
    fMCParticlesTree->Branch("fMc_kBeamNeutrino", &fMc_kBeamNeutrino, "mcc_kBeamNeutrino/O");
    fMCParticlesTree->Branch("mc_time", &fMc_Time, "mc_time/F");
    fMCParticlesTree->Branch("mc_startx", &fMc_StartX, "mc_startx/F");
    fMCParticlesTree->Branch("mc_starty", &fMc_StartY, "mc_starty/F");
    fMCParticlesTree->Branch("mc_startz", &fMc_StartZ, "mc_startz/F");
    fMCParticlesTree->Branch("mc_startx_tpc", &fMc_StartX_tpc, "mc_startx_tpc/F");
    fMCParticlesTree->Branch("mc_starty_tpc", &fMc_StartY_tpc, "mc_starty_tpc/F");
    fMCParticlesTree->Branch("mc_startz_tpc", &fMc_StartZ_tpc, "mc_startz_tpc/F");
    fMCParticlesTree->Branch("mc_startx_sce", &fMc_StartX_sce, "mc_startx_sce/F");
    fMCParticlesTree->Branch("mc_starty_sce", &fMc_StartY_sce, "mc_starty_sce/F");
    fMCParticlesTree->Branch("mc_startz_sce", &fMc_StartZ_sce, "mc_startz_sce/F");
    fMCParticlesTree->Branch("mc_endx", &fMc_EndX, "mc_endx/F");
    fMCParticlesTree->Branch("mc_endy", &fMc_EndY, "mc_endy/F");
    fMCParticlesTree->Branch("mc_endz", &fMc_EndZ, "mc_endz/F");
    fMCParticlesTree->Branch("mc_endx_tpc", &fMc_EndX_tpc, "mc_endx_tpc/F");
    fMCParticlesTree->Branch("mc_endy_tpc", &fMc_EndY_tpc, "mc_endy_tpc/F");
    fMCParticlesTree->Branch("mc_endz_tpc", &fMc_EndZ_tpc, "mc_endz_tpc/F");
    fMCParticlesTree->Branch("mc_endx_sce", &fMc_EndX_sce, "mc_endx_sce/F");
    fMCParticlesTree->Branch("mc_endy_sce", &fMc_EndY_sce, "mc_endy_sce/F");
    fMCParticlesTree->Branch("mc_endz_sce", &fMc_EndZ_sce, "mc_endz_sce/F");
    fMCParticlesTree->Branch("mc_startmomentumx", &fMc_StartMomentumX, "mc_startmomentumx/F");
    fMCParticlesTree->Branch("mc_startmomentumy", &fMc_StartMomentumY, "mc_startmomentumy/F");
    fMCParticlesTree->Branch("mc_startmomentumz", &fMc_StartMomentumZ, "mc_startmomentumz/F");
    fMCParticlesTree->Branch("mc_length", &fMc_Length, "mc_length/F");
    fMCParticlesTree->Branch("mc_length_tpc", &fMc_LengthTPC, "mc_length_tpc/F");
  }

  //// Tree for beam flashes
  fOpBeamFlashesTree = tfs->make<TTree>("OpBeamFlashes", "OpBeamFlashes Tree");
  fOpBeamFlashesTree->Branch("event", &fEvent, "event/i");
  fOpBeamFlashesTree->Branch("run", &fRun, "run/i");
  fOpBeamFlashesTree->Branch("subrun", &fSubrun, "subrun/i");
  fOpBeamFlashesTree->Branch("dataset_prescale_factor", &fDatasetPrescaleFactor, "dataset_prescale_factor/F");
  fOpBeamFlashesTree->Branch("num_flashes", &fNumOpBeamFlashes, "num_flashes/i");
  fOpBeamFlashesTree->Branch("flash_time", &fFlash_Time, "flash_time/F");
  fOpBeamFlashesTree->Branch("flash_totalPE", &fFlash_TotalPE, "flash_total_PE/i");
  fOpBeamFlashesTree->Branch("flash_z", &fFlash_Z, "flash_z/F");
  fOpBeamFlashesTree->Branch("flash_sz", &fFlash_sigmaZ, "flash_sz/F");
  fOpBeamFlashesTree->Branch("flash_y", &fFlash_Y, "flash_y/F");
  fOpBeamFlashesTree->Branch("flash_sy", &fFlash_sigmaY, "flash_sy/F");
  fOpBeamFlashesTree->Branch("flash_width", &fFlash_Width, "flash_width/F");
  fOpBeamFlashesTree->Branch("flash_abstime", &fFlash_AbsTime, "flash_abstime/F");
  fOpBeamFlashesTree->Branch("flash_num_PMT10percent", &fFlash_num10percentPMT, "flash_num_PMT10percent/i");

  fSimpleBeamFlashesTree = tfs->make<TTree>("SimpleBeamFlashes", "SimpleBeamFlashes Tree");
  fSimpleBeamFlashesTree->Branch("event", &fEvent, "event/i");
  fSimpleBeamFlashesTree->Branch("run", &fRun, "run/i");
  fSimpleBeamFlashesTree->Branch("subrun", &fSubrun, "subrun/i");
  fSimpleBeamFlashesTree->Branch("dataset_prescale_factor", &fDatasetPrescaleFactor, "dataset_prescale_factor/F");
  fSimpleBeamFlashesTree->Branch("num_flashes", &fNumSimpleBeamFlashes, "num_flashes/i");
  fSimpleBeamFlashesTree->Branch("flash_time", &fFlash_Time, "flash_time/F");
  fSimpleBeamFlashesTree->Branch("flash_totalPE", &fFlash_TotalPE, "flash_total_PE/i");
  fSimpleBeamFlashesTree->Branch("flash_z", &fFlash_Z, "flash_z/F");
  fSimpleBeamFlashesTree->Branch("flash_sz", &fFlash_sigmaZ, "flash_sz/F");
  fSimpleBeamFlashesTree->Branch("flash_y", &fFlash_Y, "flash_y/F");
  fSimpleBeamFlashesTree->Branch("flash_sy", &fFlash_sigmaY, "flash_sy/F");
  fSimpleBeamFlashesTree->Branch("flash_width", &fFlash_Width, "flash_width/F");
  fSimpleBeamFlashesTree->Branch("flash_abstime", &fFlash_AbsTime, "flash_abstime/F");
  fSimpleBeamFlashesTree->Branch("flash_num_PMT10percent", &fFlash_num10percentPMT, "flash_num_PMT10percent/i");

  //// Tree for cosmic flashes
  fOpCosmicFlashesTree = tfs->make<TTree>("OpCosmicFlashes", "OpCosmicFlashes Tree");
  fOpCosmicFlashesTree->Branch("event", &fEvent, "event/i");
  fOpCosmicFlashesTree->Branch("run", &fRun, "run/i");
  fOpCosmicFlashesTree->Branch("subrun", &fSubrun, "subrun/i");
  fOpCosmicFlashesTree->Branch("dataset_prescale_factor", &fDatasetPrescaleFactor, "dataset_prescale_factor/F");
  fOpCosmicFlashesTree->Branch("num_flashes", &fNumOpCosmicFlashes, "num_flashes/i");
  fOpCosmicFlashesTree->Branch("flash_time", &fFlash_Time, "flash_time/F");
  fOpCosmicFlashesTree->Branch("flash_totalPE", &fFlash_TotalPE, "flash_total_PE/i");
  fOpCosmicFlashesTree->Branch("flash_z", &fFlash_Z, "flash_z/F");
  fOpCosmicFlashesTree->Branch("flash_sz", &fFlash_sigmaZ, "flash_sz/F");
  fOpCosmicFlashesTree->Branch("flash_y", &fFlash_Y, "flash_y/F");
  fOpCosmicFlashesTree->Branch("flash_sy", &fFlash_sigmaY, "flash_sy/F");
  fOpCosmicFlashesTree->Branch("flash_width", &fFlash_Width, "flash_width/F");
  fOpCosmicFlashesTree->Branch("flash_abstime", &fFlash_AbsTime, "flash_abstime/F");
  fOpCosmicFlashesTree->Branch("flash_num_PMT10percent", &fFlash_num10percentPMT, "flash_num_PMT10percent/i");

  fSimpleCosmicFlashesTree = tfs->make<TTree>("SimpleCosmicFlashes", "SimpleCosmicFlashes Tree");
  fSimpleCosmicFlashesTree->Branch("event", &fEvent, "event/i");
  fSimpleCosmicFlashesTree->Branch("run", &fRun, "run/i");
  fSimpleCosmicFlashesTree->Branch("subrun", &fSubrun, "subrun/i");
  fSimpleCosmicFlashesTree->Branch("dataset_prescale_factor", &fDatasetPrescaleFactor, "dataset_prescale_factor/F");
  fSimpleCosmicFlashesTree->Branch("num_flashes", &fNumSimpleCosmicFlashes, "num_flashes/i");
  fSimpleCosmicFlashesTree->Branch("flash_time", &fFlash_Time, "flash_time/F");
  fSimpleCosmicFlashesTree->Branch("flash_totalPE", &fFlash_TotalPE, "flash_total_PE/i");
  fSimpleCosmicFlashesTree->Branch("flash_z", &fFlash_Z, "flash_z/F");
  fSimpleCosmicFlashesTree->Branch("flash_sz", &fFlash_sigmaZ, "flash_sz/F");
  fSimpleCosmicFlashesTree->Branch("flash_y", &fFlash_Y, "flash_y/F");
  fSimpleCosmicFlashesTree->Branch("flash_sy", &fFlash_sigmaY, "flash_sy/F");
  fSimpleCosmicFlashesTree->Branch("flash_width", &fFlash_Width, "flash_width/F");
  fSimpleCosmicFlashesTree->Branch("flash_abstime", &fFlash_AbsTime, "flash_abstime/F");
  fSimpleCosmicFlashesTree->Branch("flash_num_PMT10percent", &fFlash_num10percentPMT, "flash_num_PMT10percent/i");

  //// Tree for the PF particles
  fPFParticlesTree = tfs->make<TTree>("PFParticles", "PFParticles Tree");
  fPFParticlesTree->Branch("event", &fEvent, "event/i");
  fPFParticlesTree->Branch("run", &fRun, "run/i");
  fPFParticlesTree->Branch("subrun", &fSubrun, "subrun/i");
  fPFParticlesTree->Branch("num_pfp", &fNumPfp, "num_pfp/i");
  fPFParticlesTree->Branch("num_mcp", &fNumMcp, "num_mcp/i");
  fPFParticlesTree->Branch("num_mcp_saved", &fNumMcp_saved, "num_mcp_saved/i");
  fPFParticlesTree->Branch("num_flashes", &fNumSimpleBeamFlashes, "num_flashes/i");
  fPFParticlesTree->Branch("pdg_code", &fPdgCode, "pdg_code/I");
  fPFParticlesTree->Branch("num_daughters", &fNumDaughters, "num_daughters/i");
  fPFParticlesTree->Branch("is_primary", &fIsPrimary, "is_primary/O");
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
}

void Prototyping::reconfigure(fhicl::ParameterSet const &p)
{
  m_pfp_producer = p.get<std::string>("pfp_producer", "pandoraNu");
  m_cosmic_simpleflash_producer = p.get<std::string>("cosmic_simpleflash_producer", "simpleFlashCosmic");
  m_cosmic_opflash_producer = p.get<std::string>("cosmic_opflash_producer", "opflashCosmic");
  m_cosmic_ophit_producer = p.get<std::string>("cosmic_ophit_producer", "ophitCosmic");
  m_beam_simpleflash_producer = p.get<std::string>("beam_simpleflash_producer", "simpleFlashBeam");
  m_beam_opflash_producer = p.get<std::string>("beam_opflash_producer", "opflashBeam");
  m_beam_ophit_producer = p.get<std::string>("beam_ophit_producer", "ophitBeam");

  m_verb = p.get<bool>("verbose_output", true);
  m_is_lite = p.get<bool>("is_lite", true);
  m_is_data = p.get<bool>("is_data", false);
}

// Clear once per event
void Prototyping::clear()
{
  //fRun_sr = 0;
  //fSubrun_sr = 0;
  //fPot = 0;
  fDatasetPrescaleFactor = 1;
  //fNevents = 0; we do not want to clear this, just count all events.

  fRun = 0;
  fSubrun = 0;
  fEvent = 0;

  fNumPfp = 0;
  fNumMcp = 0;
  fNumMcp_saved = 0;
  fNumSimpleBeamFlashes = 0;
  fNumOpBeamFlashes = 0;
  fNumSimpleCosmicFlashes = 0;
  fNumOpCosmicFlashes = 0;
}

// Clear once per cluster
void Prototyping::clear_Cluster()
{
  // fClustersTree;
  fClusterCharge = -9999;
  fClusterWidth = -9999;
  fClusterPosition = -9999;
  fClusterNhits = 0;
  fClusterPlane = 0;
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

// Clear once per MC particle
void Prototyping::clear_MCParticle()
{
  fMc_EndX = -9999;
  fMc_EndY = -9999;
  fMc_EndZ = -9999;
  fMc_StartX = -9999;
  fMc_StartY = -9999;
  fMc_StartZ = -9999;
  fMc_Length = -9999;
  fMc_StartMomentumX = 0;
  fMc_StartMomentumY = 0;
  fMc_StartMomentumZ = 0;
  fMc_Process = 0;
  fMc_E = 0;
  fMc_PdgCode = 0;
  fMc_Time = 0;
  fMc_StatusCode = -1; // This is always 1 normally, so put it -1 if something is wrong
  fMc_EndInside = false;
  fMc_StartInside = false;
  fMc_PartInside = true; // default we say there is a part inside, it gets set to false if no intersections are found.
  fMc_kBeamNeutrino = false;
  fMc_StartX_tpc = -9999;
  fMc_StartY_tpc = -9999;
  fMc_StartZ_tpc = -9999;
  fMc_EndX_tpc = -9999;
  fMc_EndY_tpc = -9999;
  fMc_EndZ_tpc = -9999;
  fMc_LengthTPC = -9999;
}

// Clear once per beam flash
void Prototyping::clear_Flashes()
{
  fFlash_Time = -9999;
  fFlash_TotalPE = 0;
  fFlash_Z = -9999;
  fFlash_sigmaZ = -9999;
  fFlash_Y = -9999;
  fFlash_sigmaY = -9999;
  fFlash_num10percentPMT = 0;
  fFlash_Width = -9999;
  fFlash_AbsTime = -9999;
}

art::Ptr<simb::MCTruth> Prototyping::TrackIDToMCTruth(art::Event const & e, std::string _geant_producer, int geant_track_id) {

    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;

    lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geant_producer, truthToParticles, particlesToTruth);

    for (auto iter : particlesToTruth) {
      if (iter.first->TrackId() == geant_track_id) {
        return iter.second;
      }
    }

    art::Ptr<simb::MCTruth> null_ptr;
    return null_ptr;
}
