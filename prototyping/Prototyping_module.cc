#include "Prototyping.h"

void Prototyping::endSubRun(const art::SubRun &sr)
{

  fRun_sr = sr.run();
  fSubrun_sr = sr.subRun();

  art::Handle<sumdata::POTSummary> potListHandle;
  if (!m_is_data)
  {
    if (sr.getByLabel("generator", potListHandle))
      fPot = potListHandle->totpot;
    else
      fPot = 0.;
  }
  else
  {
    if (sr.getByLabel("beamdata", "bnbETOR860", potListHandle))
      fPot = potListHandle->totpot;
    else
      fPot = 0.;
  }

  fPOTTree->Fill();

  if (m_verb)
  {
    std::cout << "string_process has mamebers: " << string_process.size() << std::endl;
    for (auto elem : string_process)
    {
      std::cout << elem << ", ";
    }
  }
}

void Prototyping::analyze(art::Event const &evt)
{
  clear();
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  fNevents++;

  // Check if there is a prescale factor for a certain data trigger and stor it
  if (m_is_data)
  {
    art::Handle<raw::ubdaqSoftwareTriggerData> SWTriggerHandle;
    evt.getByLabel("daq", SWTriggerHandle);
    fDatasetPrescaleFactor = SWTriggerHandle->getPrescale("EXT_unbiased_PrescaleAlgo");
    std::cout << "[Prototyping constructor] Prescale factor for  EXT_unbiased_PrescaleAlgo is " << fDatasetPrescaleFactor << std::endl;
  }

  // Initialise the handles and associations
  auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);
  auto const &cluster_handle = evt.getValidHandle<std::vector<recob::Cluster>>(m_pfp_producer);
  auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(m_pfp_producer);

  art::ValidHandle<std::vector<recob::OpFlash>> const &simple_cosmic_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_cosmic_simpleflash_producer);
  art::ValidHandle<std::vector<recob::OpFlash>> const &op_cosmic_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_cosmic_opflash_producer);
  art::ValidHandle<std::vector<recob::OpFlash>> const &simple_beam_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_beam_simpleflash_producer);
  art::ValidHandle<std::vector<recob::OpFlash>> const &op_beam_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_beam_opflash_producer);

  art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_cluster(cluster_handle, evt, m_pfp_producer);
  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, m_pfp_producer);
  art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, m_pfp_producer);

  if (!m_is_data)
  {
    auto const &mcparticles_handle =
        evt.getValidHandle<std::vector<simb::MCParticle>>("largeant");

    fNumMcp = mcparticles_handle->size();

    for (size_t i_mcp = 0; i_mcp < fNumMcp; i_mcp++)
    {
      simb::MCParticle const &mcparticle = mcparticles_handle->at(i_mcp);
      string_process.insert(mcparticle.Process());
      if (mcparticle.E() > 0.02)
      {
        clear_MCParticle();
        if (map_process.find(mcparticle.Process()) != map_process.end())
        {
          fMc_Process = map_process[mcparticle.Process()];
        }
        else
        {
          std::cout << "[Prototyping::analyze] New MC interaction process found!" << std::endl;
        }
        fMc_Time = mcparticle.T();
        fMc_StatusCode = mcparticle.StatusCode();
        fMc_E = mcparticle.E();
        fMc_PdgCode = mcparticle.PdgCode();
        fMc_StartMomentumX = mcparticle.Px();
        fMc_StartMomentumY = mcparticle.Py();
        fMc_StartMomentumZ = mcparticle.Pz();
        fMc_StartX = mcparticle.Vx();
        fMc_StartY = mcparticle.Vy();
        fMc_StartZ = mcparticle.Vz();
        fMc_EndX = mcparticle.EndX();
        fMc_EndY = mcparticle.EndY();
        fMc_EndZ = mcparticle.EndZ();

        // Is this MC particle neutrino? This does seemingly not work!
        const art::Ptr<simb::MCTruth> mctruth = TrackIDToMCTruth(evt, "largeant", mcparticle.TrackId());
        //const art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth(mcparticle.TrackId());
        if (mctruth->Origin() == simb::kBeamNeutrino)
        {
          fMc_kBeamNeutrino = true;
          if(m_verb){
            std::cout<<"MCParticle coming from a neutrino found! " << fMc_PdgCode << std::endl; 
          }
        }

        std::vector<float> start = {fMc_StartX, fMc_StartY, fMc_StartZ};
        std::vector<float> end = {fMc_EndX, fMc_EndY, fMc_EndZ};
        fMc_StartInside = geoHelper.isActive(start);
        fMc_EndInside = geoHelper.isActive(end);
        fMc_Length = geoHelper.distance(start, end); // This is the total length, not the length in the detector!

        //Find the section that is inside the tpc
        if (!fMc_StartInside || !fMc_EndInside)
        {
          TVector3 startvec(fMc_StartX, fMc_StartY, fMc_StartZ);
          TVector3 startdir(fMc_StartMomentumX, fMc_StartMomentumY, fMc_StartMomentumZ);
          geo::TPCGeo const &thisTPC = geo->TPC();
          //std::cout << "thisTPC.ActiveHalfHeight(): " << thisTPC.ActiveHalfHeight() << std::endl;
          //std::cout << "thisTPC.HalfHeight(): " << thisTPC.HalfHeight() << std::endl;
          const geo::BoxBoundedGeo &theTpcGeo(thisTPC);
          //std::cout << "This is a geometry test: " << theTpcGeo.MaxX() << ", " << theTpcGeo.MaxY() << ", " << theTpcGeo.MaxZ() << std::endl;
          std::vector<TVector3> intersections = theTpcGeo.GetIntersections(startvec, startdir);
          uint num_intersections = intersections.size();
          if (num_intersections == 0)
          {
            fMc_PartInside = false;
            fMc_LengthTPC = 0;
          }
          else if (num_intersections == 1)
          {
            if (fMc_StartInside)
            {
              fMc_StartX_tpc = fMc_StartX;
              fMc_StartY_tpc = fMc_StartY;
              fMc_StartZ_tpc = fMc_StartZ;
              fMc_EndX_tpc = intersections[0].X();
              fMc_EndY_tpc = intersections[0].Y();
              fMc_EndZ_tpc = intersections[0].Z();
            }
            else
            {
              fMc_EndX_tpc = fMc_EndX;
              fMc_EndY_tpc = fMc_EndY;
              fMc_EndZ_tpc = fMc_EndZ;
              fMc_StartX_tpc = intersections[0].X();
              fMc_StartY_tpc = intersections[0].Y();
              fMc_StartZ_tpc = intersections[0].Z();
            }
          }
          else if (num_intersections == 2)
          {
            fMc_StartX_tpc = intersections[0].X();
            fMc_StartY_tpc = intersections[0].Y();
            fMc_StartZ_tpc = intersections[0].Z();
            fMc_EndX_tpc = intersections[1].X();
            fMc_EndY_tpc = intersections[1].Y();
            fMc_EndZ_tpc = intersections[1].Z();
          }
        }
        else //Start and end are inside
        {
          fMc_StartX_tpc = fMc_StartX;
          fMc_StartY_tpc = fMc_StartY;
          fMc_StartZ_tpc = fMc_StartZ;
          fMc_EndX_tpc = fMc_EndX;
          fMc_EndY_tpc = fMc_EndY;
          fMc_EndZ_tpc = fMc_EndZ;
        }

        if (fMc_PartInside)
        {
          std::vector<float> start_tpc = {fMc_StartX_tpc, fMc_StartY_tpc, fMc_StartZ_tpc};
          std::vector<float> end_tpc = {fMc_EndX_tpc, fMc_EndY_tpc, fMc_EndZ_tpc};
          fMc_LengthTPC = geoHelper.distance(start_tpc, end_tpc);

          //Correct the inside tpcpoints for spacecharge
          auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
          std::vector<double> sce_start = SCE->GetPosOffsets(fMc_StartX_tpc, fMc_StartY_tpc, fMc_StartZ_tpc); // Seemingly this results in 0
          std::vector<double> sce_end = SCE->GetPosOffsets(fMc_EndX_tpc, fMc_EndY_tpc, fMc_EndZ_tpc);
          //std::cout << "start_tpc " << fMc_StartX_tpc << ", " << fMc_StartY_tpc << ", " << fMc_StartZ_tpc << std::endl;
          //std::cout << "sce_start " << SCE->GetPosOffsets(fMc_StartX_tpc, fMc_StartY_tpc, fMc_StartZ_tpc)[0] << ", " << sce_start[1] << ", " << sce_start[2] << std::endl;

          fMc_StartX_sce = fMc_StartX_tpc - sce_start[0] + 0.7;
          fMc_StartY_sce = fMc_StartY_tpc + sce_start[1];
          fMc_StartZ_sce = fMc_StartZ_tpc + sce_start[2];
          fMc_EndX_sce = fMc_EndX_tpc - sce_end[0] + 0.7;
          fMc_EndY_sce = fMc_EndY_tpc + sce_end[1];
          fMc_EndZ_sce = fMc_EndZ_tpc + sce_end[2];

          fMCParticlesTree->Fill();
          fNumMcp_saved++; //Counter
        }
      }
    }
  }

  //// Filling the flashes
  fNumSimpleBeamFlashes = simple_beam_handle->size();
  std::cout << "[ProtoTyping] fNumSimpleBeamFlashes: " << fNumSimpleBeamFlashes << std::endl;
  fill_flash(simple_beam_handle, fNumSimpleBeamFlashes, fSimpleBeamFlashesTree);

  fNumOpBeamFlashes = op_beam_handle->size();
  std::cout << "[ProtoTyping] fNumOpBeamFlashes: " << fNumOpBeamFlashes << std::endl;
  fill_flash(op_beam_handle, fNumOpBeamFlashes, fOpBeamFlashesTree);

  fNumSimpleCosmicFlashes = simple_cosmic_handle->size();
  std::cout << "[ProtoTyping] fNumSimpleCosmicFlashes: " << fNumSimpleCosmicFlashes << std::endl;
  fill_flash(simple_cosmic_handle, fNumSimpleCosmicFlashes, fSimpleCosmicFlashesTree);

  fNumOpCosmicFlashes = op_cosmic_handle->size();
  std::cout << "[ProtoTyping] fNumOpCosmicFlashes: " << fNumOpCosmicFlashes << std::endl;
  fill_flash(op_cosmic_handle, fNumOpCosmicFlashes, fOpCosmicFlashesTree);

  /// PF Particle Tree
  fNumPfp = pfparticle_handle->size();
  for (uint i_pfp = 0; i_pfp < fNumPfp; i_pfp++)
  {
    clear_PFParticle();

    recob::PFParticle const &pfparticle = pfparticle_handle->at(i_pfp);
    fPdgCode = pfparticle.PdgCode();
    fNumDaughters = pfparticle.NumDaughters();
    fIsPrimary = pfparticle.IsPrimary();

    fNhits = 0;
    fNhitsU = 0;
    fNhitsV = 0;
    fNhitsY = 0;
    fNclusters = 0;

    // Clusters and Hits
    std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(i_pfp);
    fNclusters = clusters.size();
    for (art::Ptr<recob::Cluster> &cluster : clusters)
    {
      clear_Cluster();

      fClusterCharge = cluster->Integral();
      fClusterWidth = cluster->Width();
      fClusterNhits = cluster->NHits();
      fClusterPlane = cluster->View();
      fClusterPosition = (cluster->EndWire() + cluster->StartWire()) / 2.;
      fClustersTree->Fill();

      fNhits += fClusterNhits;

      std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster.at(cluster.key());
      for (art::Ptr<recob::Hit> &hit : hits)
      {
        uint plane = hit->WireID().Plane;
        if (plane == 0)
        {
          fNhitsU += 1;
        }
        else if (plane == 1)
        {
          fNhitsV += 1;
        }
        else if (plane == 2)
        {
          fNhitsY += 1;
        }
      }
    }

    // Is the Pfparticle is a TRACK
    if (fPdgCode == 13)
    {
      art::Ptr<recob::Track> const &track_obj = track_per_pfpart.at(i_pfp);
      if (track_obj.isNull())
      {
        std::cout << "[LEE Analyzer] track is Null" << std::endl;
        fTrack_Valid = false;
      }
      else
      {
        fTrack_Valid = true;
        fTrack_HasMomentum = track_obj->HasMomentum();

        fTrack_StartX = track_obj->Start().X();
        fTrack_StartY = track_obj->Start().Y();
        fTrack_StartZ = track_obj->Start().Z();
        fTrack_EndX = track_obj->Trajectory().End().X();
        fTrack_EndY = track_obj->Trajectory().End().Y();
        fTrack_EndZ = track_obj->Trajectory().End().Z();
        fTrack_Length = track_obj->Length();
        fTrack_StartMomentumX = track_obj->StartMomentumVector().X();
        fTrack_StartMomentumY = track_obj->StartMomentumVector().Y();
        fTrack_StartMomentumZ = track_obj->StartMomentumVector().Z();
        fTrack_EndMomentumX = track_obj->EndMomentumVector().X();
        fTrack_EndMomentumY = track_obj->EndMomentumVector().Y();
        fTrack_EndMomentumZ = track_obj->EndMomentumVector().Z();
        fTrack_Theta = track_obj->Theta();
        fTrack_Phi = track_obj->Phi();
        fTrack_ZenithAngle = track_obj->ZenithAngle();
        fTrack_AzimuthAngle = track_obj->AzimuthAngle();
      }
    }

    try
    {
      art::Ptr<recob::Vertex> vertex_obj = vertex_per_pfpart.at(i_pfp);
      double vertex[3];
      vertex_obj->XYZ(vertex);
      fVx = vertex[0];
      fVy = vertex[1];
      fVz = vertex[2];
    }
    catch (...)
    {
      std::cout << "No vertex found for " << fPdgCode << " with " << fNhits << std::endl;
    }
    fPFParticlesTree->Fill();
  }
  fEventTree->Fill();
}

void Prototyping::fill_flash(art::ValidHandle<std::vector<recob::OpFlash>> const &flash_handle, uint number, TTree *tree)
{
  for (uint ifl = 0; ifl < number; ++ifl)
  {
    clear_Flashes();

    recob::OpFlash const &flash = flash_handle->at(ifl);
    fFlash_TotalPE = flash.TotalPE();
    fFlash_Time = flash.Time();
    fFlash_Y = flash.YCenter();
    fFlash_Z = flash.ZCenter();
    fFlash_sigmaY = flash.YWidth();
    fFlash_sigmaZ = flash.ZWidth();
    fFlash_AbsTime = flash.AbsTime();
    fFlash_Width = flash.TimeWidth();

    if (m_verb)
    {
      std::cout << "[Prototyping::fill_flash] flash: time " << flash.Time() << "\twidth " << flash.TimeWidth() << "\tPE " << flash.TotalPE() << std::endl;
    }

    for (uint i_pmt = 0; i_pmt < 32; i_pmt++)
    {
      //std::cout << i_pmt << "\t" << fFlash_TotalPE << "\t" << fFlash_TotalPE / 10.0 << "\t" << flash.PE(i_pmt) << "\t" << fFlash_num10percentPMT << std::endl;
      if (flash.PE(i_pmt) > (fFlash_TotalPE / 10.0))
      {
        fFlash_num10percentPMT++;
      }
    }
    tree->Fill();
  }
}

DEFINE_ART_MODULE(Prototyping)
