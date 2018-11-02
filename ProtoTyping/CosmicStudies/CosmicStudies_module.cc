#include "CosmicStudies.h"

void CosmicStudies::analyze(art::Event const &evt)
{
  clear();
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  fNevents++;

  //// Filling the MC and TPCreco information -----------------------------------------------------
  if (!m_is_data)
  {
    if (!m_simmed)
    {
      pandoraHelper.Configure(evt, m_pfp_producer, m_pfp_producer, "gaushit", "largeant");
    }
    else
    {
      pandoraHelper.Configure(evt, m_pfp_producer, m_pfp_producer, "gaushit", "largeant", "gaushit");
    }
    fill_MC(evt);
  }
  fill_TPCreco(evt);
  ////---------------------------------------------------------------------------------

  //// Filling the flashes ------------------------------------------------------------
  art::ValidHandle<std::vector<recob::OpFlash>> const &simple_beam_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_beam_simpleflash_producer);
  fNumSimpleBeamFlashes = simple_beam_handle->size();
  std::cout << "[ProtoTyping] fNumSimpleBeamFlashes: " << fNumSimpleBeamFlashes << std::endl;
  fill_flash(simple_beam_handle, fNumSimpleBeamFlashes, fSimpleBeamFlashesTree);

  art::ValidHandle<std::vector<recob::OpFlash>> const &op_beam_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_beam_opflash_producer);
  fNumOpBeamFlashes = op_beam_handle->size();
  std::cout << "[ProtoTyping] fNumOpBeamFlashes: " << fNumOpBeamFlashes << std::endl;
  fill_flash(op_beam_handle, fNumOpBeamFlashes, fOpBeamFlashesTree);

  art::ValidHandle<std::vector<recob::OpFlash>> const &simple_cosmic_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_cosmic_simpleflash_producer);
  fNumSimpleCosmicFlashes = simple_cosmic_handle->size();
  std::cout << "[ProtoTyping] fNumSimpleCosmicFlashes: " << fNumSimpleCosmicFlashes << std::endl;
  fill_flash(simple_cosmic_handle, fNumSimpleCosmicFlashes, fSimpleCosmicFlashesTree);

  art::ValidHandle<std::vector<recob::OpFlash>> const &op_cosmic_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_cosmic_opflash_producer);
  fNumOpCosmicFlashes = op_cosmic_handle->size();
  std::cout << "[ProtoTyping] fNumOpCosmicFlashes: " << fNumOpCosmicFlashes << std::endl;
  fill_flash(op_cosmic_handle, fNumOpCosmicFlashes, fOpCosmicFlashesTree);
  ////---------------------------------------------------------------------------------

  fEventTree->Fill();
}

void CosmicStudies::endSubRun(const art::SubRun &sr)
{
  fRun_sr = sr.run();
  fSubrun_sr = sr.subRun();

  art::Handle<sumdata::POTSummary> potListHandle;
  if (!m_is_data)
  {
    if (sr.getByLabel("generator", potListHandle))
    {
      fPot = potListHandle->totpot;
      std::cout << "[CosmicStudies::endSubRun] POT for SubRun: " << fPot << std::endl;
    }
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

  // Reset events counter at the end of the subrun
  fNevents = 0;
}

void CosmicStudies::fill_flash(art::ValidHandle<std::vector<recob::OpFlash>> const &flash_handle, uint number, TTree *tree)
{
  float prevTime = 0; // For the first flash, the flashtime will be identical to the flashtimediff
  for (uint ifl = 0; ifl < number; ++ifl)
  {
    clear_Flashes();

    recob::OpFlash const &flash = flash_handle->at(ifl);
    fFlash_TotalPE = flash.TotalPE();
    fFlash_Time = flash.Time();
    fFlash_DiffTime = fFlash_Time - prevTime;
    prevTime = fFlash_Time;
    fFlash_Y = flash.YCenter();
    fFlash_Z = flash.ZCenter();
    fFlash_sigmaY = flash.YWidth();
    fFlash_sigmaZ = flash.ZWidth();
    fFlash_AbsTime = flash.AbsTime();
    fFlash_Width = flash.TimeWidth();

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

void CosmicStudies::fill_MC(art::Event const &evt)
{
  if (m_is_true_nu)
  {
    auto const &generator_handle = evt.getValidHandle<std::vector<simb::MCTruth>>("generator");
    auto const &generator(*generator_handle);
    fNum_nu = generator.size();

    if (m_verb)
    {
      std::cout << "[CosmicStudies] True neutrinos found: " << fNum_nu << std::endl;
    }

    for (auto &gen : generator)
    {

      if (gen.Origin() == simb::kBeamNeutrino)
      {
        fNu_pdg_code.push_back(gen.GetNeutrino().Nu().PdgCode());
        fNu_time.push_back(gen.GetNeutrino().Nu().T());
        fNu_E.push_back(gen.GetNeutrino().Nu().E());
        fNu_ccnc.push_back(gen.GetNeutrino().CCNC());
        fNu_vtx_x.push_back(gen.GetNeutrino().Nu().Vx());
        fNu_vtx_y.push_back(gen.GetNeutrino().Nu().Vy());
        fNu_vtx_z.push_back(gen.GetNeutrino().Nu().Vz());
      }
    }
  }

  auto const &mcparticles_handle = evt.getValidHandle<std::vector<simb::MCParticle>>("largeant");
  geo::TPCGeo const &thisTPC = geo_service->TPC();
  const geo::BoxBoundedGeo &theTpcGeo(thisTPC.ActiveBoundingBox());
  //std::cout << "This is a geometry test: " << theTpcGeo.MaxX() << ", " << theTpcGeo.MaxY() << ", " << theTpcGeo.MaxZ() << std::endl;
  //[-1.55254.8, 117.47], 1036.9 deviates slightly from what we expect, but ok!

  fNumMcp = mcparticles_handle->size();

  for (size_t i_mcp = 0; i_mcp < fNumMcp; i_mcp++)
  {
    simb::MCParticle const &mcparticle = mcparticles_handle->at(i_mcp);
    string_process.insert(mcparticle.Process());
    // Important, only save MC particles with energy over 100MeV, THIS WILL SAVE ALL MUONS.
    uint pdg = abs(mcparticle.PdgCode());
    bool pdg_ok = pdg == 11 or pdg == 13 or pdg == 211 or pdg == 111 or pdg == 22 or pdg == 2112 or pdg == 2212;
    if (mcparticle.E() > 0.1 && pdg_ok)
    {
      clear_MCParticle();
      if (map_process.find(mcparticle.Process()) != map_process.end())
      {
        fMc_Process = map_process[mcparticle.Process()];
      }
      else
      {
        std::cout << "[CosmicStudies::analyze] New MC interaction process found!" << std::endl;
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

      TVector3 mc_start = mcparticle.Position().Vect();

      if (m_is_true_nu)
      {
        // Is this MC particle neutrino?
        const art::Ptr<simb::MCTruth> mctruth = pandoraHelper.TrackIDToMCTruth(evt, "largeant", mcparticle.TrackId());
        if (mctruth->Origin() == simb::kBeamNeutrino)
        {
          fMc_kBeamNeutrino = true;
          if (m_verb)
          {
            std::cout << "MCParticle coming from a neutrino found! " << fMc_PdgCode << std::endl;
          }
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
        std::vector<TVector3> intersections = theTpcGeo.GetIntersections(startvec, startdir);
        uint num_intersections = intersections.size();

        //Particles completely passes the TPC without entering
        if (num_intersections == 0)
        {
          fMc_PartInside = false;
          fMc_LengthTPC = 0;
        }

        // Particle started inside the TPC and is exiting
        else if (num_intersections == 1)
        {
          fMc_StartX_tpc = fMc_StartX;
          fMc_StartY_tpc = fMc_StartY;
          fMc_StartZ_tpc = fMc_StartZ;
          fMc_EndX_tpc = intersections[0].X();
          fMc_EndY_tpc = intersections[0].Y();
          fMc_EndZ_tpc = intersections[0].Z();
        }

        // Particle crosses TPC or particle stops in TPC or particle stops before TPC
        else if (num_intersections == 2)
        {
          double len_start_cross0 = (mc_start - intersections[0]).Mag();
          double len_start_cross1 = (mc_start - intersections[1]).Mag();

          // Particle is crossing the TPC
          if (std::max(len_start_cross0, len_start_cross1) < fMc_Length)
          {
            fMc_StartX_tpc = intersections[0].X();
            fMc_StartY_tpc = intersections[0].Y();
            fMc_StartZ_tpc = intersections[0].Z();
            fMc_EndX_tpc = intersections[1].X();
            fMc_EndY_tpc = intersections[1].Y();
            fMc_EndZ_tpc = intersections[1].Z();
          }
          // Particle stops before entering
          else if (std::min(len_start_cross0, len_start_cross1) > fMc_Length)
          {
            fMc_PartInside = false;
            fMc_LengthTPC = 0;
          }
          // Particle stops inside (assume intersections are oredered)
          else
          {
            fMc_EndX_tpc = fMc_EndX;
            fMc_EndY_tpc = fMc_EndY;
            fMc_EndZ_tpc = fMc_EndZ;
            if (len_start_cross0 < fMc_Length)
            {
              fMc_StartX_tpc = intersections[0].X();
              fMc_StartY_tpc = intersections[0].Y();
              fMc_StartZ_tpc = intersections[0].Z();
            }
            else
            {
              fMc_StartX_tpc = intersections[1].X();
              fMc_StartY_tpc = intersections[1].Y();
              fMc_StartZ_tpc = intersections[1].Z();
            }
          }
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
        auto const &SCE(*lar::providerFrom<spacecharge::SpaceChargeService>());
        auto sce_start = SCE.GetPosOffsets(geo::Point_t(fMc_StartX_tpc, fMc_StartY_tpc, fMc_StartZ_tpc));
        auto sce_end = SCE.GetPosOffsets(geo::Point_t(fMc_EndX_tpc, fMc_EndY_tpc, fMc_EndZ_tpc));

        fMc_StartX_sce = fMc_StartX_tpc - sce_start.X() + 0.7;
        fMc_StartY_sce = fMc_StartY_tpc + sce_start.Y();
        fMc_StartZ_sce = fMc_StartZ_tpc + sce_start.Z();
        fMc_EndX_sce = fMc_EndX_tpc - sce_end.X() + 0.7;
        fMc_EndY_sce = fMc_EndY_tpc + sce_end.Y();
        fMc_EndZ_sce = fMc_EndZ_tpc + sce_end.Z();
      }
      fMCParticlesTree->Fill();
      fNumMcp_saved++; //Counter
    }
  }
}

void CosmicStudies::fill_TPCreco(art::Event const &evt)
{
  auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);
  std::vector<art::Ptr<recob::PFParticle>> pfp_v;
  art::fill_ptr_vector(pfp_v, pfparticle_handle);

  auto const &cluster_handle = evt.getValidHandle<std::vector<recob::Cluster>>(m_pfp_producer);
  auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(m_pfp_producer);

  art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_cluster(cluster_handle, evt, m_pfp_producer);
  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, m_pfp_producer);
  art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, m_pfp_producer);

  if (!m_is_data)
  {
    lar_pandora::PFParticlesToMCParticles matchedParticles;
    pandoraHelper.GetRecoToTrueMatches(matchedParticles);
    if (m_verb)
    {
      std::cout << "[CosmicStudies::fill_TPCreco] ";
      std::cout << "PFParticlesToMCParticles constructed: Number of PFPparticles matched: " << matchedParticles.size() << std::endl;
      std::cout << "[CosmicStudies::fill_TPCreco] ";
      std::cout << "Number of PFPparticles in event: " << pfp_v.size() << std::endl;
    }
  }
  /*
  if (m_verb)
  {
    std::cout << "[CosmicStudies::fill_TPCreco] ";
    std::cout << "PFParticlesToMCParticles constructed: Number of PFPparticles matched: " << matchedParticles.size() << std::endl;
    std::cout << "[CosmicStudies::fill_TPCreco] ";
    std::cout << "Number of PFPparticles in event: " << pfp_v.size() << std::endl;
    if (matchedParticles.size() > 0 && pfp_v.size() > 0)
    {
      art::Ptr<recob::PFParticle> firstPFP = pfp_v.at(0);
      if (matchedParticles.find(firstPFP) == matchedParticles.end())
      {
        std::cout << "PFParticle is not matched" << std::endl;
      }
      else
      {
        art::Ptr<simb::MCParticle> firstMC = matchedParticles[firstPFP];
        std::cout << "PFParticle is matched to mcparticle with pdgcode:" << firstMC->PdgCode() << std::endl;
      }
    }
  }
*/
  /// PFParticle Tree
  fNumPfp = pfp_v.size();
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
}
