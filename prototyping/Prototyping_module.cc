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
}

void Prototyping::analyze(art::Event const &evt)
{
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  fNevents++;

  // Initialise the handles and associations
  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(m_pfp_producer);
  auto const &spacepoint_handle =
      evt.getValidHandle<std::vector<recob::SpacePoint>>(m_pfp_producer);

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
      if (mcparticle.E() > 0.01)
      {
        fNumMcp10MeV++; //Counter
        clear_MCParticle();
        fMc_Process = map_process[mcparticle.Process()];
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

        std::vector<float> start = { fMc_StartX, fMc_StartY, fMc_StartZ };
        std::vector<float> end = { fMc_EndX, fMc_EndY, fMc_EndZ };        
        fMc_StartInside = geoHelper.isActive(start);
        fMc_EndInside = geoHelper.isActive(end);
        fMc_Length = geoHelper.distance(start,end); // This is the total length, not the length in the detector!
        string_process.insert(mcparticle.Process());

        fMCParticlesTree->Fill();
      }
    }

    std::cout << "string_process has mamebers: " << string_process.size() << std::endl;
    for (auto elem : string_process)
    {
      std::cout << elem << ", ";
    }
  }

  fNumPfp = pfparticle_handle->size();
  for (size_t i_pfp = 0; i_pfp < fNumPfp; i_pfp++)
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
        fPlane = hit->WireID().Plane;
        if (fPlane == 0)
        {
          fNhitsU += 1;
        }
        else if (fPlane == 1)
        {
          fNhitsV += 1;
        }
        else if (fPlane == 2)
        {
          fNhitsY += 1;
        }

        if (m_is_lite == false)
        {
          fWire = hit->WireID().Wire;
          fCharge = hit->Integral();
          fHitsTree->Fill();
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

    if (m_is_lite == false)
    {
      // Space points
      std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(i_pfp);
      for (art::Ptr<recob::SpacePoint> &sps : spcpnts)
      {
        auto xyz = sps->XYZ();
        fx = xyz[0];
        fy = xyz[1];
        fz = xyz[2];

        std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(sps.key());
        fChargeU = 0;
        fChargeV = 0;
        fChargeY = 0;
        for (art::Ptr<recob::Hit> &hit : hits)
        {
          double hit_plane = hit->WireID().Plane;
          double hit_integral = hit->Integral();
          if (hit_plane == 0)
            fChargeU += hit_integral;
          else if (hit_plane == 1)
            fChargeV += hit_integral;
          else if (hit_plane == 2)
            fChargeY += hit_integral;
          else
            std::cout << "hit plane != 0, 1, 2, but " << hit_plane << std::endl;
        }
        fSpacePointsTree->Fill();
      }
    }
  }
}

DEFINE_ART_MODULE(Prototyping)
