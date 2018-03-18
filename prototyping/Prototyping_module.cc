#include "Prototyping.h"

void ProtoTypeAnalyzer::analyze(art::Event const &evt)
{
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();

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


  for (size_t i_pfp = 0; i_pfp < pfparticle_handle->size(); i_pfp++)
  {
    recob::PFParticle const &pfparticle = pfparticle_handle->at(i_pfp);
    fPdgCode = pfparticle.PdgCode();

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

    try
    {
      art::Ptr<recob::Vertex> vertex_obj = vertex_per_pfpart.at(i_pfp);
      double neutrino_vertex[3];
      vertex_obj->XYZ(neutrino_vertex);
      fvx = neutrino_vertex[0];
      fvy = neutrino_vertex[1];
      fvz = neutrino_vertex[2];
    }
    catch (...)
    {
      fvx = 1000000.;
      fvy = 1000000.;
      fvz = 1000000.;
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

DEFINE_ART_MODULE(ProtoTypeAnalyzer)
