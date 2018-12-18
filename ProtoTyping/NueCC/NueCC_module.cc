#include "NueCC.h"

void NueCC::analyze(art::Event const &evt)
{
  clearEvent();
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  std::cout << "[NueCC::analyze]: Run " << fRun << ", Subrun " << fSubrun << ", Event " << fEvent << std::endl;

  larpandora.CollectPFParticleMetadata(evt, m_pfp_producer, pfparticles, particlesToMetadata);
  if (pfparticles.size() == 0)
    std::cout << "[NueCC::FillReconstructed] No reconstructed PFParticles in event." << std::endl;
  else
  {
    larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
    if (pfneutrinos.size() != 1)
      std::cout << "[NueCC::FillReconstructed] Number of reconstructed neutrinos in event is " << pfneutrinos.size() << std::endl;
    else // We have a reconstructed neutrino
    {
      FillReconstructed(evt);
      if (m_isData)
        FillTruth(evt);
      fEventTree->Fill();
    }
  }
}

void NueCC::FillReconstructed(art::Event const &evt)
{
  // Load associations and collections
  lar_pandora::VertexVector vertexVector_dummy;
  lar_pandora::PFParticleVector particleVector_dummy;
  larpandora.CollectVertices(evt, m_pfp_producer, vertexVector_dummy, particlesToVertices);
  larpandora.CollectPFParticles(evt, m_pfp_producer, particleVector_dummy, particlesToClusters);
  
  //lar_pandora::SpacePointVector spacePointVector_dummy;
  //lar_pandora::ClusterVector clusterVector_dummy;
  //larpandora.CollectPFParticles(evt, m_pfp_producer, particleVector_dummy, particlesToSpacePoints);
  //larpandora.CollectClusters(evt, m_pfp_producer, clusterVector_dummy, clustersToHits);
  //larpandora.CollectSpacePoints(evt, m_pfp_producer, spacePointVector_dummy, spacePointsToHits, hitsToSpacePoints);

  // Start filling information
  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
  fNu_PDG = pfnu->PdgCode();
  fNumPrimaryDaughters = pfnu->NumDaughters();
  lar_pandora::MetadataVector neutrino_metadata_vec = particlesToMetadata.at(pfnu);
  lar_pandora::VertexVector neutrino_vertex_vec = particlesToVertices.at(pfnu);
  if (neutrino_metadata_vec.size() != 1 || neutrino_vertex_vec.size() != 1)
  {
    std::cout << "[NueCC::FillReconstructed] Neutrino association problem." << std::endl;
  }
  else
  {
    const larpandoraobj::PFParticleMetadata::PropertiesMap &neutrino_properties = neutrino_metadata_vec.front()->GetPropertiesMap();
    fNu_Score = neutrino_properties.at("NuScore");
    fNu_SliceIndex = neutrino_properties.at("SliceIndex");
    const recob::Vertex::Point_t &neutrino_vtx = neutrino_vertex_vec.front()->position();
    fNu_Vx = neutrino_vtx.X();
    fNu_Vy = neutrino_vtx.Y();
    fNu_Vz = neutrino_vtx.Z();
  }
  larpandora.BuildPFParticleMap(pfparticles, particleMap);
  pandoraInterfaceHelper.CollectDownstreamPFParticles(particleMap, pfnu, pfdaughters);
  fNumDaughters = pfdaughters.size() - 1; // The neutrino itself is included here.
  std::cout << "[NueCC::FillReconstructed] neutrino found! PDG: " << fNu_PDG << ", Primary Daughters: " << fNumPrimaryDaughters;
  std::cout << ", Daughters: " << fNumDaughters << ", TopoScore: " << fNu_Score << std::endl;
  for (auto const pfp : pfdaughters)
  {
    if (!pfp->IsPrimary())
      if (!FillDaughters(pfp))
        fDaughtersStored = false;
  }

  // Get the map PFP->MCP and the set of MCPs
  if (!m_isData)
  {
    pandoraInterfaceHelper.Configure(evt, m_pfp_producer, m_pfp_producer, m_hitfinder_producer, m_geant_producer, m_hit_mcp_producer);
    pandoraInterfaceHelper.GetRecoToTrueMatches(matchedParticles);
    std::cout << "[NueCC::FillReconstructed] Reco-Truth matches constructed." << std::endl;
    for (auto it = matchedParticles.begin(); it != matchedParticles.end(); ++it)
    {
      matchedMCParticles.insert(it->second);
    }
    std::cout << "[NueCC::FillReconstructed] ";
    std::cout << "PFParticlesToMCParticles constructed: Number of PFPparticles matched: " << matchedParticles.size() << std::endl;
    std::cout << "[NueCC::FillReconstructed] ";
    std::cout << "Number of PFPparticles in event: " << pfparticles.size() << std::endl;
  }
}

bool NueCC::FillDaughters(const art::Ptr<recob::PFParticle> &pfp)
{
  clearDaughter();
  fGeneration = larpandora.GetGeneration(particleMap, pfp);
  if (pfp->PdgCode() == 13)
  {
    fIsTrack = true;
    fNumTracks++;
  }
  else //(pfp->PdgCode()==11)
  {
    fIsShower = true;
    fNumShowers++;
  }
  if (fNumPrimaryDaughters < fNumDaughters)
  {
    if (particleMap.at(pfp->Parent())->PdgCode() == 13)
    {
      fIsTrackDaughter = true;
    }
    if (pfp->NumDaughters())
    {
      for (const int daughter_id : pfp->Daughters())
      {
        if (particleMap.at(daughter_id)->PdgCode() == 13)
        {
          fHasShowerDaughter = true;
        }
      }
    }
  }

  
  const lar_pandora::ClusterVector cluster_vec = particlesToClusters.at(pfp);
  std::cout << cluster_vec.size() << std::endl;
  for (const art::Ptr<recob::Cluster> cluster : cluster_vec)
  {
    if (cluster->isValid())
    {
      std::cout << "Cluster with Plane " << cluster->View() << " has " << cluster->NHits() << " hits." << std::endl;
      if (cluster->View() == 2)
        fNhitsY += cluster->NHits();
      else if (cluster->View() == 0)
        fNhitsU += cluster->NHits();
      else
        fNhitsV += cluster->NHits();
    }
  }
  fNu_NhitsU += fNhitsU;
  fNu_NhitsV += fNhitsV;
  fNu_NhitsY += fNhitsY;
  

  /* One-to-One mapping between hits and sps, under investigation
  for(const art::Ptr<recob::SpacePoint> sps : particlesToSpacePoints.at(pfp) ){
    fNhitsSpacepoints+= spacePointsToHits.at(sps).size();
  }
  fNu_NhitsSpacepoints+=fNhitsSpacepoints;
  std::cout << "Spacepoints have " << fNhitsSpacepoints << " hits." << std::endl; 
  */
 
  if ( particlesToVertices.find(pfp) == particlesToVertices.end() )
  {
    // If a daughter has no associated vertex, count the hits to contribute to the total, but dont save the daughter
    std::cout << "[NueCC::FillReconstructed] Daughter had no associated vertex." << std::endl;
    return false;
  }
  if (particlesToMetadata.at(pfp).size() != 1 || particlesToVertices.at(pfp).size() != 1)
  {
    std::cout << "[NueCC::FillReconstructed] Daughter association problem." << std::endl;
    return false;
  }

  const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties = particlesToMetadata.at(pfp).front()->GetPropertiesMap();
  fTrackScore = pfp_properties.at("TrackScore");

  const recob::Vertex::Point_t &pfp_vtx = particlesToVertices.at(pfp).front()->position();
  fVx = pfp_vtx.X();
  fVy = pfp_vtx.Y();
  fVz = pfp_vtx.Z();
  std::cout << "Daughter with Trackscore: " << fTrackScore << ", Generation: " << fGeneration << ", Vz: " << fVz << std::endl;
  fNueDaughtersTree->Fill();
  return true;
}

void NueCC::FillTruth(art::Event const &e)
{
  std::cout << "Filling Truth info..." << std::endl;
}