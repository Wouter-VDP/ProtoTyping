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
      if (!m_isData)
      {
        FillReconTruthMatching(evt);
        FillTrueNu(evt);
        FillTrueNuDaughters(evt);
      }
      FillReconstructed(evt);
      fEventTree->Fill();
    }
  }
  std::cout << "\n\n";
}

void NueCC::FillReconstructed(art::Event const &evt)
{
  fNumPfp = pfparticles.size();
  // Load associations and collections
  lar_pandora::VertexVector vertexVector_dummy;
  lar_pandora::PFParticleVector particleVector_dummy;
  lar_pandora::SpacePointVector spacePointVector_dummy;
  larpandora.CollectVertices(evt, m_pfp_producer, vertexVector_dummy, particlesToVertices);
  larpandora.CollectPFParticles(evt, m_pfp_producer, particleVector_dummy, particlesToClusters);
  larpandora.CollectPFParticles(evt, m_pfp_producer, particleVector_dummy, particlesToSpacePoints);
  //lar_pandora::ClusterVector clusterVector_dummy;
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
    if (m_hasMCNeutrino)
    {
      fTrueNu_VtxDistance = pandoraInterfaceHelper.Distance3D(fTrueNu_VxSce, fTrueNu_VySce, fTrueNu_VzSce, fNu_Vx, fNu_Vy, fNu_Vz);
    }
  }
  larpandora.BuildPFParticleMap(pfparticles, particleMap);
  pandoraInterfaceHelper.CollectDownstreamPFParticles(particleMap, pfnu, pfdaughters);
  fNumDaughters = pfdaughters.size() - 1; // The neutrino itself is included here.
  std::cout << "[NueCC::FillReconstructed] neutrino PDG: " << fNu_PDG << ", Primary Daughters: " << fNumPrimaryDaughters;
  std::cout << ", Daughters: " << fNumDaughters << ", TopoScore: " << fNu_Score;
  std::cout << ", TrueNu_VtxDistance: " << fTrueNu_VtxDistance << std::endl;

  for (auto const pfp : pfdaughters)
  {
    if (!pfp->IsPrimary())
    {
      if (!FillDaughters(pfp))
      {
        fDaughtersStored = false;
      }
      else
      {
        if (MatchDaughter(evt, pfp))
          fNumMatchedDaughters++;
      }
    }
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
  for (const art::Ptr<recob::Cluster> cluster : cluster_vec)
  {
    if (cluster->isValid())
    {
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

  if (particlesToSpacePoints.find(pfp) == particlesToSpacePoints.end())
  {
    // If a daughter has no associated spacepoints, count the hits to contribute to the total, but dont save the daughter
    std::cout << "[NueCC::FillDaughters] Daughter had no associated vertex." << std::endl;
    return false;
  }
  fNSpacepoints = particlesToSpacePoints.at(pfp).size();
  fNu_NSpacepoints += fNSpacepoints;

  if (particlesToVertices.find(pfp) == particlesToVertices.end())
  {
    // If a daughter has no associated vertex, count the hits to contribute to the total, but dont save the daughter
    std::cout << "[NueCC::FillDaughters] Daughter had no associated vertex." << std::endl;
    return false;
  }
  if (particlesToMetadata.at(pfp).size() != 1 || particlesToVertices.at(pfp).size() != 1)
  {
    std::cout << "[NueCC::FillDaughters] Daughter association problem." << std::endl;
    return false;
  }

  const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties = particlesToMetadata.at(pfp).front()->GetPropertiesMap();
  fTrackScore = pfp_properties.at("TrackScore");

  const recob::Vertex::Point_t &pfp_vtx = particlesToVertices.at(pfp).front()->position();
  fVx = pfp_vtx.X();
  fVy = pfp_vtx.Y();
  fVz = pfp_vtx.Z();
  fVtxDistance = pandoraInterfaceHelper.Distance3D(fVx, fVy, fVz, fNu_Vx, fNu_Vy, fNu_Vz);
  std::cout << "[NueCC::FillDaughters] Trackscore: " << fTrackScore << ", Generation: " << fGeneration;
  std::cout << ", vtx distance: " << fVtxDistance;
  std::cout << ", Hits: (" << fNhitsU << "," << fNhitsV << "," << fNhitsY << ")" << std::endl;
  fNueDaughtersTree->Fill();
  return true;
}

bool NueCC::MatchDaughter(art::Event const &evt, const art::Ptr<recob::PFParticle> &pfp)
{
  if (m_isData)
    return false;
  art::Ptr<simb::MCParticle> matched_mcp;
  if (fGeneration == 2)
  {
    if (matchedParticles.find(pfp) == matchedParticles.end())
      return false;
    matched_mcp = matchedParticles.at(pfp);
  }
  else if (fGeneration == 3)
  {
    // Generation 3 particle get matched to its parent.
    const auto iter(particleMap.find(pfp->Parent()));
    if (iter == particleMap.end())
      throw cet::exception("NueCC::MatchDaughter") << "Scrambled PFParticle IDs" << std::endl;
    const art::Ptr<recob::PFParticle> &pfp_parent = iter->second;

    if (matchedParticles.find(pfp_parent) == matchedParticles.end())
      return false;
    matched_mcp = matchedParticles.at(pfp_parent);
  }
  else
  {
    std::cout << "[NueCC::MatchDaughter] Generation 4 particle is not matched." << std::endl;
    return false;
  }

  if (m_hasMCNeutrino)
  {
    // Is this MC particle neutrino?
    const art::Ptr<simb::MCTruth> mctruth = pandoraInterfaceHelper.TrackIDToMCTruth(evt, m_geant_producer, matched_mcp->TrackId());
    if (mctruth->Origin() == simb::kBeamNeutrino)
    {
      fMatchedNeutrino = true;
    }
    else
    {
      fMatchedNeutrino = false;
      fCosmicMatched = true;
    }

    fTruePDG = matched_mcp->PdgCode();
    fTrueEnergy = matched_mcp->E();
    fTrueVx = matched_mcp->Vx();
    fTrueVy = matched_mcp->Vy();
    fTrueVz = matched_mcp->Vz();
    pandoraInterfaceHelper.SCE(fTrueVx, fTrueVy, fTrueVz, matched_mcp->T(),
                               fTrueVxSce, fTrueVySce, fTrueVzSce);
    std::cout << "[NueCC::MatchDaughter] Daughter matched with PDG: " << fTruePDG << ", neutrino origin: " << fMatchedNeutrino << std::endl;
  }
  return true;
}

void NueCC::FillTrueNu(art::Event const &evt)
{
  if (m_hasMCNeutrino)
  {
    auto const &generator_handle = evt.getValidHandle<std::vector<simb::MCTruth>>("generator");
    auto const &generator(*generator_handle);
    fNumNu = generator.size();
    std::cout << "[NueCC::FillTrueNu] True neutrinos found: " << fNumNu;
    if (generator.size() > 0)
    {
      if (generator.front().Origin() != simb::kBeamNeutrino)
      {
        std::cout << "[NueCC::FillTrueNu] Origin of generator particle is not kBeamNeutrino." << std::endl;
        return;
      }
      const simb::MCNeutrino &mcnu = generator.front().GetNeutrino();

      fTrueNu_InteractionType = mcnu.InteractionType();
      fTrueNu_CCNC = mcnu.CCNC();
      fTrueNu_PDG = mcnu.Nu().PdgCode();
      fTrueNu_Energy = mcnu.Nu().E();
      fTrueNu_LeptonEnergy = mcnu.Lepton().E();
      fTrueNu_LeptonTheta = mcnu.Theta();
      fTrueNu_Time = mcnu.Nu().T();
      fTrueNu_Vx = mcnu.Nu().Vx();
      fTrueNu_Vy = mcnu.Nu().Vy();
      fTrueNu_Vz = mcnu.Nu().Vz();
      pandoraInterfaceHelper.SCE(fTrueNu_Vx, fTrueNu_Vy, fTrueNu_Vz, fTrueNu_Time,
                                 fTrueNu_VxSce, fTrueNu_VySce, fTrueNu_VzSce);
      std::cout << ", CCNC: " << fTrueNu_CCNC  << ", PDG: " << fTrueNu_PDG <<", E: " << fTrueNu_Energy << std::endl;
    }
  }
}

void NueCC::FillTrueNuDaughters(art::Event const &evt)
{
  lar_pandora::MCParticleVector mcparticles;
  larpandora.CollectMCParticles(evt, m_geant_producer, mcparticles);

  for (auto &mcparticle : mcparticles)
  {
    if (!(mcparticle->Process() == "primary" &&
          mcparticle->T() != 0 &&
          mcparticle->StatusCode() == 1))
      continue;

    const art::Ptr<simb::MCTruth> mc_truth = pandoraInterfaceHelper.TrackIDToMCTruth(evt, m_geant_producer, mcparticle->TrackId());
    if (mc_truth->Origin() == simb::kBeamNeutrino)
    {
      fTrueNu_DaughterE.push_back(mcparticle->E());
      fTrueNu_DaughterPDG.push_back(mcparticle->PdgCode());
      if(matchedMCParticles.find(mcparticle)!=matchedMCParticles.end())
      {
        fTrueNu_DaughterMatched.push_back(true);
      }
      else
      {
        fTrueNu_DaughterMatched.push_back(false);
      }
      std::cout << "[NueCC::FillTrueNuDaughters] << PDG: " << mcparticle->PdgCode() << ", E: " << mcparticle->E() << ", was matched? " << fTrueNu_DaughterMatched.back() << std::endl;
    }
  }
}

void NueCC::FillReconTruthMatching(art::Event const &evt)
{
  pandoraInterfaceHelper.Configure(evt, m_pfp_producer, m_pfp_producer, m_hitfinder_producer, m_geant_producer, m_hit_mcp_producer);
  pandoraInterfaceHelper.GetRecoToTrueMatches(matchedParticles);
  std::cout << "[NueCC::FillReconTruthMatching] ";
  std::cout << "Number of PFPparticles in event: " << pfparticles.size() << std::endl;
  for (auto it = matchedParticles.begin(); it != matchedParticles.end(); ++it)
  {
    matchedMCParticles.insert(it->second);
  }
  std::cout << "[NueCC::FillReconTruthMatching] ";
  std::cout << "PFParticlesToMCParticles constructed: Number of PFPparticles matched: " << matchedParticles.size() << std::endl;
}