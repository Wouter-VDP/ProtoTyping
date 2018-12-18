#include "NueCC.h"

void NueCC::analyze(art::Event const &evt)
{
  clearEvent();
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  std::cout << "[NueCC::analyze]: Run " << fRun << ", Subrun " << fSubrun << ", Event " << fEvent << std::endl;

  // Fill true neutrino info
  // TODO

  // -- Fill reco neutrino info --
  FillReconstructed(evt);
}

void NueCC::FillReconstructed(art::Event const &evt)
{
  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticlesToMetadata particlesToMetadata;
  larpandora.CollectPFParticleMetadata(evt, m_pfp_producer, pfparticles, particlesToMetadata);
  if (pfparticles.size() == 0)
  {
    std::cout << "[NueCC::FillReconstructed] No reconstructed PFParticles in event." << std::endl;
  }
  else
  {
    lar_pandora::PFParticleVector pfneutrinos;
    larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
    if (pfneutrinos.size() != 1)
    {
      std::cout << "[NueCC::FillReconstructed] Number of reconstructed neutrinos in event is " << pfneutrinos.size() << std::endl;
    }
    else // This means we have a reconstructed neutrino!
    {
      // Load associations and collections
      
      
      // Start filling information
      art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
      fNu_PDG = pfnu->PdgCode();
      fNumPrimaryDaughters = pfnu->NumDaughters();
      lar_pandora::MetadataVector neutrino_metadata_vec = particlesToMetadata[pfnu];
      if (neutrino_metadata_vec.size() != 1)
      {
        std::cout << "[NueCC::FillReconstructed] Neutrino metadata problem." << std::endl;
      }
      else
      {
        const larpandoraobj::PFParticleMetadata::PropertiesMap &neutrino_properties = neutrino_metadata_vec.front()->GetPropertiesMap();
        fNu_Score = neutrino_properties.at("NuScore");
        fNu_SliceIndex = neutrino_properties.at("SliceIndex");
      }

      lar_pandora::PFParticleMap particleMap;
      larpandora.BuildPFParticleMap(pfparticles, particleMap);
      lar_pandora::PFParticleVector pfdaughters;
      pandoraInterfaceHelper.CollectDownstreamPFParticles(particleMap, pfnu, pfdaughters);
      fNumDaughters = pfdaughters.size() - 1; // The neutrino itself is included here.
      std::cout << "[NueCC::FillReconstructed] neutrino found! PDG: " << fNu_PDG << ", Primary Daughters: " << fNumPrimaryDaughters;
      std::cout << ", Daughters: " << fNumDaughters << ", TopoScore: " << fNu_Score << std::endl;
      for (auto const pfp : pfdaughters)
      {
        if (!pfp->IsPrimary())
          if (!FillDaughters(pfp, particleMap, particlesToMetadata[pfp]))
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
  }
}

bool NueCC::FillDaughters(const art::Ptr<recob::PFParticle> &pfp,
                          const lar_pandora::PFParticleMap &particleMap,
                          const lar_pandora::MetadataVector &metadata_vec)
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
  if (metadata_vec.size() != 1)
  {
    std::cout << "[NueCC::FillReconstructed] Neutrino metadata problem." << std::endl;
    return false;
  }
  const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties = metadata_vec.front()->GetPropertiesMap();
  fTrackScore = pfp_properties.at("TrackScore");
  std::cout << "Daughter with Trackscore: " << fTrackScore << ", Generation: " << fGeneration << std::endl;
  return true;
}
