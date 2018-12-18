#include "NueCC.h"

void NueCC::analyze(art::Event const &evt)
{
  clear();
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
  auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);
  if (pfparticle_handle->size() == 0)
  {
    std::cout << "[NueCC::FillReconstructed] No reconstructed PFParticles in event." << std::endl;
  }
  else
  {
    auto const &pfparticles(*pfparticle_handle);
    for (auto const &pfp : pfparticles)
    {
      if (pfp.PdgCode() == 12 || pfp.PdgCode() == 14)
      {
        // Reconstructed neutrino info
        fNu_PDG = pfp.PdgCode();
        fNumDaughters = pfp.NumDaughters();
        std::cout << "[NueCC::FillReconstructed] neutrino found! PDG: " << pfp.PdgCode() << " Primary Daughters: " << fNumDaughters << std::endl;

        // Get the map PFP->MCP and the set of MCPs
        if (!m_isData)
        {
          pandoraHelper.Configure(evt, m_pfp_producer, m_pfp_producer, m_hitfinder_producer, m_geant_producer, m_hit_mcp_producer);
          pandoraHelper.GetRecoToTrueMatches(matchedParticles);
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

        // Daughters of the neutrino
        std::vector<size_t> unordered_daugthers;
        std::cout << "[NueCC::FillReconstructed] Before  traversePFParticleTree, pfp.Self(): " << pfp.Self() << std::endl;
        pandoraHelper.traversePFParticleTree(pfparticle_handle, pfp.Self(), unordered_daugthers, m_pfp_producer);
        std::cout << "[NueCC::FillReconstructed] All Daughters: " << unordered_daugthers.size() << std::endl;
        for (const size_t pfp_id : unordered_daugthers)
        {
          FillDaughters(pfparticles.at(pfp_id));
        }
      }
    }
  }
}

bool NueCC::FillDaughters(recob::PFParticle const &pfp)
{
  std::cout << "Daughter with PDG: " << pfp.PdgCode() << std::endl;
  return false;
}
