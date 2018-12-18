#ifndef PANDORAINTERFACEHELPER_H
#define PANDORAINTERFACEHELPER_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

namespace lar_pandora
{
typedef std::map<art::Ptr<recob::PFParticle>, art::Ptr<simb::MCParticle>> PFParticlesToMCParticles;
}

typedef std::map<art::Ptr<recob::PFParticle>, unsigned int> RecoParticleToNMatchedHits;
typedef std::map<art::Ptr<simb::MCParticle>, RecoParticleToNMatchedHits> ParticleMatchingMap;
typedef std::set<art::Ptr<recob::PFParticle>> PFParticleSet;
typedef std::set<art::Ptr<simb::MCParticle>> MCParticleSet;

class PandoraInterfaceHelper
{
public:
  PandoraInterfaceHelper();
  ~PandoraInterfaceHelper() = default;

  void reconfigure(fhicl::ParameterSet const &pset);

  void get_daughter_tracks(std::vector<size_t> pf_ids, const art::Event &evt,
                           std::vector<art::Ptr<recob::Track>> &tracks,
                           std::string m_pfp_producer);

  void get_daughter_showers(std::vector<size_t> pf_ids, const art::Event &evt,
                            std::vector<art::Ptr<recob::Shower>> &showers,
                            std::string m_pfp_producer);

  /**
    *  @brief Returns matching between true and reconstructed particles
    *
    *  @param matchedParticles the output matches between reconstructed and true particles
    */
  void GetRecoToTrueMatches(lar_pandora::PFParticlesToMCParticles &matchedParticles);

  /**
     *  @brief Configure function parameters (call this function first)
     *
     *  @param e the art::Event
     *  @param m_pfp_producer the PFParticle producer label
     *  @param m_spacepoint_producer the SpacePoint producer label
     *  @param m_hitfinder_producer the Hit producer label
     *  @param m_geant_producer The Geant4 producer label
     */
  void Configure(art::Event const &e,
                 std::string m_pfp_producer,
                 std::string m_spacepoint_producer,
                 std::string m_hitfinder_producer,
                 std::string m_geant_producer,
                 std::string m_hit_mcp_producer);

  void Configure(art::Event const &e,
                 std::string m_pfp_producer,
                 std::string m_spacepoint_producer,
                 std::string m_hitfinder_producer,
                 std::string m_geant_producer);

  art::Ptr<simb::MCTruth> TrackIDToMCTruth(art::Event const &e, std::string _geant_producer, int geant_track_id);
  /**
     *  @brief Collect a vector of MCParticle objects from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param particleVector the output vector of MCParticle objects
     */
  void CollectMCParticles(const art::Event &evt,
                          const std::string &label,
                          lar_pandora::MCParticleVector &particleVector);
  /**
     *  @brief Collect truth information from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param truthToParticles output map from MCTruth to MCParticle objects
     *  @param particlesToTruth output map from MCParticle to MCTruth objects
     */
  void CollectMCParticles(const art::Event &evt,
                          const std::string &label,
                          lar_pandora::MCTruthToMCParticles &truthToParticles,
                          lar_pandora::MCParticlesToMCTruth &particlesToTruth);

  /**
 *  @brief  Collect all downstream particles of those in the input vector
 *
 *  @param  pfParticleMap the mapping from PFParticle ID to PFParticle
 *  @param  parentPFParticles the input vector of PFParticles
 *  @param  downstreamPFParticle the output vector of PFParticles including those downstream of the input
 */
  void CollectDownstreamPFParticles(const lar_pandora::PFParticleMap &pfParticleMap,
                                    const art::Ptr<recob::PFParticle> &particle,
                                    lar_pandora::PFParticleVector &downstreamPFParticles) const;

protected:
  lar_pandora::HitsToMCParticles m_hit_to_mcps_map; ///< A map from recon hits to MCParticles
  lar_pandora::PFParticlesToHits m_pfp_to_hits_map; ///< A map from PFParticles to recon hits
private:
  bool m_configured = false;
  bool m_isOverlaidSample = false;
  bool m_debug = false;
  bool m_verbose = false;
};

#endif