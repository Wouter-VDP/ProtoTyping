#ifndef MCPARTICLEHELPER_H
#define MCPARTICLEHELPER_H

#include "larcorealg/Geometry/TPCGeo.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcore/Geometry/Geometry.h"

#include "GeometryHelper.h"
#include "TVector3.h"

struct MCParticleInfo
{
    uint process; // std::string
    //bool pfp_matched; // if the mc particle is primary, is it matched with a pfp?
    bool startInside;
    bool endInside;
    bool partInside;  // This means that the track is crossing, starts/ends inside, is completely inside.
    float startX_tpc; // if the track starts outside but crosses, this will be stored here.
    float startY_tpc;
    float startZ_tpc;
    float startX_sce; // spacecharge corrected version of the tpc edge or the inside start end point.
    float startY_sce;
    float startZ_sce;
    float endX_tpc;
    float endY_tpc;
    float endZ_tpc;
    float endX_sce;
    float endY_sce;
    float endZ_sce;
    float length;
    float lengthTPC;
    float length_sce;
};

class MCParticleHelper
{
  public:
    MCParticleHelper();
    ~MCParticleHelper() = default;

    std::set<std::string> getProcesses() { return string_process; };
    MCParticleInfo fillMCP(simb::MCParticle const &mcparticle);

  private:
    geo::BoxBoundedGeo theTpcGeo;
    art::ServiceHandle<geo::Geometry> geo_service;
    GeometryHelper geoHelper;

    std::set<std::string> string_process; // This variable counts the different processes invloved.
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
            {"primary", 23}};
};

MCParticleHelper::MCParticleHelper()
{
    //// Set up the boxboundedgeo:
    geo::TPCGeo const &thisTPC = geo_service->TPC();
    theTpcGeo = thisTPC.ActiveBoundingBox();
    std::cout << std::endl;
    std::cout << "[MCParticleHelper constructor] ActivePTC volume max: " << theTpcGeo.MaxX() << ", " << theTpcGeo.MaxY() << ", " << theTpcGeo.MaxZ() << std::endl;
    std::cout << "[MCParticleHelper constructor] ActivePTC volume min: " << theTpcGeo.MinX() << ", " << theTpcGeo.MinY() << ", " << theTpcGeo.MinZ() << std::endl;
    //// Check if spacecharge correction is working
    auto const &SCE(*lar::providerFrom<spacecharge::SpaceChargeService>());
    auto scecorr = SCE.GetPosOffsets(geo::Point_t(25, 110, 250));
    double xOffset = scecorr.X();
    double yOffset = scecorr.Y();
    double zOffset = scecorr.Z();
    std::cout << "[CosmicStudies constructor] Spacecharge correction at test point " << xOffset << ", " << yOffset << ", " << zOffset << std::endl;
    std::cout << std::endl;
}

struct MCParticleInfo MCParticleHelper::fillMCP(simb::MCParticle const &mcparticle)
{
    MCParticleInfo this_mcp;

    string_process.insert(mcparticle.Process());
    if (map_process.find(mcparticle.Process()) != map_process.end())
    {
        this_mcp.process = map_process[mcparticle.Process()];
    }
    else
    {
        this_mcp.process = 999;
        std::cout << "[CosmicStudies::analyze] New MC interaction process found!" << std::endl;
    }

    TVector3 mc_start = mcparticle.Position().Vect();
    TVector3 mc_end = mcparticle.EndPosition().Vect();

    std::vector<double> start = {mc_start.X(), mc_start.Y(), mc_start.Z()};
    std::vector<double> end = {mc_end.X(), mc_end.Y(), mc_end.Z()};
    this_mcp.startInside = geoHelper.isActive(start);
    this_mcp.endInside = geoHelper.isActive(end);
    this_mcp.length = geoHelper.distance(start, end); // This is the total length, not the length in the detector!

    this_mcp.partInside=true;
    //Find the section that is inside the tpc
    if (!this_mcp.startInside || !this_mcp.endInside)
    {
        TVector3 mc_startdir = mcparticle.Momentum().Vect();
        std::vector<TVector3> intersections = theTpcGeo.GetIntersections(mc_start, mc_startdir);
        uint num_intersections = intersections.size();

        //Particles completely passes the TPC without entering
        if (num_intersections == 0)
        {
            this_mcp.partInside = false;
            this_mcp.lengthTPC = 0;
            this_mcp.length_sce = 0;
        }

        // Particle started inside the TPC and is exiting
        else if (num_intersections == 1)
        {
            this_mcp.startX_tpc = mc_start.X();
            this_mcp.startY_tpc = mc_start.Y();
            this_mcp.startZ_tpc = mc_start.Z();
            this_mcp.endX_tpc = intersections[0].X();
            this_mcp.endY_tpc = intersections[0].Y();
            this_mcp.endZ_tpc = intersections[0].Z();
        }

        // Particle crosses TPC or particle stops in TPC or particle stops before TPC
        else if (num_intersections == 2)
        {
            float len_start_cross0 = (mc_start - intersections[0]).Mag();
            float len_start_cross1 = (mc_start - intersections[1]).Mag();

            // Particle is crossing the TPC
            if (std::max(len_start_cross0, len_start_cross1) < this_mcp.length)
            {
                this_mcp.startX_tpc = intersections[0].X();
                this_mcp.startY_tpc = intersections[0].Y();
                this_mcp.startZ_tpc = intersections[0].Z();
                this_mcp.endX_tpc = intersections[1].X();
                this_mcp.endY_tpc = intersections[1].Y();
                this_mcp.endZ_tpc = intersections[1].Z();
            }
            // Particle stops before entering
            else if (std::min(len_start_cross0, len_start_cross1) > this_mcp.length)
            {
                this_mcp.partInside = false;
                this_mcp.lengthTPC = 0;
                this_mcp.length_sce = 0;
            }
            // Particle stops inside
            else
            {
                this_mcp.endX_tpc = mc_end.X();
                this_mcp.endY_tpc = mc_end.Y();
                this_mcp.endZ_tpc = mc_end.Z();
                if (len_start_cross0 < this_mcp.length)
                {
                    this_mcp.startX_tpc = intersections[0].X();
                    this_mcp.startY_tpc = intersections[0].Y();
                    this_mcp.startZ_tpc = intersections[0].Z();
                }
                else
                {
                    this_mcp.startX_tpc = intersections[1].X();
                    this_mcp.startY_tpc = intersections[1].Y();
                    this_mcp.startZ_tpc = intersections[1].Z();
                }
            }
        }
    }
    else //Start and end are inside
    {
        this_mcp.startX_tpc = mc_start.X();
        this_mcp.startY_tpc = mc_start.Y();
        this_mcp.startZ_tpc = mc_start.Z();
        this_mcp.endX_tpc = mc_end.X();
        this_mcp.endY_tpc = mc_end.Y();
        this_mcp.endZ_tpc = mc_end.Z();
    }

    if (this_mcp.partInside)
    {
        std::vector<float> start_tpc = {this_mcp.startX_tpc, this_mcp.startY_tpc, this_mcp.startZ_tpc};
        std::vector<float> end_tpc = {this_mcp.endX_tpc, this_mcp.endY_tpc, this_mcp.endZ_tpc};
        this_mcp.lengthTPC = geoHelper.distance(start_tpc, end_tpc);

        //Correct the inside tpcpoints for spacecharge
        auto const &SCE(*lar::providerFrom<spacecharge::SpaceChargeService>());
        auto sce_start = SCE.GetPosOffsets(geo::Point_t(this_mcp.startX_tpc, this_mcp.startY_tpc, this_mcp.startZ_tpc));
        auto sce_end = SCE.GetPosOffsets(geo::Point_t(this_mcp.endX_tpc, this_mcp.endY_tpc, this_mcp.endZ_tpc));

        this_mcp.startX_sce = this_mcp.startX_tpc - sce_start.X() + 0.7;
        this_mcp.startY_sce = this_mcp.startY_tpc + sce_start.Y();
        this_mcp.startZ_sce = this_mcp.startZ_tpc + sce_start.Z();
        this_mcp.endX_sce = this_mcp.endX_tpc - sce_end.X() + 0.7;
        this_mcp.endY_sce = this_mcp.endY_tpc + sce_end.Y();
        this_mcp.endZ_sce = this_mcp.endZ_tpc + sce_end.Z();

        std::vector<float> start_sce = {this_mcp.startX_sce, this_mcp.startY_sce, this_mcp.startZ_tpc};
        std::vector<float> end_sce = {this_mcp.endX_sce, this_mcp.endY_sce, this_mcp.endZ_sce};
        this_mcp.length_sce = geoHelper.distance(start_sce, end_sce);
    }
    return this_mcp;
}

#endif