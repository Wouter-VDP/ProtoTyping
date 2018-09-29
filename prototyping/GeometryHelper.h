#ifndef GEOMETRYHELPER_H
#define GEOMETRYHELPER_H

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"

class GeometryHelper
{
  public:
    GeometryHelper();
    ~GeometryHelper() = default;

    bool isActive(const std::vector<float> &x) const;
    bool isActive(const float x[3]) const;

    float distance(const std::vector<float> &a, const std::vector<float> &b) const;

  private:
};

bool GeometryHelper::isActive(const std::vector<float> &x) const
{
    if (x.size() != 3)
    {
        return false;
    }

    return this->isActive(&x[0]);
}

GeometryHelper::GeometryHelper()
{
    //// Check if things are set up properly:
    std::cout << std::endl;
    std::cout << "[GeometryHelper constructor] Checking set-up" << std::endl;
    art::ServiceHandle<geo::Geometry> geo;
    std::cout << "[GeometryHelper constructor] Detector dimensions from geo: "
              << 2.0 * geo->DetHalfWidth() << ", " << geo->DetHalfHeight() << ", " << geo->DetLength() << std::endl;
}

bool GeometryHelper::isActive(const float x[3]) const
{

    art::ServiceHandle<geo::Geometry> geo;
    std::vector<double> bnd = {
        0., 2.0 * geo->DetHalfWidth(),
        -1 * geo->DetHalfHeight(), geo->DetHalfHeight(),
        0., geo->DetLength()};

    bool is_x = x[0] > bnd[0] && x[0] < bnd[1];
    bool is_y = x[1] > bnd[2] && x[1] < bnd[3];
    bool is_z = x[2] > bnd[4] && x[2] < bnd[5];
    return is_x && is_y && is_z;
}

float GeometryHelper::distance(const std::vector<float> &a,
                               const std::vector<float> &b) const
{
    if (a.size() != 3 || b.size() != 3)
    {
        return -1;
    }

    float d = 0;

    for (int i = 0; i < 3; i++)
    {
        d += pow((a[i] - b[i]), 2);
    }

    return sqrt(d);
}

#endif