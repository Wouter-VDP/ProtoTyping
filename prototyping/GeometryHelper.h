#ifndef GEOMETRYHELPER_H
#define GEOMETRYHELPER_H

#include "larcore/Geometry/Geometry.h"

class GeometryHelper
{
  public:
    GeometryHelper() = default;
    ~GeometryHelper() = default;

    bool isActive(const std::vector<double> &x) const;
    bool isActive(const double x[3]) const;

    double distance(const std::vector<double> &a, const std::vector<double> &b) const;

  private:
};

bool GeometryHelper::isActive(const std::vector<double> &x) const
{
    if (x.size() != 3)
    {
        return false;
    }

    return this->isActive(&x[0]);
}

bool GeometryHelper::isActive(const double x[3]) const
{

    art::ServiceHandle<geo::Geometry> geo;
    std::vector<double> bnd = {
        0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(),
        0., geo->DetLength()};

    bool is_x = x[0] > bnd[0] && x[0] < bnd[1];
    bool is_y = x[1] > bnd[2] && x[1] < bnd[3];
    bool is_z = x[2] > bnd[4] && x[2] < bnd[5];
    return is_x && is_y && is_z;
}

double GeometryHelper::distance(const std::vector<double> &a,
                                const std::vector<double> &b) const
{
    if (a.size() != 3 || b.size() != 3)
    {
        return -1;
    }

    double d = 0;

    for (int i = 0; i < 3; i++)
    {
        d += pow((a[i] - b[i]), 2);
    }

    return sqrt(d);
}

#endif