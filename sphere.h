#ifndef __SPHERE_PT
#define __SPHERE_PT

#include <cmath>

#include "path.h"
#include "color.h"

enum Refl_t { DIFF, SPEC, REFR };

class Sphere{
    private:

    public:
        double radius;
        Vector position;
        Color emission;
        Color col;
        Refl_t material;

        Sphere(double r, Vector pos) : radius(r), position(pos) {}
        Sphere(double r, Vector pos, Color e, Color c) : radius(r), position(pos), emission(e), col(c) {}
        Sphere(double r, Vector pos, Color e, Color c, Refl_t m) : radius(r), position(pos), emission(e), col(c), material(m) {}

        double intersect(Path r);
};

#endif /* __SPHERE_PT */
