#ifndef __SPHERE_PT
#define __SPHERE_PT

#include <cmath>

#include "path.h"
#include "color.h"
#include "utils.h"

class Sphere{
    private:

    public:
        double radius;
        Vector position;
        Color emission;
        Color col;
        Refl_t material;
        double haziness;

        //Sphere(double r, Vector pos) : radius(r), position(pos) {}
        //Sphere(double r, Vector pos, Color e, Color c) : radius(r), position(pos), emission(e), col(c) {}
        //Sphere(double r, Vector pos, Color e, Color c, double m) : radius(r), position(pos), emission(e), col(c), material(m) {}
        Sphere(double r, Vector pos, Color e, Color c, Refl_t m) : radius(r), position(pos), emission(e), col(c), material(m) {}
        Sphere(double r, Vector pos, Color e, Color c, Refl_t m, double h) : radius(r), position(pos), emission(e), col(c), material(m), haziness(h) {}

        double intersect(Path r);
};

#endif /* __SPHERE_PT */
