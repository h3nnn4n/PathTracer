#ifndef __PLANE_PT
#define __PLANE_PT

#include "vector.h"
#include "color.h"
#include "path.h"
#include "utils.h"

class Plane{
    private:

    public:
        Vector normal;
        Vector position;
        Color  col;
        Color  emission;
        Refl_t material;

        Plane(Vector p, Vector n) : normal(n), position(p) {}
        Plane(Vector p, Vector n, Color cc) : normal(n), position(p), col(cc) {}
        Plane(Vector p, Vector n, Color cc, Color ee) : normal(n), position(p), col(cc), emission(ee) {}
        Plane(Vector p, Vector n, Color cc, Color ee, Refl_t mm) : normal(n), position(p), col(cc), emission(ee), material(mm) {}

        double intersect(Path p);
};

#endif /* __PLANE_PT */
