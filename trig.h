#ifndef __TRIG_PT
#define __TRIG_PT

#include "vector.h"
#include "path.h"
#include "color.h"
#include "utils.h"

class Triangle {
    private:

    public:
        Vector v1, v2, v3;
        Color emission;
        Color col;
        Refl_t material;

        Triangle(Vector a1, Vector a2, Vector a3): v1(a1), v2(a2), v3(a3) {}
        Triangle(Vector a1, Vector a2, Vector a3, Color e, Color c): v1(a1), v2(a2), v3(a3), emission(e), col(c) {}
        Triangle(Vector a1, Vector a2, Vector a3, Color e, Color c, Refl_t mm): v1(a1), v2(a2), v3(a3), emission(e), col(c), material(mm) {}

        double intersect(Path p);
};

#endif /* __TRIG_PT */

