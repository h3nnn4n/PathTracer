#ifndef __TRIG_PT
#define __TRIG_PT

#include "vector.h"
#include "path.h"
#include "color.h"

class Triangle {
    private:

    public:
        Vector v1, v2, v3;
        Color emision;
        Color col;

        Triangle(Vector a1, Vector a2, Vector a3): v1(a1), v2(a2), v3(a3) {}
        Triangle(Vector a1, Vector a2, Vector a3, Color e, Color c): v1(a1), v2(a2), v3(a3), emision(e), col(c) {}

        double intersects(Path p);
};

#endif /* __TRIG_PT */

