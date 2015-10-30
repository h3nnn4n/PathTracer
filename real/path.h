#ifndef __PATH_PT
#define __PATH_PT

#include "vector.h"

class Path{
    private:

    public:
        Vector begin, end;

        Path(Vector b, Vector e) : begin(b), end(e) {}

        Vector getPoint(double t);
        Vector getDir();
};

#endif /* __PATH_PT */
