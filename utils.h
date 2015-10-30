#ifndef __UTILS_PT
#define __UTILS_PT

#include <cmath>

enum Refl_t { DIFF, SPEC, REFR };

inline double clamp(double x) {
    double w;

    w = x;

    if ( w < 0.0 ){
        w = 0.0;
    } else if ( w > 1.0 ){
        w = 1.0;
    }

    return w;
}

inline int toInt(double x){
    return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

#endif /* __UTILS_PT */
