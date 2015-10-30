#include "vector.h"
#include "color.h"
#include "plane.h"
#include "path.h"

double Plane::intersect(Path p){
    Vector u = p.getDir();      // B - A
    Vector w = p.begin - position;

    double D =  normal.dot(u);
    double N = -normal.dot(w);

    if ( D == 0 ){
        if ( N == 0 ){
            return 0.0;
        } else {
            return 0.0;
        }
    }

    double intersec = N / D;

    if ( intersec < 0 ){
        return 0.0;
    }

    return intersec;
}

