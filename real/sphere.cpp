#include "vector.h"
#include "sphere.h"

double Sphere::intersect(path r){
    Vector op = position - r.begin;
    double t;
    double eps = 1e-4;
    double b = op.dot(r.end);
    double det = b * b - op.dot(op) + radius * radius;

    if (det<0){
        return 0;
    } else {
        det = sqrt(det);
    }

    return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
}
