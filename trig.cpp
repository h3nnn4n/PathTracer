#include "trig.h"
#include "vector.h"
#include "path.h"
#include "color.h"

#define EPSILON 0.000001

double Triangle::intersect(Path p){
        Vector O = p.begin;
        Vector D = p.end;

        Vector edge1, edge2;
        Vector P, Q, T;
        float det, inv_det, u, v;
        float t;

        //Find vectors for two edges sharing v1
        edge1 = v2 - v1;
        edge2 = v3 - v1;

        //Begin calculating determinant - also used to calculate u parameter
        P = D.cross(edge2);
        //if determinant is near zero, ray lies in plane of triangle
        det = edge1.dot(P);

        //NOT CULLING
        if(det > -EPSILON && det < EPSILON) return 0;
        inv_det = 1.f / det;

        //calculate distance from v1 to ray origin
        T =  O - v1;

        //Calculate u parameter and test bound
        u = T.dot(P) * inv_det;
        //The intersection lies outside of the triangle
        if(u < 0.f || u > 1.f) return 0;

        //Prepare to test v parameter
        Q =  T.cross(edge1);

        //Calculate V parameter and test bound
        v = D.dot(Q) * inv_det;
        //The intersection lies outside of the triangle
        if(v < 0.f || u + v  > 1.f) return 0;

        t = edge2.dot(Q) * inv_det;

        if(t > EPSILON) { //ray intersection
            return t;
        }

        // No hit, no win
    return 0.0;
}
