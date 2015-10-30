#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "vector.h"
#include "sphere.h"
#include "color.h"
#include "path.h"

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

Sphere spheres[] = {
    Sphere(1e5 , Vector(50, 1e5, 81.6),      Color(.0, .0, .0), Color(.75, .75, .75), DIFF), // Floor

    Sphere(16.5, Vector(27, 16.5, 47),       Color(.0, .0, .0), Color(1. , .75, 1. ), SPEC), // Mirror ball
    Sphere(99.9, Vector(-50, 70,-170),       Color(.0, .0, .0), Color(1. , .99, 1. ), SPEC), // Big Mirror ball

    Sphere(16.5, Vector(73, 16.5, 78),       Color(.0, .0, .0), Color(.15, .15, .75), DIFF), // Difuse ball front

    Sphere(16.5, Vector(113,16.5,-10),       Color(.0, .0, .0), Color(.95, .0,   .0), DIFF), // Difuse ball behind

    Sphere(16.5, Vector(69, 16.5,-30),       Color(1.,1.,1.31), Color(1. , 1. , 1. ), DIFF) // Light ball
};

class Plane{ // Nao funciona :)
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

        double intersect(Path p){
            Vector u = p.getDir();      // B - A
            Vector w = p.begin - position;

            double D =  normal.dot(u);
            double N = -normal.dot(w);

            if ( D == 0 ){
                if ( N == 0 ){
                    std::cout << "Inside" << '\n';
                    return 0.0;
                } else {
                    std::cout << "Parallel" << '\n';
                    return 0.0;
                }
            }

            double intersec = N / D;

            if ( intersec < 0 ){
                return 0.0;
            }

            return intersec;
        }
};

Plane planes[] = {
    //Plane(Vector(0, 0, 0), Vector(69, 16.5, -30), Color(0, 0, 0), Color(0, 0.1, 0), SPEC)
};

bool intersects_plane(Path p, int *n, double *dist){
    double size = sizeof(planes) / sizeof(Plane);
    double distance;
    double dold = 1<<20;
          *dist = 1<<20;

    for(int i = 0 ; i < size ; i++){
        dold = planes[i].intersect(p);
        if (dold > 0 && dold < *dist){
            *n = i;
            *dist = dold;
        }
    }

    return *dist < 1<<20;
}

bool intersects(Path p, int *n, double *dist){
    double size = sizeof(spheres) / sizeof(Sphere);
    double distance;
    double dold = 1<<20;
          *dist = 1<<20;

    for(int i = 0 ; i < size ; i++){
        dold = spheres[i].intersect(p);
        if (dold > 0 && dold < *dist){
            *n = i;
            *dist = dold;
        }
    }

    return *dist < 1<<20;
}

Color tracer(Path ray, int iter){
    double distance_plane;
    double distance;
    int id = 0;
    int id_plane = 0;

    bool a, b;

    a = intersects(ray, &id, &distance);

    b = intersects_plane(ray, &id_plane, &distance_plane);

    if ( !a && !b){
        return Color();
    }

    if (distance_plane < distance){ // Nao Funciona
        distance = distance_plane;
        id = id_plane;

        Plane *target = &planes[id];

        Vector x  = ray.begin + ray.end * distance;
        Vector n  = (x - target->normal).normalized();
        Vector nl = n.dot(ray.begin) < 0 ? n : n * -1;
        Color  f  = target->col;

        if ( n.dot(ray.end) < 0 ){
            nl = n;
        } else {
            nl = n * -1.0;
        }

        double p = f.r > f.g && f.r > f.b ? f.r : f.g > f.b ? f.g : f.b;

        if (iter > 2){
            return target->emission;
        } else if (++iter > 5){
            if ( drand48() < p){
                f = f * (1 / p);
            } else {
                return target->emission;
            }
        }

        if (target->material == SPEC){
                return target->emission + f * (tracer(Path(x, ray.end - n * 2.0 * n.dot(ray.end)), iter));
        }

        return Color();
    } else {
        Sphere *target = &spheres[id];
        Vector x  = ray.begin + ray.end * distance;
        Vector n  = (x - target->position).normalized();
        Vector nl = n.dot(ray.begin) < 0 ? n : n * -1;
        Color  f  = target->col;

        if ( n.dot(ray.end) < 0 ){
            nl = n;
        } else {
            nl = n * -1.0;
        }

        double p = f.r > f.g && f.r > f.b ? f.r : f.g > f.b ? f.g : f.b;

        if (iter > 75){
            return target->emission;
        } else if (++iter > 10){
            if ( drand48() < p){
                f = f * (1 / p);
            } else {
                return target->emission;
            }
        }

        if (target->material == DIFF) {
            double r1  = 2 * M_PI * drand48();
            double r2  = drand48();
            double r2s = sqrt(r2);

            Vector w = nl;
            Vector u = ((fabs(w.x) > .1 ? Vector(0,1) : Vector(1)).cross(w)).normalized();
            Vector v = w.cross(u);
            Vector d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalized();

            return target->emission + f * (tracer(Path(x, d), iter));

        } else if (target->material == SPEC){
                return target->emission + f * (tracer(Path(x, ray.end - n * 2.0 * n.dot(ray.end)), iter));
        }
    }
}

int main(){
    int screenx = 200 * 2;
    int screeny = 150 * 2;

    Path cam(Vector(50,  52      , 295.6),
             Vector( 0, -0.042612, -1).normalized());

    Vector dx = Vector(screenx * .5135 / screeny);
    Vector dy = (dx.cross(cam.end)).normalized() * .5135;

    Color r;

    Color *image = new Color[screenx * screeny];
    int i;

    int samps = 1000;

#pragma omp parallel for schedule(dynamic, 1) private(r)

    for (int y = 0; y < screeny; y++){
        if ( y %5 == 0){fprintf(stdout, "\r%.2f%%", ((double)y/screeny)*100.0);fflush(stdout);}
        for (int x = 0; x < screenx; x++){
            r = Color();
            for (int sy = 0, i = (screeny - y - 1) * screenx + x; sy < 2; sy++ ){
                for (int sx = 0; sx < 2; sx++){
                    for (int s = 0; s < samps; s++ ){
                        double r1 = 2 * drand48(), ddx = r1<1 ? sqrt(r1) - 1: 1 - sqrt(2 - r1);
                        double r2 = 2 * drand48(), ddy = r2<1 ? sqrt(r2) - 1: 1 - sqrt(2 - r2);
                        Vector d  = cam.end + dx * (((sx + .5  +  ddx)/2 + x)/screenx  -  .5)  +
                                              dy * (((sy + .5  +  ddy)/2 + y)/screeny  -  .5)  ;
                        r  =  r + tracer(Path(cam.begin + d * 140, d.normalized()), 0) * (1. / samps);
                    }
                    image[i] = image[i] + r.truncated() * .25;
                }
            }
        }
    }

    std::cout << "\rDone\n";

    FILE *f = fopen("image.ppm", "w");         // PPM cancer
    fprintf(f, "P3\n%d %d\n%d\n", screenx, screeny, 255);

    for (int i = 0; i < screenx * screeny; i++) {
        //fprintf(f,"%d %d %d ", (image[i].gammaCorrection().r),
                               //(image[i].gammaCorrection().g),
                               //(image[i].gammaCorrection().b));
        fprintf(f,"%d %d %d ", toInt(image[i].r),
                               toInt(image[i].g),
                               toInt(image[i].b));
    }

    return 0;
}
