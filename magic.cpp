#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

#include "utils.h"
#include "vector.h"
#include "sphere.h"
#include "color.h"
#include "plane.h"
#include "path.h"
#include "trig.h"

#ifndef MAX_BOUNCE
#define MAX_BOUNCE 16
#endif /* MAX_BOUNCE */

Plane planes[] = {
};

//Sphere spheres[] = {//Scene: radius, position, emission, color, material
  //Sphere(1e5, Vector( 1e5+1,40.8,81.6), Color()         ,Color(.75,.25,.25),DIFF),//Left
  //Sphere(1e5, Vector(-1e5+99,40.8,81.6),Color()         ,Color(.25,.25,.75),DIFF),//Rght
  //Sphere(1e5, Vector(50,40.8, 1e5),     Color()         ,Color(.75,.75,.75),DIFF),//Back
  //Sphere(1e5, Vector(50,40.8,-1e5+170), Color()         ,Color(),           DIFF),//Frnt
  //Sphere(1e5, Vector(50, 1e5, 81.6),    Color()         ,Color(.75,.75,.75),DIFF),//Botm
  //Sphere(1e5, Vector(50,-1e5+81.6,81.6),Color()         ,Color(.75,.75,.75),DIFF),//Top
  //Sphere(16.5,Vector(27,16.5,47),       Color()         ,Color(0,1,1)*.999, DIFF),//Mirr
////Sphere(16.5,Vector(73,16.5,78),       Color()         ,Color(1,1,1)*.999, SPEC),//Glas
  //Sphere(600, Vector(50,681.6-.27,81.6),Color(12,12,12) ,Color(          ), DIFF) //Lite
//};
//int w_octa = 20;
//int x_octa = 73;
//int y_octa = 20;
//int z_octa = 78;

//Triangle triangles[] = {
    //Triangle(Vector( -w_octa + x_octa,  0 + y_octa,  0 + z_octa), Vector( 0 + x_octa,  w_octa + y_octa,  0 + z_octa), Vector( 0 + x_octa,  0 + y_octa,  w_octa + z_octa), Color(0, 0, 0), Color(.25, .75, .25) * .99, SPEC),
    //Triangle(Vector(  w_octa + x_octa,  0 + y_octa,  0 + z_octa), Vector( 0 + x_octa,  w_octa + y_octa,  0 + z_octa), Vector( 0 + x_octa,  0 + y_octa, -w_octa + z_octa), Color(0, 0, 0), Color(.25, .75, .25) * .99, SPEC),
    //Triangle(Vector( -w_octa + x_octa,  0 + y_octa,  0 + z_octa), Vector( 0 + x_octa,  w_octa + y_octa,  0 + z_octa), Vector( 0 + x_octa,  0 + y_octa, -w_octa + z_octa), Color(0, 0, 0), Color(.25, .75, .25) * .99, SPEC),
    //Triangle(Vector(  w_octa + x_octa,  0 + y_octa,  0 + z_octa), Vector( 0 + x_octa, -w_octa + y_octa,  0 + z_octa), Vector( 0 + x_octa,  0 + y_octa,  w_octa + z_octa), Color(0, 0, 0), Color(.25, .75, .25) * .99, SPEC),
    //Triangle(Vector( -w_octa + x_octa,  0 + y_octa,  0 + z_octa), Vector( 0 + x_octa, -w_octa + y_octa,  0 + z_octa), Vector( 0 + x_octa,  0 + y_octa,  w_octa + z_octa), Color(0, 0, 0), Color(.25, .75, .25) * .99, SPEC),
    //Triangle(Vector(  w_octa + x_octa,  0 + y_octa,  0 + z_octa), Vector( 0 + x_octa, -w_octa + y_octa,  0 + z_octa), Vector( 0 + x_octa,  0 + y_octa, -w_octa + z_octa), Color(0, 0, 0), Color(.25, .75, .25) * .99, SPEC),
    //Triangle(Vector( -w_octa + x_octa,  0 + y_octa,  0 + z_octa), Vector( 0 + x_octa, -w_octa + y_octa,  0 + z_octa), Vector( 0 + x_octa,  0 + y_octa, -w_octa + z_octa), Color(0, 0, 0), Color(.25, .75, .25) * .99, SPEC)
//};

double light_x = 69;
double light_y = 16.5;
double light_z = -30;

Sphere spheres[] = {
    Sphere(1e5 , Vector(50, 1e5, 81.6),             Color(.0, .0, .0), Color(.99, .99, .99),    DIFF), // Floor
    Sphere(16.5, Vector(27, 16.5, 47),              Color(.0, .0, .0), Color(.00 ,.75, .99),    SPEC), // Mirror ball
    Sphere(16.5, Vector(73, 16.5, 78),              Color(.0, .0, .0), Color(.15, .15, .75),    DIFF), // Difuse ball front
    Sphere(16.5, Vector(113,16.5,-10),              Color(.0, .0, .0), Color(.95, .0,   .0),    DIFF), // Difuse ball behind

    Sphere(16.5 + 80, Vector(133, 16.5 + 80, -100), Color(.0, .0, .0), Color(.99, .99, .99),    SPEC), // Big mirror ball

    Sphere(16.5, Vector(light_x,
                        light_y,
                        light_z),            Color(9., 9., 9.), Color(.99 , .99 , .99 ),        DIFF) // Light ball
};

int x  = 40;
int xx = 40;

Triangle triangles[] = {
    Triangle(Vector(  x,  0 + xx,  0), Vector( 0,  x + xx,  0), Vector( 0,  0 + xx,  x), Color(0, 0, 0), Color(1, 1, 1) * .99, DIFF),
    Triangle(Vector( -x,  0 + xx,  0), Vector( 0,  x + xx,  0), Vector( 0,  0 + xx,  x), Color(0, 0, 0), Color(1, 1, 1) * .99, DIFF),
    Triangle(Vector(  x,  0 + xx,  0), Vector( 0,  x + xx,  0), Vector( 0,  0 + xx, -x), Color(0, 0, 0), Color(1, 1, 1) * .99, DIFF),
    Triangle(Vector( -x,  0 + xx,  0), Vector( 0,  x + xx,  0), Vector( 0,  0 + xx, -x), Color(0, 0, 0), Color(1, 1, 1) * .99, DIFF),

    Triangle(Vector(  x,  0 + xx,  0), Vector( 0, -x + xx,  0), Vector( 0,  0 + xx,  x), Color(0, 0, 0), Color(1, 1, 1) * .99, DIFF),
    Triangle(Vector( -x,  0 + xx,  0), Vector( 0, -x + xx,  0), Vector( 0,  0 + xx,  x), Color(0, 0, 0), Color(1, 1, 1) * .99, DIFF),
    Triangle(Vector(  x,  0 + xx,  0), Vector( 0, -x + xx,  0), Vector( 0,  0 + xx, -x), Color(0, 0, 0), Color(1, 1, 1) * .99, DIFF),
    Triangle(Vector( -x,  0 + xx,  0), Vector( 0, -x + xx,  0), Vector( 0,  0 + xx, -x), Color(0, 0, 0), Color(1, 1, 1) * .99, DIFF)
};

bool intersects_trig(Path p, int *n, double *dist){
    double size = sizeof(triangles) / sizeof(Triangle);
    double dold = 1<<20;
          *dist = 1<<20;

    for(int i = 0 ; i < size ; i++){
        dold = triangles[i].intersect(p);
        if (dold > 0 && dold < *dist){
            *n = i;
            *dist = dold;
        }
    }

    return *dist < 1<<20;
}

bool intersects(Path p, int *n, double *dist){
    double size = sizeof(spheres) / sizeof(Sphere);
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

    b = intersects_trig(ray, &id_plane, &distance_plane);

    if ( !a && !b){
        return Color();
    }

    if (distance_plane < distance){
        distance = distance_plane;                                           // Distance betway where the ray was cast and interception point
        id = id_plane;                                                       // Position in the array where the object is

        Triangle *target = &triangles[id];                                   // Intercepted triangle

        Vector x  = ray.begin + ray.end * distance;                          // Intersection point
        Vector n  = (x - ((target->v2-x).cross(target->v3-x))).normalized(); // Normal at reflection point
        Color  f  = target->col;                                             // Color to be reflected
        Vector nl;                                                           // This will be the 'absolute value'

        if ( n.dot(ray.end) < 0 ){                                           // If the dot product is negative, turns the ray around
            nl = n;                                                          //
        } else {                                                             //
            nl = n * -1.0;                                                   //
        }                                                                    //

        double p = 0.0;
        
        if ( f.r > p ){                                                      // Get the greater color intersity
            p = f.r;                                                         // 
        } else if ( f.g > p ){                                               // 
            p = f.g;                                                         // 
        } else if ( f.b > p ){                                               // 
            p = f.b;                                                         // 
        }                                                                    // 

        if (iter > MAX_BOUNCE){                                              // If ray was reflected MAX_BOUNCE times stop and return face emission
            return target->emission;
        } else if (++iter > 5){                                              // If more than 5 reflections were done star russian roulette
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
            return target->emission + f * (tracer(Path(x, ray.end - n * 2.0 * n.dot(ray.end)), iter));    // Ideal reflection
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

        if (iter > MAX_BOUNCE){
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

    return Color();
}

int render(int i){
    int screenx = 200 * 2;
    int screeny = 150 * 2;

    Path cam(Vector(50,  52      , 295.6),
             Vector( 0, -0.042612, -1).normalized());

    Vector dx = Vector(screenx * .5135 / screeny);
    Vector dy = (dx.cross(cam.end)).normalized() * .5135;

    Color r;

    Color *image = new Color[screenx * screeny];

    int samps = 50;

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

    std::cout << "\rDone   \n";

    char name[256];

    sprintf(name, "out_%03d.ppm", i);

    //FILE *f = fopen("out.ppm", "w");         // PPM cancer
    FILE *f = fopen(name, "w");         // PPM cancer
    fprintf(f, "P3\n%d %d\n%d\n", screenx, screeny, 255);

    for (int i = 0; i < screenx * screeny; i++) {
        //fprintf(f,"%d %d %d ", (image[i].gammaCorrection().r),
                               //(image[i].gammaCorrection().g),
                               //(image[i].gammaCorrection().b));
        fprintf(f,"%d %d %d ", toInt(image[i].r),
                               toInt(image[i].g),
                               toInt(image[i].b));
        }

    fclose(f);

    return 0;
}

int main(){

    int i;

    int ay = 0;

    light_y = 16.5 + 35;

    for ( i = 0; i < 17; i++ ){
        printf("\n %d %f\n", i, light_y);
        ay += 1;
        light_y -= ay;

        if ( light_y == 16.5 ) {
            ay *= -1;
        } else if ( light_y < 16.5 ) {

            spheres[5] =  Sphere(16.5, Vector(light_x,
                                              16.5   ,
                                              light_z ), Color(9., 9., 9.), Color(.99 , .99 , .99 ), DIFF); // Light ball
            render(i);

            ay += 1;
            ay *= -1;
            continue;
        }

        spheres[5] =  Sphere(16.5, Vector(light_x,
                                          light_y,
                                          light_z ), Color(9., 9., 9.), Color(.99 , .99 , .99 ), DIFF); // Light ball
        render(i);
    }

    return 0;
}
