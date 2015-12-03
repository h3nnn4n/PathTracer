#ifndef __scene
#define __scene

#include "sphere.h"
#include "trig.h"

Plane planes[] = {

};

/*Sphere spheres[] = { //Scene: radius, position, emission, color, material*/
    /*Sphere(1e5, Vector( 1e5+1,40.8,81.6)   ,Color()         ,Color(.75,.25,.25),  DIFF),       //Left*/
    /*Sphere(1e5, Vector(-1e5+99,40.8,81.6)  ,Color()         ,Color(.25,.25,.75),  DIFF),       //Rght*/
    /*Sphere(1e5, Vector(50,40.8, 1e5)       ,Color()         ,Color(.75,.75,.75),  DIFF),       //Back*/
    /*Sphere(1e5, Vector(50,40.8,-1e5+170)   ,Color()         ,Color(           ),  DIFF),       //Frnt*/
    /*Sphere(1e5, Vector(50, 1e5, 81.6)      ,Color()         ,Color(.75,.75,.75),  DIFF),       //Botm*/
    /*Sphere(1e5, Vector(50,-1e5+81.6,81.6)  ,Color()         ,Color(.75,.75,.75),  DIFF),       //Top*/

    /*Sphere(16.5,Vector(27,16.5,47)         ,Color()         ,Color(.99,.99,.99),  SPEC),       //Mirr*/

    /*Sphere(16.5, Vector(73, 16.5, 78)      ,Color(.0,.0,.0) ,Color(.05,.75,.99),  DIFF),       //Difuse ball front*/

    /*Sphere(600, Vector(50,681.6-.27,81.6)  ,Color(12,12,12) ,Color(            ), DIFF)        //Lite*/
/*};*/

/*int w_octa = 20;*/
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
/*};*/

double light_x = 69;
double light_y = 16.5;
double light_z = -30;

Sphere spheres[] = {
    Sphere(1e5 , Vector(50, 1e5, 81.6),             Color(.0, .0, .0), Color(.99, .99, .99),    DIFF), // Floor
    Sphere(16.5, Vector(27, 16.5, 47),              Color(.0, .0, .0), Color(.00, .75, .99),    SPEC), // Mirror ball
    Sphere(16.5, Vector(73, 16.5, 78),              Color(.0, .0, .0), Color(.15, .15, .75),    DIFF), // Difuse ball front
    Sphere(16.5, Vector(113,16.5,-10),              Color(.0, .0, .0), Color(.95, .0 , .0 ),    DIFF), // Difuse ball behind

    Sphere(16.5 + 80, Vector(133, 16.5 + 80, -100), Color(.0, .0, .0), Color(.99, .99, .99),    SPEC), // Big mirror ball

    Sphere(16.5, Vector(light_x,
                        light_y,
                        light_z),                   Color(9., 9., 9.), Color(.99 ,.99 ,.99),    DIFF) // Light ball
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

#endif /* __scene */
