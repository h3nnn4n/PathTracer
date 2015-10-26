#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

int gammaCorrection(double);

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

class vector{
    private:

    public:
        double x, y, z;

        vector(double xi = 0, double yi = 0, double zi = 0){
            x = xi;
            y = yi;
            z = zi;
        }

        vector(const vector &u){
            x = u.x;
            y = u.y;
            z = u.z;
        }

        vector operator + (vector v){
            return vector(x + v.x, y + v.y, z + v.z);
        }

        vector operator + (double k){
            return vector(x + k, y + k, z + k);
        }

        vector operator - (double k){
            return vector(x - k, y - k, z - k);
        }

        vector operator - (vector v){
            return vector(x - v.x, y - v.y, z - v.z);
        }

        vector operator / (double k){
            return vector(x / k, y / k, z / k);
        }

        vector operator * (double k){
            return vector(x * k, y * k, z * k);
        }

        vector operator * (vector v){
            return vector(x * v.x, y * v.y, z * v.z);
        }

        double dot(vector u){
            return x * u.x + y * u.y + z * u.z;
        }

        double norm(){
            return sqrt( x * x + y * y + z * z);
        }

        vector copy(){
            return vector(x, y, z);
        }

        void normalize(){
            *this = *this / norm();
            return;
        }

        vector cross (vector u){
            return vector(y * u.z - z * u.y,
                          z * u.x - x * u.z,
                          x * u.y - y * u.x);
        }

        vector normalized(){
            vector w = *this / norm();
            return w.copy();
        }

        void print(){ std::cout << "( " << x << ", " << y << ", " << z << " )\n";
        }
};

class color{
    private:

    public:
        double r, g, b;
        color(double rr = 0, double gg = 0, double bb = 0){
            r = rr;
            g = gg;
            b = bb;
        }

        void print(){
            std::cout << "( " << r << ", " << g << ", " << b << " )\n";
        }

        color operator + (color v){
            return color(r + v.r, g + v.g, b + v.b);
        }

        color norm(){
            return color(gammaCorrection(r), gammaCorrection(g), gammaCorrection(b));
        }

        color operator * (double k){
            return color(r * k, g * k, b * k);
        }

        color operator * (color v){
            return color(r * v.r, g * v.g, b * v.b);
        }
};

class path{
    private:

    public:
        vector begin, end;

        path(vector b, vector e) : begin(b), end(e) {}

        vector getPoint(double t){
            return begin + (end - begin) * t;
        }

        vector getDir(){
            return end - begin;
        }
};

enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

class sphere{
    private:

    public:
        double radius;
        vector position;
        color emission;
        color col;
        Refl_t material;

        sphere(double r, vector pos) : radius(r), position(pos) {}
        sphere(double r, vector pos, color e, color c) : radius(r), position(pos), emission(e), col(c) {}
        sphere(double r, vector pos, color e, color c, Refl_t m) : radius(r), position(pos), emission(e), col(c), material(m) {}

        double intersect(path r){
            vector op = position - r.begin; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
            // Vec op = p        - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
            double t;
          //double t;
            double eps = 1e-4;
          //double eps = 1e-4;
            double b = op.dot(r.end);
          //double b = op.dot(r.d);
            double det = b * b - op.dot(op) + radius * radius;
          //double det = b * b - op.dot(op) + rad * rad;

            if (det<0)
                return 0;
            else det = sqrt(det);

            return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
        }
};

sphere spheres[] = {//Scene: radius, position, emission, color, material
    sphere(1e5, vector( 1e5+1,40.8,81.6), color(),color(.75,.25,.25),DIFF),//Left
    sphere(1e5, vector(-1e5+99,40.8,81.6),color(),color(.25,.25,.75),DIFF),//Rght
    sphere(1e5, vector(50,40.8, 1e5),     color(),color(.75,.75,.75),DIFF),//Back
    sphere(1e5, vector(50,40.8,-1e5+170), color(),color(),           DIFF),//Frnt
    sphere(1e5, vector(50, 1e5, 81.6),    color(),color(.75,.75,.75),DIFF),//Botm
    sphere(1e5, vector(50,-1e5+81.6,81.6),color(),color(.75,.75,.75),DIFF),//Top

    //sphere(16.5,vector(27,16.5,47),       color(),color(1,1,1)*.999, SPEC),//Mirr

    //sphere(16.5,vector(73,46.5,68),       vector(),vector(1,1,1)*.999, REFR),//Glas
    //sphere(16.5,vector(53,46.5,68),       vector(),vector(1,1,1)*.999, REFR),//Glas

    //sphere(16.5,vector(63,41.5,73),       vector(),vector(1,1,1)*.999, REFR),//Glas

    //sphere(16.5,vector(73,16.5,78),       color(),color(1,1,1)*.999, SPEC),//Glas
    //sphere(16.5,vector(53,36.5,78),       color(),color(1,1,1)*.999, REFR),//Glas

    sphere(600, vector(50,681.6-.27,81.6),color(12,12,12),  color(), DIFF) //Lite
};

class plane{
    private:

    public:
        vector normal;
        vector point;
        color  c;
        color  emission;

        plane(vector p, vector n) : normal(n), point(p) {}
        plane(vector p, vector n, color cc) : normal(n), point(p), c(cc) {}
        plane(vector p, vector n, color cc, color ee) : normal(n), point(p), c(cc), emission(ee) {}

        double intersect(path p){
            vector u = p.getDir();      // B - A
            vector w = p.begin - point;

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
                //std::cout<<"aksadjkas\n";
                return 0.0;
            }

            return intersec;
        }
};

plane planes[] = {
    plane( vector( 0, 0, 0), vector(1, 0, 0), color(.9, .0, .0), color(0.329, 0.249, 0.959)),
    plane( vector( 0, 0, 0), vector(0, 1, 0), color(.0, .9, .0), color(0.529, 0.429, 0.429)),
    plane( vector( 0, 0, 0), vector(0, 0, 1), color(.0, .0, .9), color(0.299, 0.949, 0.209))
};


bool intersects(path p, int *n, double *dist){
    double size = sizeof(spheres) / sizeof(sphere);
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

    *dist = dold;

    return *dist < 1<<20;
}

double truncate(double w){
    if ( w < 0.0 ){
        return 0.0;
    } else if (w > 1.0) {
        return 1.0;
    } else {
        return w;
    }
}

int gammaCorrection(double x){
    return (pow( truncate(x), 1.0 / 2.2) * 255.0 + .5 );
}

int toInt(double x){
    return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

color tracer(path ray, int iter){
    double distance;
    int id = 0;

    if (!intersects(ray, &id, &distance)){
        return color();
    }

    sphere *target = &spheres[id];

    //vector x     = ray.begin + ray.end * distance;
    //vector norma = (x - target->point).normalized();
    //vector norma = (target->normal);
    //color f      = target->c;
    //vector nl;

    //if ( norma.dot(ray.end) < 0 ){
        //nl = norma;
    //} else {
        //nl = norma * -1.0;
    //}

    //double p = f.r > f.g && f.r > f.b ? f.r : f.g > f.b ? f.g : f.b; // max refl

    //if (iter-- > 0){
        //if ( drand48() < p ){
            //f = f * (1 / p);
        //} else {
            //return target->emission;
        //}
    //}

    //if ( iter < -15 ){
        //return color();
    //}

    vector x = ray.begin + ray.end * distance;
    vector n = (x - target->position).normalized();
    vector nl = n.dot(ray.begin) < 0 ? n : n * -1;
    color f = target->col;

    double p = f.r > f.g && f.r > f.b ? f.r : f.g > f.b ? f.g : f.b; // mar refl

    if (++iter > 5){ // Quits if done more than 5 recursions
        if ( drand48() < p){
            f = f * (1 / p);
        } else {
            return target->emission; //R.R.
        }
    }

    if (target->material == DIFF) {                  // Ideal DIFFUSE reflection
        double r1  = 2 * M_PI * drand48();
        double r2  = drand48();
        double r2s = sqrt(r2);

        vector w = nl;
        vector u = ((fabs(w.x) > .1 ? vector(0,1) : vector(1)).cross(w)).normalized();
        vector v = w.cross(u);
        vector d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalized();

        return target->emission + f * (tracer(path(x, d), iter));

    } else if (target->material == SPEC){           // Ideal SPECULAR reflection
        return target->emission + f * (tracer(path(x, ray.end - n * 2 * n.dot(ray.end)), iter));
    } else if (target->material == REFR){

        path reflpath(x, ray.end-n*2*n.dot(ray.end));     // Ideal dielectric REFRACTION

        bool into  = n.dot(nl) > 0;                // path from outside going in?
        double nc  = 1;
        double nt  = 1.5;
        double nnt = into?nc/nt:nt/nc;
        double ddn = ray.end.dot(nl);
        double cos2t;

        if ((cos2t = 1-nnt * nnt * (1-ddn * ddn)) < 0) {   // Total internal reflection
            return target->emission + f * (tracer(reflpath,iter));
        }

        vector tdir  = (ray.end * nnt - n * ((into ? 1: -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
        double a  = nt-nc;
        double b  = nt + nc; 
        double R0 = a * a/(b * b); 
        double c  =  1 -(into ? -ddn : tdir.dot(n));
        double Re = R0 + (1-R0) * c * c * c * c * c;
        double Tr = 1-Re;
        double P = .25 + .5 * Re;
        double RP = Re/P;
        double TP = Tr/(1-P);

        return target->emission + f * (iter>2 ? (drand48()<P ?   // Russian roulette
            tracer(reflpath, iter) * RP : tracer(path(x, tdir), iter) * TP):
            tracer(reflpath, iter) * Re + tracer(path(x, tdir), iter) * Tr);

    }


    // Mirror reflection not working
    //return target->emission + f * (tracer(path(x, ((ray.begin *-1) + (x * 2) +target->normal * (ray.begin - x).dot(target->normal) * 2.0 ).normalized() ), iter));

    //double r1  = 2 * M_PI * drand48();
    //double r2  = drand48();
    //double r2s = sqrt(r2);

    //vector w = nl;
    //vector u = ((fabs(w.x) > .1 ? vector(0,1) : vector(1)).cross(w)).norm();
    //vector v = w.cross(u);
    //vector d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

    //return target->emission + f * (tracer(path(x, d), iter));

    //return target->emission + f * (tracer(path(x, ray.end - target->normal * 2 * (target->normal).dot(ray.end)), iter));
    //return target->emission + f * (tracer(path(x, ray.end - target->normal * 2 * (target->normal).dot(ray.end) + vector(drand48(),drand48(),drand48())), iter));
}

int main(){
    int screenx = 200;
    int screeny = 150;

    int steps = 10;

    path cam(vector(50,52,295.6),
             vector(0,-0.042612,-1).normalized()); // cam pos, dir

    vector dx = vector(screenx * .5135 / screeny); // horizontal increment
    vector dy = dx.cross(cam.end).normalized() * .5135;         // Vertical increment

    color r;

    color *image = new color[screenx * screeny];
    int i;

    int samps = 10;

    for (int y = 0; y < screeny; y++){
        int x;
        if ( y %5 == 0){std::cout<<'\r'<<y;fflush(stdout);}
        for (x = 0; x < screenx; x++){
            r = color();
            for (int sy = 0, i = (screeny - y - 1) * screenx + x; sy < 2; sy++ ){    // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++){        // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++ ){
                        double r1 = 2 * drand48(), ddx = r1<1 ? sqrt(r1) - 1: 1 - sqrt(2 - r1);
                        double r2 = 2 * drand48(), ddy = r2<1 ? sqrt(r2) - 1: 1 - sqrt(2 - r2);
                        vector d  = cam.end + dx * ( ( (sx + .5  +  ddx)/2  +  x)/screenx  -  .5)  +
                                              dy * ( ( (sy + .5  +  ddy)/2  +  y)/screeny  -  .5)  ;
                        r  =  r + tracer(path(cam.begin + d * 140, d.normalized()), 0) * (1. / samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    image[i] = image[i] + color(clamp(r.r),clamp(r.g),clamp(r.b)) * .25;
                    //image[i].print();
                }
            }
            //vector ray = dx * ((double)x/screenx + 0.5) +
                         //dy * ((double)y/screeny + 0.5) + cam.end;

            //int id = -1;
            //double dist = 1<<20 ;
            //bool ya = intersects(path(cam.begin, ray), &id, &dist);
            //std::cout << " id = " << id << " dist = " << dist << '\n';

            //continue;

            //r = r + tracer(path(cam.begin + ray*140, ray.normalized()), 5);
            //i = (screeny - y - 1) * screeny + x;
            //r.print();
            //image[i] = image[i] + r.norm();
            //image[x + y*screenx] = image[x + y*screenx] + r.norm();
            //image[x + y*screenx].print();
        }

        //printf("%f %f %f \n", ((image[i].r)),
                              //((image[i].g)),
                              //((image[i].b)));

    }

    std::cout << "\rDone\n";

    FILE *f = fopen("image.ppm", "w");         // PPM cancer
    fprintf(f, "P3\n%d %d\n%d\n", screenx, screeny, 255);

    for (int i = 0; i < screenx * screeny; i++) {
        fprintf(f,"%d %d %d ", (toInt(image[i].r)),
                               (toInt(image[i].g)),
                               (toInt(image[i].b)));
        //fprintf(f,"%d %d %d ", ((image[i].r)),
                               //((image[i].g)),
                               //((image[i].b)));
    }

    return 0;
}
