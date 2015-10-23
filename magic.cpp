#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

int gammaCorrection(double);

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

        void print(){
            std::cout << "( " << x << ", " << y << ", " << z << " )\n";
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
    double size = sizeof(planes) / sizeof(plane);
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

color tracer(path ray, int iter){
    double distance;
    int id = 0;

    if (!intersects(ray, &id, &distance)){
        return color();
    }

    plane *target = &planes[id];

    vector x     = ray.begin + ray.end * distance;
    vector norma = (x - target->point).normalized();
    //vector norma = (target->normal);
    color f      = target->c;
    vector nl;

    if ( norma.dot(ray.end) < 0 ){
        nl = norma;
    } else {
        nl = norma * -1.0;
    }

    double p = f.r > f.g && f.r > f.b ? f.r : f.g > f.b ? f.g : f.b; // max refl

    if (iter-- > 0){
        if ( drand48() < p ){
            f = f * (1 / p);
        } else {
            return target->emission;
        }
    }

    if ( iter < -15 ){
        return color();
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
    return target->emission + f * (tracer(path(x, x - target->normal * 2 * (target->normal).dot(x) + vector(drand48(),drand48(),drand48())), iter));
}

int main(){
    //plane  p(vector(0, 0, 0), vector(1, 0, 0));
    //vector wx(0.2, 0.2, 0.2);
    //path   u(vector(1, 1, 1)+wx, vector(0, 0, 0)+wx);

    //vector x = u.getPoint(p.intersect(u));

    //x.print();

    //vector w = u.end - p.normal * 2* p.normal.dot(u.end);

    //w.print();

    //return 0;

    int screenx = 80;
    int screeny = 60;

    int steps = 10;

    path cam(vector(5, 5, 1.6),
             vector(1, 1, 1).normalized()); // cam pos, dir


    //path cam(vector(15, 15, 15), vector(10, 10, 10).normalized()); // Camera position: pos and direction
    vector dx = vector(screenx * .5135 / screeny); // horizontal increment
    vector dy = dx.cross(cam.end) * .5135;         // Vertical increment

    color r;

    color *image = new color[screenx * screeny];

    for (int y = 0; y < screeny; y++){
        int x;
        for (x = 0; x < screenx; x++){
            r = color();
            vector ray = dx * ((double)x/screenx + 0.5) +
                         dy * ((double)y/screeny + 0.5) + cam.end;

            int id = -1;
            double dist = 1<<20 ;
            bool ya = intersects(path(cam.begin, ray), &id, &dist);
            std::cout << " id = " << id << " dist = " << dist << '\n';

            //continue;

            r = r + tracer(path(cam.begin + ray*140, ray.normalized()), 5);
            //r.norm().print();
            image[x + y*screenx] = image[x + y*screenx] + r.norm();
            //image[x + y*screenx].print();
        }

    }

    std::cout << "Done\n";

    FILE *f = fopen("image.ppm", "w");         // PPM cancer
    fprintf(f, "P3\n%d %d\n%d\n", screenx, screeny, 255);

    for (int i = 0; i < screenx * screeny; i++) {
        fprintf(f,"%d %d %d ", (image[i].r),
                               (image[i].g),
                               (image[i].b));
    }

    return 0;
}
