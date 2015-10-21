#include <iostream>
#include <cmath>

class vector{
    private:
        double x, y, z;

    public:
        vector(double xi, double yi, double zi) : x(xi), y(yi), z(zi) {}
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
            return x * u.x + y + u.y + z * u.z;
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
        color(double rr, double gg, double bb) : r(rr), g(gg), b(bb) {}

};

class path{
    private:

    public:
        vector begin, end;

        path(vector b, vector e) : begin(b), end(e) {}
};

class plane{
    private:
        vector normal;
        vector point;
        color  c;
    public:
        plane(vector p, vector n, color cc) : normal(n), point(p), c(cc) {}

        double intersect(path p){
            return ((point - p.begin).dot(normal)) / p.end.dot(normal);
        }
};

plane planes[] = {
    plane( vector(-1, 0, -1), vector(0, 1, 0), color(0.999, 0.999, 0.909)),
    plane( vector(-1, 0, -1), vector(0, 1, 0), color(0.999, 0.999, 0.909))
};

bool intersects( path p, double t, int *n, double *dist){
    double size = sizeof(planes) / sizeof(plane);
    double distance;
    double dold = 1<<20;

    for(int i = 0 ; i < size ; i++){
        dold = planes[i].intersect(p);
        if (dold && dold < t){
            *n = i;
            t = dold;
        }
    }

    *dist = dold;

    return t < 1<<20;
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

int main(){

    return 0;
}
