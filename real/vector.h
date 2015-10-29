#ifndef __VECTOR_PT
#define __VECTOR_PT

class Vector{
    private:

    public:
        double x, y, z;

        Vector(double xi = 0, double yi = 0, double zi = 0){
            x = xi;
            y = yi;
            z = zi;
        }

        Vector(const Vector &u){
            x = u.x;
            y = u.y;
            z = u.z;
        }

        Vector operator + (Vector v);
        Vector operator + (double k);
        Vector operator - (double k);
        Vector operator - (Vector v);
        Vector operator / (double k);
        Vector operator * (double k);
        Vector operator * (Vector v);
        double dot(Vector u);
        double norm();
        Vector copy();
        Vector cross (Vector u);
        Vector normalized();
        void   normalize();
        void   print();
};

#endif /* __VECTOR_PT */
