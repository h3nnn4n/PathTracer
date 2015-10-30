#ifndef __COLOR_PT
#define __COLOR_PT

class Color{
    private:

    public:
        double r, g, b;
        Color(double rr = 0, double gg = 0, double bb = 0){
            r = rr;
            g = gg;
            b = bb;
        }

        Color  operator + (Color v);
        Color  operator * (Color v);
        Color  operator * (double k);
        Color  gammaCorrection();
        double truncate(double);
        Color  truncated();
        void   print();
        Color  norm();
};

#endif /* __COLOR_PT */
