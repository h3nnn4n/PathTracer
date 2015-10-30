#include <iostream>
#include <cmath>

#include "color.h"

void Color::print(){
    std::cout << "( " << r << ", " << g << ", " << b << " )\n";
}

Color Color::truncated(){
    return Color((*this).truncate(r),
                 (*this).truncate(g),
                 (*this).truncate(b));
}

double Color::truncate(double w){
    double q;

    if ( w < 0.0 ){
        q = 0.0;
    } else if (w > 1.0) {
        q = 1.0;
    } else {
        q = w;
    }

    return q;
}

Color Color::gammaCorrection(){
    Color p = (*this).truncated();

    return Color((pow( (*this).truncate(r), 1.0 / 2.2) * 255.0 + .5 ),
                 (pow( (*this).truncate(g), 1.0 / 2.2) * 255.0 + .5 ),
                 (pow( (*this).truncate(b), 1.0 / 2.2) * 255.0 + .5 ));
}

Color Color::operator + (Color v){
    return Color(r + v.r,
                 g + v.g,
                 b + v.b);
}

Color Color::operator * (double k){
    return Color(r * k,
                 g * k,
                 b * k);
}

Color Color::operator * (Color v){
    return Color(r * v.r,
                 g * v.g,
                 b * v.b);
}
