#include "color.h"
#include <iostream>

int gammaCorrection(double);

void Color::print(){
    std::cout << "( " << r << ", " << g << ", " << b << " )\n";
}

Color Color::operator + (Color v){
    return Color(r + v.r,
            g + v.g,
            b + v.b);
}

Color Color::norm(){
    return Color(gammaCorrection(r),
                 gammaCorrection(g),
                 gammaCorrection(b));
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
