#include <iostream>
#include <cmath>
#include "vector.h"

Vector Vector::operator + (Vector v){
    return Vector(x + v.x,
                    y + v.y,
                    z + v.z);
}

Vector Vector::operator + (double k){
    return Vector(x + k,
                    y + k,
                    z + k);
}

Vector Vector::operator - (double k){
    return Vector(x - k,
                    y - k,
                    z - k);
}

Vector Vector::operator - (Vector v){
    return Vector(x - v.x,
                    y - v.y,
                    z - v.z);
}

Vector Vector::operator / (double k){
    return Vector(x / k,
                    y / k,
                    z / k);
}

Vector Vector::operator * (double k){
    return Vector(x * k,
                    y * k,
                    z * k);
}

Vector Vector::operator * (Vector v){
    return Vector(x * v.x,
                    y * v.y,
                    z * v.z);
}

double Vector::dot(Vector u){ // Produto escalar
    return x * u.x +
            y * u.y +
            z * u.z;
}

double Vector::norm(){
    return sqrt( x * x +
                    y * y +
                    z * z);
}

Vector Vector::copy(){
    return Vector(x,
                    y,
                    z);
}

void Vector::normalize(){
    *this = *this / norm();
    return;
}

Vector Vector::cross (Vector u){        // Produto vetorial
    return Vector(y * u.z - z * u.y,
                    z * u.x - x * u.z,
                    x * u.y - y * u.x);
}

Vector Vector::normalized(){
    Vector w = *this / norm();
    return w.copy();
}

void Vector::print(){
    std::cout << "( " << x << ", " << y << ", " << z << " )\n";
}
