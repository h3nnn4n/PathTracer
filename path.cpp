#include "path.h"

Vector Path::getPoint(double t){
    return begin + (end - begin) * t;
}

Vector Path::getDir(){
    return end - begin;
}
