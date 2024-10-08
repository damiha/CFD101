
#include "../include/Globals.h"

float clamp(float val, float min, float max){
    return val > max ? max : (val < min ? min : val);
}

int clamp(int val, int min, int max){
    return val > max ? max : (val < min ? min : val);
}