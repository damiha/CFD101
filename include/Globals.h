
#ifndef GLOBALS_H
#define GLOBALS_H

const int WINDOW_SIZE = 800;
const float SIMULATION_SIZE = 1.0f;

const int numIterationsGaussSeidel = 120;

const float inVelocity = 3.0;

const float overrelaxationConstant = 1.9;

// number of cells
const int nC = 160;

const float h = SIMULATION_SIZE * 1.0 / nC;

const int nX = nC + 2;
const int nY = nC + 2;

const float cScale = WINDOW_SIZE * 1.0 / SIMULATION_SIZE;

const int nCells = nX * nY;

const int cellSizeInPixels = WINDOW_SIZE / nC;

enum Field : unsigned int{
    U_FIELD = 0,
    V_FIELD = 1,
    D_FIELD = 2
};

float clamp(float val, float min, float max);

int clamp(int val, int min, int max);

#endif