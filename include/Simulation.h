
#ifndef SIMULATION_H
#define SIMULATION_H

#include "Globals.h"
#include <SFML/Graphics.hpp>
#include <stdio.h>

class Simulation{

public:

    float u[nY][nX];
    float v[nY][nX];
    float density[nY][nX];
    int s[nY][nX];

    sf::RenderWindow& window;

    Simulation(sf::RenderWindow& window);

    void update();

    void setInitialCondition();

    // this is called every frame to renew the boundary conditions
    void setBoundaryConditions();

    void solveIncompressibility(float dt);

    void advectVelocity(float dt);

    void advectSmoke(float dt);

    void extrapolate();

    float sampleField(float x, float y, Field field);

    // drawing stuff
    void drawDensity();
};

#endif