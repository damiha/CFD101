
#ifndef SIMULATION_H
#define SIMULATION_H

#include "Globals.h"
#include <SFML/Graphics.hpp>
#include <stdio.h>
#include "Circle.h"

class Simulation{

public:

    float u[nY][nX];
    float v[nY][nX];
    float density[nY][nX];
    int s[nY][nX];

    std::vector<Circle> circles;

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

    float sampleField(float x, float y, Field field, bool show);

    float totalAbsoluteDivergence();

    // drawing stuff
    void drawDensity();

    void drawBoundary();

    void drawDivergence();

    void drawObstacles();

    void moveObstacles(sf::Vector2f& mousePos, sf::Vector2f& dPos);

    sf::Color fieldToColor(Field field);
    sf::Vector2f getSimPos(int x, int y, Field field);
    sf::Vector2f simPosToScreenPos(sf::Vector2f simPos);
    void drawText(sf::RenderWindow& window, const std::string& text, float x, float y, sf::Color color);
    void drawGrid(sf::RenderWindow& window);

    std::string formatFloat(float value);
};

#endif