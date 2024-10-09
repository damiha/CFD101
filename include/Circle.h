
#ifndef CIRCLE_H
#define CIRLCE_H

#include <SFML/System/Vector2.hpp>

class Circle{

public:
    sf::Vector2f center;
    float radius;

    Circle(sf::Vector2f center, float radius) : center(center), radius(radius) {};
};

#endif