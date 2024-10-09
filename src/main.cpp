#include <SFML/Graphics.hpp>
#include "../include/Globals.h"
#include "../include/Simulation.h"

int main()
{
    sf::RenderWindow window(sf::VideoMode(WINDOW_SIZE, WINDOW_SIZE), "CFD 101");

    Simulation simulation(window);

    simulation.circles.push_back(Circle(sf::Vector2f(SIMULATION_SIZE / 2, SIMULATION_SIZE / 2), SIMULATION_SIZE / 8));

    simulation.setBoundaryConditions();

    sf::Vector2f mousePos;
    sf::Vector2f lastMousePos;
    bool mousePressed = false;
    bool simulationStarted = false; 

    while (window.isOpen())
    {

        sf::Vector2i mousePosInteger = sf::Mouse::getPosition(window);
        mousePos = window.mapPixelToCoords(mousePosInteger);

        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();

            if(event.type == sf::Event::MouseButtonPressed){
                lastMousePos = mousePos;
                mousePressed = true;
            }

            if(event.type ==sf::Event::MouseButtonReleased){
                mousePressed = false;
            }

            if(event.type == sf::Event::KeyPressed){
                simulationStarted = true;
            }
        }

        window.clear(sf::Color::White);

        if(mousePressed){
            sf::Vector2f dPos = mousePos - lastMousePos;

            simulation.moveObstacles(mousePos, dPos);
        }

        if(simulationStarted){
            simulation.update();
        }

        simulation.drawDensity();

        simulation.drawBoundary();
        
        window.display();

        lastMousePos = mousePos;
    }

    return 0;
}