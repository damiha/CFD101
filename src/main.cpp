#include <SFML/Graphics.hpp>
#include "../include/Globals.h"
#include "../include/Simulation.h"

int main()
{
    sf::RenderWindow window(sf::VideoMode(WINDOW_SIZE, WINDOW_SIZE), "CFD 101");

    Simulation simulation(window);

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        simulation.update();

        window.clear(sf::Color::White);

        simulation.drawDensity();
        
        window.display();
    }

    return 0;
}