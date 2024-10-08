
#include "../include/Simulation.h"

Simulation::Simulation(sf::RenderWindow& window): window(window){
    setInitialCondition();
}

void Simulation::setInitialCondition(){
    
    // set everything to zero at first
    for(int y = 0; y < nY; y++){
        for(int x = 0; x < nX; x++){

            s[y][x] = 1;

            // the top, bottom and left should be solid
            if(x == 0 || y == 0 || y == nY - 1){
                s[y][x] = 0;
            }

            u[y][x] = 0;
            v[y][x] = 0;
            density[y][x] = 0;
        }
    }

    printf("INFO: Fields initialized\n");
}

// this is called every frame to renew the boundary conditions
void Simulation::setBoundaryConditions(){

    // YOUR BOUNDARY CONDITIONS COULD BE HERE
    float relativeStreamWidth = 0.2;

    int streamWidthInCells = int(nC * relativeStreamWidth);

    int startStreamAt = int(nC * 0.5) - streamWidthInCells / 2;
    int endStreamAt = startStreamAt + streamWidthInCells;

    for(int y = startStreamAt; y <= endStreamAt; y++){
        density[y][1] = 1.0;
    }

    // wind blows from left to right
    for(int y = 0; y < nY; y++){
        u[y][1] = 1.0;
    }
}

void Simulation::solveIncompressibility(float dt){

    for(int i = 0; i < numIterationsGaussSeidel; i++){
        
        for(int y = 1; y < nY - 1; y++){
            for(int x = 1; x < nX - 1; x++){
                
                if(s[y][x] == 0){
                    continue;
                }

                int sx0 = s[y][x - 1];
                int sx1 = s[y][x + 1];
                int sy0 = s[y - 1][x];
                int sy1 = s[y + 1][x];

                int sTotal = sx0 + sx1 + sy0 + sy1;

                // surrounded by walls?
                if(sTotal == 0){
                    continue;
                }

                float div = (u[y][x + 1] * sx1 - u[y][x] * sx0) + (v[y + 1][x] * sy1 - v[y][x] * sy0);

                // every value in staggered grid that can be corrected is corrected by the same amount
                float p = -div / sTotal;

                // contributes negatively so needs to be positive
                u[y][x] -= p * sx0;

                u[y + 1][x] += p * sx1;

                v[y][x] -= p * sy0;

                v[y + 1][x] += p * sy1;
            }
        }
    }
}

void Simulation::advectVelocity(float dt){

    float newU[nY][nX];
    float newV[nY][nX];

    // advect velocity by going back
    for(int y = 1; y < nY; y++){
        for(int x = 1; x < nX; x++){
            
            // advect u
            if(s[y][x] != 0 && s[y][x - 1] != 0 && y < nY - 1){
                
                float posX = x * h;
                float posY = (y + 0.5) * h;

                float avgV = (v[y][x] + v[y][x - 1] + v[y + 1][x - 1] + v[y + 1][x]) * 0.25f;

                float prevPosX = posX - u[y][x] * dt;
                float prevPosY = posY - avgV  * dt;

                float interpolated = sampleField(prevPosX, prevPosY, U_FIELD);

                newU[y][x] = interpolated;
            }

            if(s[y][x] != 0 && s[y - 1][x] != 0 && x < nX - 1){
                float posX = (x + 0.5) * h;
                float posY = y * h;

                float avgU = (u[y-1][x] + u[y - 1][x] + u[y][x] + u[y][x + 1]) * 0.25f;

                float prevPosX = posX - avgU * dt;
                float prevPosY = posY - v[y][x];

                float interpolated = sampleField(prevPosX, prevPosY, V_FIELD);

                newV[y][x] = interpolated;
            }
        }
    }

    // replace old by new (only where new values are set)
    for(int y = 1; y < nY; y++){
        for(int x = 1; x < nX; x++){
            u[y][x] = newU[y][x];
            v[y][x] = newV[y][x];
        }
    }
}

void Simulation::advectSmoke(float dt){


    float newDensity[nY][nX];

    // advect density by going back
    for(int y = 1; y < nY - 1; y++){
        for(int x = 1; x < nX - 1; x++){
            
            // advect u
            if(s[y][x] != 0){
                
                float posX = (x + 0.5) * h;
                float posY = (y + 0.5) * h;

                float uAtMid = (u[y][x] + u[y][x + 1]) * 0.5;
                float vAtMid = (v[y][x] + v[y + 1][x]) * 0.5;

                float prevPosX = posX - uAtMid * dt;
                float prevPosY = posY - vAtMid  * dt;

                float interpolated = sampleField(prevPosX, prevPosY, D_FIELD);

                newDensity[y][x] = interpolated;
            }
        }
    }

    // replace old by new (only where new values are set)
    for(int y = 1; y < nY - 1; y++){
        for(int x = 1; x < nX - 1; x++){
            density[y][x] = newDensity[y][x];
        }
    }
}

void Simulation::extrapolate(){

}

float Simulation::sampleField(float x, float y, Field field){

    x = clamp(x, h, nX * h);
    y = clamp(y, h, nY * h);

    float dx = 0;
    float dy = 0;

    if(field == U_FIELD){
        dy = h * 0.5;
    }
    else if(field == V_FIELD){
        dx = h * 0.5;
    }
    else{
        dx = h * 0.5;
        dy = h * 0.5;
    }

    int x0 = int((x - dx) / h);
    x0 = clamp(x0, 0, nX - 1);
    int x1 = clamp(x0 + 1, 0, nX - 1);

    int y0 = int((y - dy) / h);
    y0 = clamp(y0, 0, nY - 1);
    int y1 = clamp(y0 + 1, 0, nY - 1);

    // bilinear interpolation
    float tx = ((x - dx) - x0 * h) / h;
    float ty = ((y - dy) - y0 * h) / h;

    float sx = 1.0 - tx;
    float sy = 1.0 - ty;

    // set these subsequently
    float topLeft = -1;
    float topRight = -1;
    float bottomLeft = -1;
    float bottomRight = -1;

    if(field == U_FIELD){
        topLeft = u[y0][x0];
        topRight = u[y0][x1];
        bottomLeft = u[y1][x0];
        bottomRight = u[y1][x1];
    }
    else if(field == V_FIELD){
        topLeft = v[y0][x0];
        topRight = v[y0][x1];
        bottomLeft = v[y1][x0];
        bottomRight = v[y1][x1];
    }
    else if(field == D_FIELD){
        topLeft = density[y0][x0];
        topRight = density[y0][x1];
        bottomLeft = density[y1][x0];
        bottomRight = density[y1][x1];
    }

    return (sx * sy * topLeft) + (tx * sy * topRight) + (sx * ty * bottomLeft) + (tx * ty * bottomRight); 
}

void Simulation::update(){
    setBoundaryConditions();
}

void Simulation::drawDensity(){

    sf::RectangleShape shape(sf::Vector2f(cellSizeInPixels, cellSizeInPixels));

    // skip the first and last
    for(int y = 1; y < nY - 1; y++){
        for(int x = 1; x < nX - 1; x++){

            int cVal = int((1. - density[y][x]) * 255);
            
            shape.setPosition(sf::Vector2f(x * cellSizeInPixels, y * cellSizeInPixels));
            shape.setFillColor(sf::Color(cVal, cVal, cVal, 255));
            window.draw(shape);
        }
    }
}