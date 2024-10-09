
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
    for(int y = 1; y < nY - 1; y++){
        u[y][1] = inVelocity;
    }

    for(int y = 1; y < nY - 1; y++){
        for(int x = 1; x < nX - 1; x++){

            float xPos = (x + 0.5) * h;
            float yPos = (y + 0.5)* h;

            s[y][x] = 1.0;

            for(Circle& circle : circles){

                float dx = circle.center.x - xPos;
                float dy = circle.center.y - yPos;

                if(dx * dx + dy * dy <= circle.radius * circle.radius){
                    s[y][x] = 0;
                    u[y][x] = 0;
                    v[y][x] = 0;
                    density[y][x] = 0;
                }
            }
        }
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

                float div = (u[y][x + 1] - u[y][x]) + (v[y + 1][x] - v[y][x]);

                if(std::abs(div) > 1000.f){
                    printf("abnormal divergence at x: %d, y: %d\n", x, y);
                    printf("u[y][x]: %.3f, u[y][x + 1]: %.3f, v[y][x]: %.3f, v[y + 1][x]: %.3f\n", u[y][x], u[y][x + 1], v[y][x], v[y + 1][x]);
                }

                // every value in staggered grid that can be corrected is corrected by the same amount
                float p = -div / sTotal;
                p *= overrelaxationConstant;

                // contributes negatively so needs to be positive
                u[y][x] -= p * sx0;

                u[y][x + 1] += p * sx1;

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

            newU[y][x] = u[y][x];
            newV[y][x] = v[y][x];
            
            // advect u
            if(s[y][x] != 0 && s[y][x - 1] != 0 && y < nY - 1){
                
                float posX = x * h;
                float posY = (y + 0.5) * h;

                float avgV = (v[y][x] + v[y + 1][x] + v[y][x - 1] + v[y + 1][x - 1]) * 0.25f;

                float prevPosX = posX - u[y][x] * dt;
                float prevPosY = posY - avgV  * dt;

                float interpolated = sampleField(prevPosX, prevPosY, U_FIELD, false);

                newU[y][x] = interpolated;
            }

            if(s[y][x] != 0 && s[y - 1][x] != 0 && x < nX - 1){
                float posX = (x + 0.5) * h;
                float posY = y * h;

                float avgU = (u[y][x] + u[y][x + 1] + u[y - 1][x] + u[y - 1][x + 1]) * 0.25f;

                float prevPosX = posX - avgU * dt;
                float prevPosY = posY - v[y][x] * dt;

                float interpolated = sampleField(prevPosX, prevPosY, V_FIELD, false);

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

void Simulation::moveObstacles(sf::Vector2f& mousePos, sf::Vector2f& dPos){

    sf::Vector2f mousePosInSimCoords = mousePos / cScale;
    sf::Vector2f dPosInSimCoords = dPos / cScale;

    for(Circle& circle : circles){

        sf::Vector2f distVec = circle.center - mousePosInSimCoords;

        if(distVec.x * distVec.x + distVec.y * distVec.y <= circle.radius * circle.radius){

            circle.center += dPosInSimCoords;
        }
    }
}


void Simulation::advectSmoke(float dt){


    float newDensity[nY][nX];

    // advect density by going back
    for(int y = 1; y < nY - 1; y++){
        for(int x = 1; x < nX - 1; x++){

            newDensity[y][x] = density[y][x];
            
            // advect u
            if(s[y][x] != 0){
                
                float posX = (x + 0.5) * h;
                float posY = (y + 0.5) * h;

                float uAtMid = (u[y][x] + u[y][x + 1]) * 0.5;
                float vAtMid = (v[y][x] + v[y + 1][x]) * 0.5;

                float prevPosX = posX - uAtMid * dt;
                float prevPosY = posY - vAtMid  * dt;

                //printf("prevPosX: %.3f, prevPosY: %.3f\n", posX, posY);

                float interpolated = sampleField(prevPosX, prevPosY, D_FIELD, false);

                //printf("interpolated: %.3f\n", interpolated);

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

float Simulation::sampleField(float x, float y, Field field, bool show){

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

    float interpolated = (sx * sy * topLeft) + (tx * sy * topRight) + (sx * ty * bottomLeft) + (tx * ty * bottomRight); 

    if (show) {
        sf::Vector2f topLeftScreenPos = simPosToScreenPos(getSimPos(x0, y0, field));
        sf::Vector2f topRightScreenPos = simPosToScreenPos(getSimPos(x1, y0, field));
        sf::Vector2f bottomLeftScreenPos = simPosToScreenPos(getSimPos(x0, y1, field));
        sf::Vector2f bottomRightScreenPos = simPosToScreenPos(getSimPos(x1, y1, field));

        sf::Color c = fieldToColor(field);

        sf::CircleShape circle(8);
        circle.setFillColor(sf::Color::Transparent);
        circle.setOutlineColor(c);
        circle.setOutlineThickness(2);

        circle.setPosition(topLeftScreenPos);
        window.draw(circle);
        circle.setPosition(topRightScreenPos);
        window.draw(circle);
        circle.setPosition(bottomLeftScreenPos);
        window.draw(circle);
        circle.setPosition(bottomRightScreenPos);
        window.draw(circle);

        sf::Vector2f screenPos = simPosToScreenPos(sf::Vector2f(x, y));
        drawText(window, formatFloat(interpolated), screenPos.x, screenPos.y, c);
    }

    return interpolated;
}

void Simulation::update(){
    setBoundaryConditions();

    float dt = 1.0 / 120;

    float absDivergenceBefore = totalAbsoluteDivergence();

    solveIncompressibility(dt);

    float absDivergenceAfter = totalAbsoluteDivergence();

    //printf("Abs divergence before: %.3f, Abs divergence after: %.3f\n", absDivergenceBefore, absDivergenceAfter);

    advectVelocity(dt);

    advectSmoke(dt);
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

void Simulation::drawObstacles(){

    sf::CircleShape shape(10);
    shape.setFillColor(sf::Color(100, 50, 150, 255));

    for(Circle& circle : circles){
        
        shape.setPosition(circle.center * cScale);
        shape.setOrigin(sf::Vector2f(circle.radius * cScale, circle.radius * cScale));
        shape.setRadius(circle.radius * cScale);
        window.draw(shape);
    }
}

void Simulation::drawBoundary(){

    sf::RectangleShape shape(sf::Vector2f(cellSizeInPixels, cellSizeInPixels));

    // skip the first and last
    for(int y = 0; y < nY; y++){
        for(int x = 0; x < nX; x++){

            if(s[y][x] == 0){
            
            shape.setPosition(sf::Vector2f(x * cellSizeInPixels, y * cellSizeInPixels));
            shape.setFillColor(sf::Color(100, 50, 150, 255));
            window.draw(shape);
            }
        }
    }
}

void Simulation::drawDivergence(){

    sf::RectangleShape shape(sf::Vector2f(cellSizeInPixels, cellSizeInPixels));

    // skip the first and last
    for(int y = 1; y < nY - 1; y++){
        for(int x = 1; x < nX - 1; x++){

            float div = (u[y][x + 1] - u[y][x]) + (v[y + 1][x] - v[y][x]);

            int cVal = int((clamp(div, -1.0f, 1.0f) + 1.0) * 0.5f * 255);
            
            shape.setPosition(sf::Vector2f(x * cellSizeInPixels, y * cellSizeInPixels));
            shape.setFillColor(sf::Color(cVal, cVal, cVal, 255));
            window.draw(shape);
        }
    }
}

void Simulation::drawGrid(sf::RenderWindow& window) {
    // Horizontal grid lines
    for (int y = 0; y <= nY; ++y) {
        sf::Vertex line[] = {
            sf::Vertex(sf::Vector2f(0, y * cellSizeInPixels)),
            sf::Vertex(sf::Vector2f(WINDOW_SIZE, y * cellSizeInPixels))
        };
        line[0].color = sf::Color::Black;
        line[1].color = sf::Color::Black;
        window.draw(line, 2, sf::Lines);
    }

    // Vertical grid lines
    for (int x = 0; x <= nX; ++x) {
        sf::Vertex line[] = {
            sf::Vertex(sf::Vector2f(x * cellSizeInPixels, 0)),
            sf::Vertex(sf::Vector2f(x * cellSizeInPixels, WINDOW_SIZE))
        };
        line[0].color = sf::Color::Black;
        line[1].color = sf::Color::Black;
        window.draw(line, 2, sf::Lines);
    }

    // Draw staggered and colocated
    for (int y = 0; y < nY; ++y) {
        for (int x = 0; x < nX; ++x) {
            sf::CircleShape circle(4);
            circle.setOrigin(4, 4);

            // U field
            sf::Vector2f uPos = getSimPos(x, y, U_FIELD);
            sf::Vector2f uPosOnScreen = simPosToScreenPos(uPos);
            circle.setPosition(uPosOnScreen);
            circle.setFillColor(fieldToColor(U_FIELD));
            window.draw(circle);
            drawText(window, formatFloat(u[y][x]), uPosOnScreen.x, uPosOnScreen.y, sf::Color::Black);

            // V field
            sf::Vector2f vPos = getSimPos(x, y, V_FIELD);
            sf::Vector2f vPosOnScreen = simPosToScreenPos(vPos);
            circle.setPosition(vPosOnScreen);
            circle.setFillColor(fieldToColor(V_FIELD));
            window.draw(circle);
            drawText(window, formatFloat(v[y][x]), vPosOnScreen.x, vPosOnScreen.y, sf::Color::Black);

            // F field (density)
            sf::Vector2f pPos = getSimPos(x, y, D_FIELD);
            sf::Vector2f pPosOnScreen = simPosToScreenPos(pPos);
            circle.setPosition(pPosOnScreen);
            circle.setFillColor(fieldToColor(D_FIELD));
            window.draw(circle);
            drawText(window, formatFloat(density[y][x]), pPosOnScreen.x, pPosOnScreen.y, sf::Color::Black);
        }
    }
}

sf::Color Simulation::fieldToColor(Field field) {
    switch (field) {
        case U_FIELD: return sf::Color::Red;
        case V_FIELD: return sf::Color::Blue;
        case D_FIELD: return sf::Color::Green;
        default: return sf::Color::White;
    }
}

sf::Vector2f Simulation::getSimPos(int x, int y, Field field) {
    float dx = 0, dy = 0;
    switch (field) {
        case U_FIELD: dy = h / 2; break;
        case V_FIELD: dx = h / 2; break;
        case D_FIELD: dx = h / 2; dy = h / 2; break;
    }
    return sf::Vector2f(x * h + dx, y * h + dy);
}

sf::Vector2f Simulation::simPosToScreenPos(sf::Vector2f simPos) {
    return sf::Vector2f(simPos.x * cScale, simPos.y * cScale);
}

void Simulation::drawText(sf::RenderWindow& window, const std::string& text, float x, float y, sf::Color color) {
    static sf::Font font;
    static bool fontLoaded = false;
    if (!fontLoaded) {
        if (!font.loadFromFile("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf")) {
            // Handle font loading error
            return;
        }
        fontLoaded = true;
    }

    sf::Text sfText(text, font, 12);
    sfText.setPosition(x, y);
    sfText.setFillColor(color);
    window.draw(sfText);
}

std::string Simulation::formatFloat(float value) {
    char buffer[10];
    snprintf(buffer, sizeof(buffer), "%.2f", value);
    return std::string(buffer);
}

float Simulation::totalAbsoluteDivergence(){

    float absDivergence = 0;

    // skip the first and last
    for(int y = 1; y < nY - 1; y++){
        for(int x = 1; x < nX - 1; x++){

            float div = (u[y][x + 1] - u[y][x]) + (v[y + 1][x] - v[y][x]);

            absDivergence += std::abs(div);
        }
    }

    return absDivergence;
}