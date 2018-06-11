//
//  visualisation.cpp
//  ParticleBoxSimulation
//
//  Created by Diamor Marke on 06/06/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#include "visualisation.hpp"

Screen::Screen(string title){
    if(SDL_Init(SDL_INIT_VIDEO)<0){
        cout << "graphics failed";
        throw exception();
    }
    
    m_window = SDL_CreateWindow("Particle Experiment", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    
    if (m_window == NULL){
        SDL_Quit();
    }
    
    m_renderer = SDL_CreateRenderer(m_window, -1, SDL_RENDERER_PRESENTVSYNC);
    m_texture = SDL_CreateTexture(m_renderer,SDL_PIXELFORMAT_RGBA8888,SDL_TEXTUREACCESS_STATIC,SCREEN_WIDTH,SCREEN_HEIGHT);
    
    if(m_renderer == NULL){
        close();
    }
    
    m_buffer = new Uint32[SCREEN_WIDTH*SCREEN_HEIGHT];
}

void Screen::update(){
    SDL_UpdateTexture(m_texture, NULL, m_buffer, SCREEN_WIDTH*sizeof(Uint32));
    SDL_RenderClear(m_renderer);
    SDL_RenderCopy(m_renderer,m_texture,NULL,NULL);
    SDL_RenderPresent(m_renderer);
}

bool Screen::processEvents(){
    SDL_Event event;
    
    while(SDL_PollEvent(&event)){
        if(event.type == SDL_QUIT){
            return true;
        }
    }
    
    return false;
}

void Screen::close(){
    SDL_DestroyRenderer(m_renderer);
    SDL_DestroyTexture(m_texture);
    SDL_DestroyWindow(m_window);
    SDL_Quit();
}

void Screen::fillScreen(Uint8 shade){
    memset(m_buffer,shade,SCREEN_WIDTH*SCREEN_HEIGHT * sizeof(Uint32));
}

void Screen::drawCircle(int x, int y, int radius,Uint8 R, Uint8 G, Uint8 B, bool filled){
    if(filled){
        for (int _x = x - radius; _x <= x + radius;_x++){
            if (_x < 0 or _x >= SCREEN_WIDTH){
                continue;
            }
            for (int _y = y - radius; _y <= y + radius; _y++){
                if (_y < 0 or _y >= SCREEN_HEIGHT){
                    continue;
                }
                else if ((_x-x)*(_x-x) + (_y-y)*(_y-y) <= radius*radius){
                    setPixel(_x, _y, R, G, B);
                }
                
            }
        }
    }
}

void Screen::drawDot(int x, int y,Uint8 R, Uint8 G, Uint8 B){
    setPixel(x, y, R, G, B);
}

void Screen::drawLine2D(int x0, int y0, int x1, int y1,Uint8 R, Uint8 G, Uint8 B){
    if (x0 == x1){
        for (int i = min(y0,y1);i <= max(y0,y1);i++){
            setPixel(x0,i,R,G,B);
        }
    }
    else if(y0 == y1){
        for (int i = min(x0,x1);i <= max(x0,x1);i++){
            setPixel(i,y0,R,G,B);
        }
    }
    else{
        int x = min(x0,x1);
        int y = (x == x0) ? y0 : y1;
        
        int x_bound = (x == x0) ? x1 : x0;
        int y_bound = (y == y0) ? y1 : y0;
        
        int a = x0 - x1;
        int b = y0 - y1;
        double m = double(b)/double(a);
        double c = y0 - m * x0;
        int direction = signof(m);
        
        setPixel(x,y,R,G,B);
        
        while ((x <= x_bound)&&(direction*y<=y_bound)){
            if((error_from_line(x+1, y+direction, m, -1, c) <= error_from_line(x, y+direction, m, -1, c))&&(error_from_line(x+1, y+direction, m, -1, c) <= error_from_line(x+1, y, m, -1, c))){
                x+=1;
                y+=direction;
                setPixel(x,y,R,G,B);
            }
            else if(error_from_line(x+1, y, m, -1, c) <= error_from_line(x, y+direction, m, -1, c)){
                x+=1;
                setPixel(x,y,R,G,B);
            }
            else{
                y+=direction;
                setPixel(x,y,R,G,B);
            }
        }
    }
}

int Screen::square_distance(int x1,int x2,int y1,int y2){
    return (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
}

double Screen::error_from_line(int x0,int y0,double a,double b,double c){
    return abs(a * x0 + b * y0 + c);
}

int Screen::signof(double i){
    if (i == 0){return 0;}
    else if( i > 0){return 1;}
    else{return -1;}
}

void Screen::setPixel(int x, int y, Uint8 R, Uint8 G, Uint8 B){
    if (x <SCREEN_WIDTH & y < SCREEN_HEIGHT & x >= 0 & y >= 0){
        m_buffer[(y * SCREEN_WIDTH)+x] = getColor(R,G,B);
    }
}

Uint32 Screen::getColor(Uint8 R, Uint8 G, Uint8 B){
    Uint32 color = 0;
    color += R;
    color <<= 8;
    color += G;
    color <<= 8;
    color += B;
    color <<= 8;
    color += 0xFF;
    return color;
}

ParticleScreen::ParticleScreen(string title):Screen(title){
    ticks = SDL_GetTicks();
}

void ParticleScreen::update(ParticleBox* P){
    int x = 0;
    int y = 0;
    
    //click to rotate viewing angle
    if(SDL_GetMouseState(&x, &y)&SDL_BUTTON(1)){;
    
    double theta1 = (x - SCREEN_WIDTH/2)*0.005;
    double theta2 = (y - SCREEN_HEIGHT/2)*0.005;
    
    ModelView = Matrix::rotate_x(theta2) * Matrix::rotate_y(theta1);
    }
    //updates positions
    for (int i = 0; i < timejumps; i++){
        P->update();
    }
    //clears buffer
    memset(m_buffer,0,SCREEN_WIDTH*SCREEN_HEIGHT * sizeof(Uint32));
    
    //draws particles
    vector<Particle> particles = P->getParticles();
    
    //sorts particles by proximity to front of view
    sort(particles.begin(), particles.end(),*this);
    
    drawAxis();
    drawContainer(P->getWidth());
    
    Vector pos;

    for(Particle p : particles){
        pos = ModelView * p.x;
        
        x = SCREEN_WIDTH/2 + pixels_per_metre*pos(1);
        y = SCREEN_HEIGHT/2 + pixels_per_metre*pos(2);
        this->setPixel(x, y, p.R, p.G, p.B);
    }
    
    if (plot_centroid){
        pos = ModelView * P->getCentroid();
        x = SCREEN_WIDTH/2 + pixels_per_metre*pos(1);
        y = SCREEN_HEIGHT/2 + pixels_per_metre*pos(2);
        this->drawLine2D(x-5, y-5, x+5, y+5, 255, 0, 0);
        this->drawLine2D(x+5, y-5, x-5, y+5, 255, 0, 0);
    }

    
    //updates screen
    this->Screen::update();
    
    //calculates fps by comparing the number of ticks with the number of frames in a given time period
    ++frames;
    if (SDL_GetTicks() - ticks > 1000*fps_update){
        double time = (SDL_GetTicks() - ticks);
        time = time / 1000;
        int fps = frames/time;
        string title = "Particle Experiment; FPS: " + to_string(fps);
        SDL_SetWindowTitle(m_window, title.c_str());
        ticks = SDL_GetTicks();
        frames = 0;
    }
}

//comparison for particles to determine which particle is infront of another
bool ParticleScreen::operator ()(const Particle first, const Particle second){
    return ((ModelView * first.x)(3) > (ModelView * second.x)(3));
}

void ParticleScreen::drawAxis(){
    Vector x = 0.05*Vector({1,0,0});
    Vector y = 0.05*Vector({0,1,0});
    Vector z = 0.05*Vector({0,0,1});
    
    drawLine3D(Vector(3), x, 255, 0, 0);
    drawLine3D(Vector(3), y, 0, 255, 0);
    drawLine3D(Vector(3), z, 0, 0, 255);
    drawLine3D(Vector(3), (-1)*z, 0, 0, 255);
    drawLine3D(Vector(3), (-1)*y, 0, 255, 0);
    drawLine3D(Vector(3), (-1)*x, 255, 0, 0);
}

//draws a line in 3D space that changes with the viewing angle
void ParticleScreen::drawLine3D(Vector vertex1, Vector vertex2, Uint8 R, Uint8 G, Uint8 B){
    Vector pos1 = ModelView * vertex1;
    Vector pos2 = ModelView * vertex2;
    
    if ((pos1(1)!=pos2(1))||(pos1(2)!=pos2(2))){
        drawLine2D(SCREEN_WIDTH/2 + pixels_per_metre*pos1(1),
                   SCREEN_HEIGHT/2 + pixels_per_metre*pos1(2),
                   SCREEN_WIDTH/2 + pixels_per_metre*pos2(1),
                   SCREEN_HEIGHT/2 + pixels_per_metre*pos2(2),
                   R,G,B);
    }
}

//draws a polygon with edges defined by adjacent vertices
void ParticleScreen::drawPolygon(vector<Vector> vertices, Uint8 R, Uint8 G, Uint8 B){
    if (vertices.size() <= 1) return;
    for (int i = 0; i< vertices.size()-1; ++i){
        drawLine3D(vertices[i], vertices[i+1], R, G, B);
    }
    drawLine3D(vertices[vertices.size()-1], vertices[0], R, G, B);
}

//draws cube that contains particles
void ParticleScreen::drawContainer(double size){
    vector<Vector> front;
    vector<Vector> back;
    vector<Vector> other;
    
    front.push_back(size*Vector({1, 1, 1}));
    front.push_back(size*Vector({1, -1, 1}));
    front.push_back(size*Vector({1, -1, -1}));
    front.push_back(size*Vector({1, 1, -1}));
    drawPolygon(front, 255, 255, 255);
    
    back.push_back(size*Vector({-1, 1, 1}));
    back.push_back(size*Vector({-1, -1, 1}));
    back.push_back(size*Vector({-1, -1, -1}));
    back.push_back(size*Vector({-1, 1, -1}));
    drawPolygon(back, 255, 255, 255);
    
    for (int i = 0; i <= 1; ++i){
        for (int j = 0; j <= 1; ++j){
            other.clear();
            other.push_back(size*Vector({1,
                static_cast<double>(2*i-1),
                static_cast<double>(2*j-1)}));
            other.push_back(size*Vector({-1,
                static_cast<double>(2*i-1),
                static_cast<double>(2*j-1)}));
            drawPolygon(other, 255, 255, 255);
        }
    }
    
}
