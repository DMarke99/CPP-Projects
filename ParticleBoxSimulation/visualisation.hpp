//
//  visualisation.hpp
//  ParticleBoxSimulation
//
//  Created by Diamor Marke on 06/06/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#ifndef visualisation_hpp
#define visualisation_hpp

#include <SDL2/SDL.h>
#include "dynamics.hpp"
using namespace std;

//class used to draw in SDL
class Screen{
public:
    const static int SCREEN_WIDTH = 600;
    const static int SCREEN_HEIGHT = 600;
protected:
    SDL_Window *m_window;
    SDL_Renderer *m_renderer;
    SDL_Texture *m_texture;
    Uint32 *m_buffer;
    
public:
    Screen(string title);
    void update();
    bool processEvents();
    void close();
    
    void fillScreen(Uint8 shade);
    void drawCircle(int x, int y, int radius,Uint8 R, Uint8 G, Uint8 B, bool filled = true);
    void drawDot(int x, int y,Uint8 R, Uint8 G, Uint8 B);
    void drawLine2D(int x0, int y0, int x1, int y1,Uint8 R, Uint8 G, Uint8 B);
    
    int square_distance(int x1,int x2,int y1,int y2);
    double error_from_line(int x0,int y0,double a,double b,double c);
    int signof(double i);
    
    void setPixel(int x, int y, Uint8 R, Uint8 G, Uint8 B);
    Uint32 getColor(Uint8 R, Uint8 G, Uint8 B);
};

//class used to render ParticleBox
class ParticleScreen: public Screen{
private:
    int timejumps = 1;
    int pixels_per_metre = 200;
    Matrix ModelView = Matrix::identity(3);      //the matrix that represents the current orientation of the viewing angle
    
    int frames = 0;                      //frame counter
    double fps_update = 0.5;             //the minimum time in which the fps counter should be updated
    int ticks;                           //tick counter
    bool plot_centroid = false;
public:
    ParticleScreen(string title);
    void update(ParticleBox* P);
    bool operator ()(const Particle first, const Particle second);
    void drawAxis();
    void drawLine3D(Vector vertex1, Vector vertex2, Uint8 R, Uint8 G, Uint8 B);
    void drawPolygon(vector<Vector> vertices, Uint8 R, Uint8 G, Uint8 B);
    void drawContainer(double size);
};
#endif /* visualisation_hpp */
