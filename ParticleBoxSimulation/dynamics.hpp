//
//  dynamics.hpp
//  ParticleBoxSimulation
//
//  Created by Diamor Marke on 06/06/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#ifndef dynamics_hpp
#define dynamics_hpp
#include <SDL2/SDL.h>
#include <thread>
#include "linalg.hpp"

const double G = pow(10,-10);                  //gravitational constant
const double epsilon = pow(10,-30);            //determines how much to move a particle when two particles occupy the same position
const double particle_radius = pow(10,-3);    //the minimum distance at which forces are calculated between particles
const int MAX_THREADS = 4;

//defines an individual particle
class Particle{
public:
    Particle(Vector x = Vector({0,0,0}), Vector v = Vector({0,0,0}));
    static int curr_index;
    int index = 0;                        //unique ID
    double m = 1;                         //mass
    Vector x = Vector({0,0,0});           //position
    Vector temp_x;                        //temp position used in calculation
    Vector v = Vector({0,0,0});           //velocity
    Vector a = Vector({0,0,0});           //acceleration
    Vector p();                           //linear momentum
    Vector L(Vector origin = Vector(3));  //angular momentum
    double s();                           //speed
    double T();                           //kinetic energy
    
    
    Uint8 R = 255;                        //particle color
    Uint8 G = 255;
    Uint8 B = 255;
};

const Vector unit_vec(const Vector& r1, const Vector& r2);                      //unit vector
const Vector newton_grav(const double& M, const Vector& r1, const Vector& r2);  //newton's law of gravitation
double triangle(double x, double period); //defines the triangle wave. Used for totally elastic collisions with wall
int flips(double x, double width); //determines the number of times to flip velocity at each iteration

//defines the container in which particles reside
class ParticleBox{
private:
    vector<Particle> particles;
    double width;                           //distance between centre and wall
    double dt = 0.05;                       //time interval
    double threshold = 1.1;                   //threshold for force approximation
    Vector centroid = Vector({0,0,0});
public:
    ParticleBox(double width);
    void add_particles(const int& n, const double& mass, const double& v,
                       Uint8 R=255, Uint8 G=255, Uint8 B=255);
    vector<Particle> getParticles();
    Vector getCentroid() const;
    double getWidth() const;
    void update();                          //updates positions based on forces calculated by barnes hut algorithm and ODEs solved by the runge-kutta 4th order approximation
};

//defines the tree structure used for force calculations
class Octree{
private:
    int n = 0;
    double mass = 0;                        //total mass
    Vector center = Vector(3);
    Vector centroid;                        //centre of mass
    double width = 0;
    Particle first;
    vector<vector<vector<Octree>>> subspaces;
    vector<vector<vector<bool>>> initialized;
public:
    Octree(){};
    Octree(Particle first);
    Octree(Vector center, double width, Particle first);
    Octree(Vector center, double width, vector<Particle> particles);
    int count() const;
    void insert(Particle* p);
    Vector getCentroid();
    Vector grav(const int& id, const Vector& x, const double& threshold) const;  //calculates the acceration of a particle at x due to the particles in the tree
};

void update_particle(vector<Particle>& p, const int start, const int end, const double dt, const double threshold, const Octree p_tree);
#endif /* dynamics_hpp */
