//
//  dynamics.cpp
//  ParticleBoxSimulation
//
//  Created by Diamor Marke on 06/06/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#include "dynamics.hpp"

//particle variables
Particle::Particle(Vector x, Vector v){
    index = ++curr_index;
    this->x = x;
    this->v = v;
}

int Particle::curr_index = -1;

Vector Particle::p(){return m*v;}
Vector Particle::L(Vector origin){return m*(x - origin)^v;}
double Particle::s(){return sqrt(v*v);}
double Particle::T(){return 0.5*m*v*v;}

//defines the triangle wave function for a given period
double triangle(double x, double period){
    if (x >=0){return period/4 - abs(fmod(x+ period/4, period)- period/2);}
    else { return - triangle(-x, period);}
}

int flips(double x, double width){
    return floor((x+width/2)/width);
}

const Vector unit_vec(const Vector& r1, const Vector& r2) {
    Vector r = r1 - r2;
    return (r == Vector({0,0,0})) ? Vector({0,0,0}) : r/sqrt(r*r);
}

 const Vector newton_grav(const double& M, const Vector& r1, const Vector& r2){
    double inv_sq = (r1-r2)*(r1-r2);
    double dist = sqrt(inv_sq);
    if (inv_sq > particle_radius*particle_radius){return (G* (M)/(inv_sq+0.0001)) * unit_vec(r1, r2);}
    else{return (G* (M)/(particle_radius*particle_radius+0.0001)) * unit_vec(r1, r2) * (dist/particle_radius);}
};

ParticleBox::ParticleBox(double width){
    this->width = width;
}

double ParticleBox::getWidth() const {
    return width;
}

Vector ParticleBox::getCentroid() const {
    return centroid;
}

void ParticleBox::add_particles(const int& n, const double& mass, const double& v, Uint8 R, Uint8 G, Uint8 B){
    random_device rand_dev;
    mt19937 generator(rand_dev());
    std::uniform_real_distribution<double>  distr(-width, width);
    
    for (int i = 0; i< n; ++i){
        Particle p;
        p.v = v*Vector::rand();
        
        p.R = R;
        p.G = G;
        p.B = B;
        
        p.m = mass;
        
        p.x(1) = distr(generator);
        p.x(2) = distr(generator);
        p.x(3) = distr(generator);
        
        particles.push_back(p);
    }
}

//used to split calculations into different threads
void update_particle(vector<Particle>& particles,const int start, const int end, const double dt, const double threshold, const Octree p_tree){
    Vector k1, k2, k3;
    for (int i = start; i < end; ++i){
    k1 = dt * p_tree.grav(particles[i].index, particles[i].x, threshold);
    k2 = dt * p_tree.grav(particles[i].index, particles[i].x + 0.5 * dt * particles[i].v
                                 + 0.125 * dt * k1, threshold);
    k3 = dt * p_tree.grav(particles[i].index, particles[i].x + dt * particles[i].v
                                 + 0.5 * dt * k2, threshold);
    
    particles[i].temp_x = particles[i].x + dt * (particles[i].v + (k1 + 2*k2)/6);
    particles[i].v = particles[i].v + k1/6 + 2*k2/3 + k3/6;
    }
};

//uses runge kutta approximations to update particles
void ParticleBox::update(){
    Vector k1, k2, k3;
    vector<thread> threads;
    
    Octree particle_tree = Octree(Vector({0,0,0}), width, particles);
    for (int i = 0; i < MAX_THREADS; ++i){
        int start = (particles.size()*i)/MAX_THREADS;
        int end = (particles.size()*(i+1))/MAX_THREADS;
        threads.push_back(thread(update_particle, ref(particles), start, end, dt, threshold, particle_tree));
    }
    
    for (int i = 0; i < MAX_THREADS; ++i){
        threads[i].join();
    }
    
    for (int i = 0; i < particles.size(); ++i){
        particles[i].x = particles[i].temp_x;
    }
    
    int flip = 0;
    
    for (int i = 0; i < particles.size(); ++i){
        
        //effectively causes a completely elastic bounce off the wall
        for (int j = 1; j<= 3; ++j){
            flip = flips(particles[i].x(j), 2*width);
            particles[i].x(j) = triangle(particles[i].x(j), 4*width);
            particles[i].v(j) = particles[i].v(j)*pow(-1,flip);
            
            
        }
    }
    centroid = particle_tree.getCentroid();
}

vector<Particle> ParticleBox::getParticles(){
    return particles;
};

Octree::Octree(Particle first) : first(first){
    this->mass = first.m;
    this->centroid = first.x;
}

Octree::Octree(Vector center, double width, Particle first){
    this->center = center;
    this->width = width;
    this->first = first;
    this->mass = first.m;
    this->centroid = first.x;
    this->n = 1;
}

Octree::Octree(Vector center, double width, vector<Particle> particles){
    this->center = center;
    this->width = width;
    this->n = 1;
    this->first = (particles[0]);
    this->mass = particles[0].m;
    this->centroid = particles[0].x;
    for (int i = 1; i < particles.size(); ++i){
        this->insert(&particles[i]);
    }
}

int Octree::count() const {
    return n;
}

Vector Octree::getCentroid(){
    return centroid;
}



void Octree::insert(Particle* p){
    ++n;
    
    if(n == 1){centroid = p->x;}
    else{centroid = (mass*centroid+p->x*p->m)/(mass + p->m);}
    
    mass += p->m;
    
    //inserts element if there are more than two elements
    if (n > 2){
        int x = p->x(1) > center(1);
        int y = p->x(2) > center(2);
        int z = p->x(3) > center(3);
        
        if (initialized[x][y][z]) {
            subspaces[x][y][z].insert(p);
        }
        else{
            Vector augment = Vector({static_cast<double>(x-0.5),
                static_cast<double>(y-0.5),
                static_cast<double>(z-0.5)});
            
            subspaces[x][y][z] = Octree(center+width*augment, width/2, *p);
            initialized[x][y][z] = true;
        }
    }
    
    //holds onto particle if it's the ony l particle in the tree
    else if (n == 1){
        first = *p;
    }
    
    //initialises the subspaces within the octree when there are at least two elements
    else if (n == 2){
        Vector new_centre = center - (width/2)*Vector({1,1,1});
        subspaces.resize(2);
        initialized.resize(2);
        for (int i = 0; i < 2; ++i){
            subspaces[i].resize(2);
            initialized[i].resize(2);
            for (int j = 0; j < 2; ++j){
                subspaces[i][j].resize(2);
                initialized[i][j].resize(2);
                for (int k = 0; k < 2; ++k){
                    Vector augment = Vector({static_cast<double>(i),
                                             static_cast<double>(j),
                                             static_cast<double>(k)});
                    initialized[i][j][k] = false;
                }
            }
        }
        
        //inserts first particle
        int x = first.x(1) > center(1);
        int y = first.x(2) > center(2);
        int z = first.x(3) > center(3);

        if (initialized[x][y][z]) {
            subspaces[x][y][z].insert(&first);
        }
        else{
            Vector augment = Vector({static_cast<double>(x-0.5),
                                     static_cast<double>(y-0.5),
                                     static_cast<double>(z-0.5)});
            
            subspaces[x][y][z] = Octree(center+width*augment, width/2, first);
            initialized[x][y][z] = true;
        }
        
        //inserts second particle, handling particles being too close to each other
        if (width < epsilon){
            random_device rand_dev;
            mt19937 generator(rand_dev());
            std::uniform_int_distribution<double>  distr(0, 1);
            x = distr(generator);
            y = distr(generator);
            z = distr(generator);
        }
        else {
            if (first.x == p->x){
                p->x = p->x + min(epsilon,width/4) * Vector::rand();
            }
            x = (p->x(1) > center(1));
            y = (p->x(2) > center(2));
            z = (p->x(3) > center(3));
        }
        
        if (initialized[x][y][z]) {
            subspaces[x][y][z].insert(p);
        }
        else{
            Vector augment = Vector({static_cast<double>(x-0.5),
                                     static_cast<double>(y-0.5),
                                     static_cast<double>(z-0.5)});
            
            subspaces[x][y][z] = Octree(center+width*augment, width/2, *p);
            initialized[x][y][z] = true;
        }

    }
}

Vector Octree::grav(const int& id, const Vector& x, const double& threshold) const {
    Vector r = x - center;
    double dist = sqrt(r*r);
    if (count() == 0){
        return Vector(3);
    }
    else if (count() == 1){
        Vector dist = first.x - x;
        if (sqrt(dist*dist)<particle_radius or id ==first.index){
            return Vector(3);
        }
        else{
            return newton_grav(mass, centroid, x);
        }
    }
    else if(dist/width > threshold){
        return newton_grav(mass, centroid, x);
    }
    else{
        Vector acc = Vector(3);
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 2; ++j){
                for (int k = 0; k < 2; ++k){
                    if (initialized[i][j][k] and subspaces[i][j][k].count() > 0){
                        acc = acc + subspaces[i][j][k].grav(id, x, threshold);
                    }
                }
            }
        }
        return acc;
    }
}
