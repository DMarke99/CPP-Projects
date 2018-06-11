//
//  main.cpp
//  ParticleBoxSimulation
//
//  Created by Diamor Marke on 05/06/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#include "linalg.hpp"
#include "dynamics.hpp"
#include "visualisation.hpp"
const int MAX_FPS = 30;

int main(int argc, const char * argv[]) {

    ParticleBox Simulation(0.8);
    Simulation.add_particles(1000, pow(10,7), 1);
    
    
    //interesting setting
    //n=1000, m=5*10^6, v=1; "galaxy" formation observed in long term behaviour
    
    ParticleScreen Screen("Particle Simulation");
    Screen.update(&Simulation);
    
    //event loop
    while(true){
        long time = clock();
        Screen.update(&Simulation);
        if(Screen.processEvents()){break;}
        while ((double) CLOCKS_PER_SEC/(clock()-time) > MAX_FPS){}
    }
    
    Screen.close();
    
    return 0;
}
