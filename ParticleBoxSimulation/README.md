# Particle Gravity Simulation

This is a Particle Gravity Simulation implemented in C++ over 4 days during Trinity Term.

## Prerequisites

To use the project on your computer you will need:

```
C++ Compiler that supports C++11/C++14
SDL 2.0
```
## Setup

1) Ensure you have a working C++11+ compiler.
2) Download the [Particle Box Simulation Folder](https://github.com/DMarke99/CPP-Projects/tree/master/ParticleBoxSimulation).
3) Execute by compiling [main.cpp](https://github.com/DMarke99/CPP-Projects/blob/master/ParticleBoxSimulation/main.cpp).

## Usage

- Click and hold screen to rotate viewing angle.
- Change settings before the main loop to alter mass per particle and number of particles.

## Implementation

### Acceleration and Velocity Calculation

For calculating acceleration and velocity, I used the Runge-Kutta 4th order method, as I found that this bounded total energy better than linear approximations due to being correct up to order 4 compared to order 2. Directly calculating each force interaction scales like O(n<sup>2</sup>), so I instead used the Barnst-Hut Method, by using an Octree to efficiently approximate forces for clusters of distant bodies, which has time complexity O(n log n).

### Multithreading

I multithreaded the process of velocity and acceleration update, since it requires a large number of calculations being computed which do not depend on each other at the time of calculation. I found a significant increase in performace going from one to 4 threads, however more threads didn't yield a greater increase in performance. This was on a Macbook Pro however, so if this were ran on a computer with more cores and more processing power the performance gains would have been much greater.

### Rendering

I used the SDL framework to render the particle container. I handled rendering in two separate ways; 2D and 3D objects. 2D objects are invariant to the viewing angle, and 3D objects aren't. The coordinates in the chosen viewing angle were found by pre-multiplying the coordinate vector in 3D by an orthogonal matrix then projecting onto R2.

## Contributers

This project was created entirely by myself over 4 days during Trinity Term. The most interesting parts of this project were multithreading the application to improve performance, and learning about efficient computational methods for approximation and calculation of acceleration and velocity.
