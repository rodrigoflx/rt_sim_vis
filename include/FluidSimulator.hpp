#ifndef FLUID_SIMULATOR_HPP
#define FLUID_SIMULATOR_HPP

#include <vector>

namespace FluidSim {

struct Particle {
    float x, y;
};

class FluidSimulator {
public:
    FluidSimulator(int nx, int ny, float dt, float dx);
    
    // Advances the simulation by one time step.
    void step();
    
    const std::vector<Particle>& getParticles() const;

private:
    // Simulation parameters.
    int nx, ny;
    float dt, dx;
    static constexpr int pressureIterations = 40;
    static constexpr float gravity = -9.8f;

    // Staggered grid for velocities:
    // u: horizontal velocity at vertical cell faces (size: (nx+1)*ny)
    // v: vertical velocity at horizontal cell faces (size: nx*(ny+1))
    std::vector<float> u, v;
    std::vector<float> u0, v0; 

    // Pressure and divergence at cell centers (size: nx*ny)
    std::vector<float> pressure;
    std::vector<float> divergence;

    // Marker particles for tracking the fluid.
    std::vector<Particle> particles;

    int idxU(int i, int j) const; // For u (i in [0, nx], j in [0, ny-1])
    int idxV(int i, int j) const; // For v (i in [0, nx-1], j in [0, ny])
    int idxP(int i, int j) const; // For pressure (i in [0, nx-1], j in [0, ny-1])
    int idxDiv(int i, int j) const; // For divergence (i in [0, nx-1], j in [0, ny-1])

    // Simulation steps.
    void addForces();
    void advectVelocities();
    void computeDivergence();
    void solvePressure();
    void subtractPressureGradient();
    void advectParticles();
    
    // Bilinear interpolation for velocities.
    float interpolateU(float x, float y);
    float interpolateV(float x, float y);

    // Simple boundary conditions.
    void setBoundary();
};

} // namespace FluidSim

#endif // FLUID_SIMULATOR_HPP
