#ifndef FLUID_SIMULATOR_HPP
#define FLUID_SIMULATOR_HPP

#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>

namespace FluidSim {

// A simple particle structure.
struct Particle {
    float x, y;
};

// FluidState exposes just the data needed for solver routines.
struct FluidState {
    std::vector<float>& u;          // horizontal velocity field
    std::vector<float>& v;          // vertical velocity field
    std::vector<float>& u0;         // temporary horizontal velocity field
    std::vector<float>& v0;         // temporary vertical velocity field
    std::vector<float>& pressure;   // pressure field
    std::vector<float>& divergence; // divergence field
    int nx, ny;                     // grid dimensions
    float dt, dx;                   // time step and cell size
};

class FluidSimulator {
public:
    FluidSimulator(int nx, int ny, float dt, float dx);

    // Set the solver callbacks.
    void setAdvectionSolver(std::function<void(FluidState&)> solver);
    void setPressureSolver(std::function<void(FluidState&)> solver);

    // Advance the simulation one time step.
    void step();

    // Accessor for particles.
    const std::vector<Particle>& getParticles() const;

private:
    int nx, ny;
    float dt, dx;
    static constexpr int pressureIterations = 40;
    static constexpr float gravity = -9.8f;

    // Fluid grid data.
    std::vector<float> u, v;
    std::vector<float> u0, v0;
    std::vector<float> pressure;
    std::vector<float> divergence;

    // Marker particles.
    std::vector<Particle> particles;

    // Solver callbacks.
    std::function<void(FluidState&)> advectionSolverFunc;
    std::function<void(FluidState&)> pressureSolverFunc;

    // Simulation sub-steps.
    void addForces();
    void computeDivergence();
    void subtractPressureGradient();
    void setBoundary();
    void advectParticles();

    // Bilinear interpolation functions.
    float interpolateU(float x, float y);
    float interpolateV(float x, float y);
};

} // namespace FluidSim

#endif // FLUID_SIMULATOR_HPP
