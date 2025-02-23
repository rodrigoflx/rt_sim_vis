#include "../include/FluidSimulator.hpp"
#include <iostream>

using namespace FluidSim;

int main() {
    // Create a FluidSimulator instance with a 64x64 grid.
    FluidSimulator simulator(64, 64, 0.1f, 1.0f);

    constexpr int steps = 100;
    for (int i = 0; i < steps; ++i) {
        simulator.step();
        if (i % 10 == 0) {
            const auto& particles = simulator.getParticles();
            if (!particles.empty()) {
                std::cout << "Step " << i << ", Particle[0]: ("
                          << particles[0].x << ", "
                          << particles[0].y << ")\n";
            }
        }
    }

    return 0;
}
