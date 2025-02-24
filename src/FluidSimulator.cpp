#include "../include/FluidSimulator.hpp"

#include <algorithm>
#include <memory>
#include <vector>
#include <iostream>

namespace FluidSim {

FluidSimulator::FluidSimulator(int nx, int ny, float dt, float dx)
    : nx(nx), ny(ny), dt(dt), dx(dx),
      u((nx + 1) * ny, 0.0f),
      v(nx * (ny + 1), 0.0f),
      u0((nx + 1) * ny, 0.0f),
      v0(nx * (ny + 1), 0.0f),
      pressure(nx * ny, 0.0f),
      divergence(nx * ny, 0.0f)
{
    // Set default advection solver (semi-Lagrangian, simplified).
    advectionSolverFunc = [this](FluidState& state) {
        // Save the current velocity fields into temporary storage.
        std::copy(state.u.begin(), state.u.end(), state.u0.begin());
        std::copy(state.v.begin(), state.v.end(), state.v0.begin());

        for (int j = 0; j < state.ny; j++) {
            for (int i = 0; i < state.nx + 1; i++) {
                float x = i * state.dx;
                float y = (j + 0.5f) * state.dx;
                float u_val = interpolateU(x, y);
                float v_val = interpolateV(x, y);
                float x_back = x - state.dt * u_val;
                float y_back = y - state.dt * v_val;
                x_back = std::clamp(x_back, 0.0f, state.nx * state.dx);
                y_back = std::clamp(y_back, 0.0f, state.ny * state.dx);
                state.u[i + j * (state.nx + 1)] = interpolateU(x_back, y_back);
            }
        }
    
        // Advect v (vertical velocity).
        for (int j = 0; j < state.ny + 1; j++) {
            for (int i = 0; i < state.nx; i++) {
                float x = (i + 0.5f) * state.dx;
                float y = j * state.dx;
                float u_val = interpolateU(x, y);
                float v_val = interpolateV(x, y);
                float x_back = x - state.dt * u_val;
                float y_back = y - state.dt * v_val;
                x_back = std::clamp(x_back, 0.0f, state.nx * state.dx);
                y_back = std::clamp(y_back, 0.0f, state.ny * state.dx);
                state.v[i + j * state.nx] = interpolateV(x_back, y_back);
            }
        }
    };

    // Set default pressure solver (Jacobi iteration, simplified).
    pressureSolverFunc = [this](FluidState& state) {
        std::fill(state.pressure.begin(), state.pressure.end(), 0.0f);
        const int iterations = 40; // or use state parameter if desired.
        for (int iter = 0; iter < iterations; ++iter) {
            std::vector<float> p_new(state.pressure.size(), 0.0f);
            for (int j = 0; j < state.ny; ++j) {
                for (int i = 0; i < state.nx; ++i) {
                    float p_left  = (i > 0) ? state.pressure[i - 1 + j * state.nx] : 0;
                    float p_right = (i < state.nx - 1) ? state.pressure[i + 1 + j * state.nx] : 0;
                    float p_down  = (j > 0) ? state.pressure[i + (j - 1) * state.nx] : 0;
                    float p_up    = (j < state.ny - 1) ? state.pressure[i + (j + 1) * state.nx] : 0;
                    p_new[i + j * state.nx] = (p_left + p_right + p_down + p_up -
                        state.dx * state.dx * state.divergence[i + j * state.nx]) / 4.0f;
                }
            }
            state.pressure = p_new;
        }
    };

    // Initialize particles in a central square.
    for (int j = ny / 4; j < 3 * ny / 4; ++j) {
        for (int i = nx / 4; i < 3 * nx / 4; ++i) {
            Particle p { (i + 0.5f) * dx, (j + 0.5f) * dx };
            particles.push_back(p);
        }
    }
}

void FluidSimulator::addForces() {
    // Apply gravity to vertical velocity.
    // (Assuming v is defined at cell faces vertically.)
    for (int j = 0; j < ny + 1; ++j) {
        for (int i = 0; i < nx; ++i) {
            int index = i + j * nx;
            if (index < static_cast<int>(v.size()))
                v[index] += dt * gravity;
        }
    }
}

void FluidSimulator::computeDivergence() {
    // Compute divergence at cell centers.
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            float u_right = u[(i + 1) + j * (nx + 1)];
            float u_left  = u[i + j * (nx + 1)];
            float v_top   = v[i + (j + 1) * nx];
            float v_bottom= v[i + j * nx];
            divergence[i + j * nx] = (u_right - u_left + v_top - v_bottom) / dx;
        }
    }
}

void FluidSimulator::subtractPressureGradient() {
    // Adjust horizontal velocity (u).
    for (int j = 0; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            float gradP = (pressure[i + j * nx] - pressure[i - 1 + j * nx]) / dx;
            u[i + j * (nx + 1)] -= dt * gradP;
        }
    }
    // Adjust vertical velocity (v).
    for (int j = 1; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            float gradP = (pressure[i + j * nx] - pressure[i + (j - 1) * nx]) / dx;
            v[i + j * nx] -= dt * gradP;
        }
    }
}

void FluidSimulator::setBoundary() {
    // Enforce zero velocity at the boundaries.
    for (int j = 0; j < ny; ++j) {
        u[0 + j * (nx + 1)] = 0;
        u[nx + j * (nx + 1)] = 0;
    }
    for (int i = 0; i < nx; ++i) {
        v[i + 0 * nx] = 0;
        v[i + ny * nx] = 0;
    }
}

void FluidSimulator::advectParticles() {
    for (auto &p : particles) {
        float u_val = interpolateU(p.x, p.y);
        float v_val = interpolateV(p.x, p.y);
        p.x += dt * u_val;
        p.y += dt * v_val;
        p.x = std::clamp(p.x, 0.0f, nx * dx);
        p.y = std::clamp(p.y, 0.0f, ny * dx);
    }
}

float FluidSimulator::interpolateU(float x, float y) {
    // Convert physical coordinates to grid coordinates for u.
    float fx = x / dx;             // x coordinate in grid units.
    float fy = (y / dx) - 0.5f;      // shift y by 0.5.
    int i0 = static_cast<int>(std::floor(fx));
    int j0 = static_cast<int>(std::floor(fy));
    // Clamp indices.
    if (i0 < 0) { i0 = 0; fx = 0.0f; }
    if (i0 > nx - 1) { i0 = nx - 1; fx = static_cast<float>(nx - 1); }
    if (j0 < 0) { j0 = 0; fy = 0.0f; }
    if (j0 > ny - 2) { j0 = ny - 2; fy = static_cast<float>(ny - 2); }
    float s = fx - i0;
    float t = fy - j0;
    // u is stored in a vector of size (nx+1)*ny with index: i + j*(nx+1).
    int idx00 = i0 + j0 * (nx + 1);
    int idx10 = (i0 + 1) + j0 * (nx + 1);
    int idx01 = i0 + (j0 + 1) * (nx + 1);
    int idx11 = (i0 + 1) + (j0 + 1) * (nx + 1);
    float u00 = u[idx00];
    float u10 = u[idx10];
    float u01 = u[idx01];
    float u11 = u[idx11];
    return (1 - t) * ((1 - s) * u00 + s * u10) + t * ((1 - s) * u01 + s * u11);
}

float FluidSimulator::interpolateV(float x, float y) {
    // Convert physical coordinates to grid coordinates for v.
    float fx = (x / dx) - 0.5f;  // shift x by 0.5.
    float fy = y / dx;
    int i0 = static_cast<int>(std::floor(fx));
    int j0 = static_cast<int>(std::floor(fy));
    // Clamp indices.
    if (i0 < 0) { i0 = 0; fx = 0.0f; }
    if (i0 > nx - 2) { i0 = nx - 2; fx = static_cast<float>(nx - 2); }
    if (j0 < 0) { j0 = 0; fy = 0.0f; }
    if (j0 > ny - 1) { j0 = ny - 1; fy = static_cast<float>(ny - 1); }
    float s = fx - i0;
    float t = fy - j0;
    // v is stored in a vector of size nx*(ny+1) with index: i + j*nx.
    int idx00 = i0 + j0 * nx;
    int idx10 = (i0 + 1) + j0 * nx;
    int idx01 = i0 + (j0 + 1) * nx;
    int idx11 = (i0 + 1) + (j0 + 1) * nx;
    float v00 = v[idx00];
    float v10 = v[idx10];
    float v01 = v[idx01];
    float v11 = v[idx11];
    return (1 - t) * ((1 - s) * v00 + s * v10) + t * ((1 - s) * v01 + s * v11);
}


void FluidSimulator::step() {
    addForces();

    // Create a FluidState view.
    FluidState state { u, v, u0, v0, pressure, divergence, nx, ny, dt, dx };

    if (advectionSolverFunc)
        advectionSolverFunc(state);

    setBoundary();
    computeDivergence();

    if (pressureSolverFunc)
        pressureSolverFunc(state);

    subtractPressureGradient();
    setBoundary();
    advectParticles();
}

const std::vector<Particle>& FluidSimulator::getParticles() const {
    return particles;
}

} // namespace FluidSim
