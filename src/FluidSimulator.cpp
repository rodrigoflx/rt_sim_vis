#include "../include/FluidSimulator.hpp"
#include <algorithm>
#include <cmath>
#include <cassert>

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
    // Initialize particles in a central square.
    for (int j = ny / 4; j < 3 * ny / 4; ++j) {
        for (int i = nx / 4; i < 3 * nx / 4; ++i) {
            Particle p{ (i + 0.5f) * dx, (j + 0.5f) * dx };
            particles.push_back(p);
        }
    }
}

int FluidSimulator::idxU(int i, int j) const {
    return i + j * (nx + 1);
}

int FluidSimulator::idxV(int i, int j) const {
    return i + j * nx;
}

int FluidSimulator::idxP(int i, int j) const {
    return i + j * nx;
}

int FluidSimulator::idxDiv(int i, int j) const {
    return i + j * nx;
}

void FluidSimulator::addForces() {
    // Apply gravity to vertical velocity v.
    for (int j = 0; j < ny + 1; j++) {
        for (int i = 0; i < nx; i++) {
            int index = idxV(i, j);
            v[index] += dt * gravity;
        }
    }
}

void FluidSimulator::advectVelocities() {
    // Save current velocities.
    std::copy(u.begin(), u.end(), u0.begin());
    std::copy(v.begin(), v.end(), v0.begin());

    // Advect u (horizontal velocity) at (i*dx, (j+0.5)*dx)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx + 1; i++) {
            float x = i * dx;
            float y = (j + 0.5f) * dx;
            float u_val = interpolateU(x, y);
            float v_val = interpolateV(x, y);
            float x_back = x - dt * u_val;
            float y_back = y - dt * v_val;
            // Clamp traced position using std::clamp (C++17).
            x_back = std::clamp(x_back, 0.0f, nx * dx);
            y_back = std::clamp(y_back, 0.0f, ny * dx);
            u[idxU(i, j)] = interpolateU(x_back, y_back);
        }
    }

    // Advect v (vertical velocity) at ((i+0.5)*dx, j*dx)
    for (int j = 0; j < ny + 1; j++) {
        for (int i = 0; i < nx; i++) {
            float x = (i + 0.5f) * dx;
            float y = j * dx;
            float u_val = interpolateU(x, y);
            float v_val = interpolateV(x, y);
            float x_back = x - dt * u_val;
            float y_back = y - dt * v_val;
            x_back = std::clamp(x_back, 0.0f, nx * dx);
            y_back = std::clamp(y_back, 0.0f, ny * dx);
            v[idxV(i, j)] = interpolateV(x_back, y_back);
        }
    }
}

void FluidSimulator::computeDivergence() {
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            float u_right = u[idxU(i + 1, j)];
            float u_left  = u[idxU(i, j)];
            float v_top   = v[idxV(i, j + 1)];
            float v_bottom = v[idxV(i, j)];
            divergence[idxDiv(i, j)] = (u_right - u_left + v_top - v_bottom) / dx;
        }
    }
}

void FluidSimulator::solvePressure() {
    std::fill(pressure.begin(), pressure.end(), 0.0f);
    for (int iter = 0; iter < pressureIterations; iter++) {
        std::vector<float> p_new(pressure.size(), 0.0f);
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                float p_left  = (i > 0) ? pressure[idxP(i - 1, j)] : 0;
                float p_right = (i < nx - 1) ? pressure[idxP(i + 1, j)] : 0;
                float p_down  = (j > 0) ? pressure[idxP(i, j - 1)] : 0;
                float p_up    = (j < ny - 1) ? pressure[idxP(i, j + 1)] : 0;
                p_new[idxP(i, j)] = (p_left + p_right + p_down + p_up -
                                      dx * dx * divergence[idxDiv(i, j)]) / 4.0f;
            }
        }
        pressure = p_new;
    }
}

void FluidSimulator::subtractPressureGradient() {
    // Update horizontal velocity u.
    for (int j = 0; j < ny; j++) {
        for (int i = 1; i < nx; i++) {
            float gradP = (pressure[idxP(i, j)] - pressure[idxP(i - 1, j)]) / dx;
            u[idxU(i, j)] -= dt * gradP;
        }
    }
    // Update vertical velocity v.
    for (int j = 1; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            float gradP = (pressure[idxP(i, j)] - pressure[idxP(i, j - 1)]) / dx;
            v[idxV(i, j)] -= dt * gradP;
        }
    }
}

void FluidSimulator::advectParticles() {
    for (auto& p : particles) {
        float u_val = interpolateU(p.x, p.y);
        float v_val = interpolateV(p.x, p.y);
        p.x += dt * u_val;
        p.y += dt * v_val;
        p.x = std::clamp(p.x, 0.0f, nx * dx);
        p.y = std::clamp(p.y, 0.0f, ny * dx);
    }
}

float FluidSimulator::interpolateU(float x, float y) {
    // Map physical coordinates to grid coordinates for u.
    float fx = x / dx;
    float fy = (y / dx) - 0.5f;
    int i = std::clamp(static_cast<int>(fx), 0, nx);
    int j = std::clamp(static_cast<int>(fy), 0, ny - 1);
    float s = fx - i;
    float t = fy - j;
    float u00 = u0[idxU(i, j)];
    float u10 = (i < nx) ? u0[idxU(i + 1, j)] : u00;
    float u01 = u0[idxU(i, std::min(j + 1, ny - 1))];
    float u11 = (i < nx) ? u0[idxU(i + 1, std::min(j + 1, ny - 1))] : u01;
    return (1 - t) * ((1 - s) * u00 + s * u10) + t * ((1 - s) * u01 + s * u11);
}

float FluidSimulator::interpolateV(float x, float y) {
    // Map physical coordinates to grid coordinates for v.
    float fx = (x / dx) - 0.5f;
    float fy = y / dx;
    int i = std::clamp(static_cast<int>(fx), 0, nx - 1);
    int j = std::clamp(static_cast<int>(fy), 0, ny);
    float s = fx - i;
    float t = fy - j;
    float v00 = v0[idxV(i, j)];
    float v10 = (i < nx - 1) ? v0[idxV(i + 1, j)] : v00;
    float v01 = (j < ny) ? v0[idxV(i, j + 1)] : v00;
    float v11 = (i < nx - 1 && j < ny) ? v0[idxV(i + 1, j + 1)] : v01;
    return (1 - t) * ((1 - s) * v00 + s * v10) + t * ((1 - s) * v01 + s * v11);
}

void FluidSimulator::setBoundary() {
    // Zero-out velocities at the domain boundaries.
    for (int j = 0; j < ny; j++) {
        u[idxU(0, j)] = 0;
        u[idxU(nx, j)] = 0;
    }
    for (int i = 0; i < nx; i++) {
        v[idxV(i, 0)] = 0;
        v[idxV(i, ny)] = 0;
    }
}

void FluidSimulator::step() {
    addForces();
    advectVelocities();
    setBoundary();
    computeDivergence();
    solvePressure();
    subtractPressureGradient();
    setBoundary();
    advectParticles();
}

const std::vector<Particle>& FluidSimulator::getParticles() const {
    return particles;
}

} // namespace FluidSim
