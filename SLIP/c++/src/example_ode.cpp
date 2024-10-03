#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <chrono>

using namespace boost::numeric::odeint;

// Define the state type (a vector for the system [r, theta, r_dot, theta_dot])
typedef std::vector<double> state_type;

// Structure to hold system parameters
struct Params {
    double m;   // mass
    double g;   // gravitational acceleration
    double l0;  // natural leg length
};

// Function to compute spring stiffness (you need to define this)
double spring_stiffness(double t, const state_type &x_polar, const Params &params) {
    // Example placeholder stiffness function (you can modify it as needed)
    // Replace this with your actual spring stiffness computation
    double k = 1000.0;  // example constant stiffness
    return k;
}

// SLIP ground dynamics function
void dynamics_g(const state_type &x_polar, state_type &xdot, double t, const Params &params) {
    // Unpack parameters
    double m = params.m;
    double g = params.g;
    double l0 = params.l0;

    // Get spring stiffness
    double k = spring_stiffness(t, x_polar, params);

    // Polar state: x_polar = [r, theta, r_dot, theta_dot]
    double r = x_polar[0];
    double theta = x_polar[1];
    double r_dot = x_polar[2];
    double theta_dot = x_polar[3];

    // Dynamics equations
    xdot[0] = r_dot;
    xdot[1] = theta_dot;
    xdot[2] = r * theta_dot * theta_dot - g * cos(theta) + (k / m) * (l0 - r);
    xdot[3] = -(2 / r) * r_dot * theta_dot + (g / r) * sin(theta);
}

int main() {

    //  to measure time 
    auto start = std::chrono::high_resolution_clock::now();

    // Initial conditions: [r, theta, r_dot, theta_dot]
    state_type x_polar(4);
    x_polar[0] = 1.0;   // r (initial leg length)
    x_polar[1] = M_PI / 4.0;  // theta (initial leg angle)
    x_polar[2] = 0.0;   // r_dot (initial radial velocity)
    x_polar[3] = 0.0;   // theta_dot (initial angular velocity)

    // Define system parameters
    Params params;
    params.m = 80.0;  // mass (in kg)
    params.g = 9.81;  // gravitational acceleration (in m/s^2)
    params.l0 = 1.0;  // natural leg length (in meters)

    // Integrate the system using odeint
    double t_start = 0.0;
    double t_end = 2.0;
    double dt = 0.01;  // Time step

    // Define the stepper (Runge-Kutta 4th order)
    runge_kutta4<state_type> stepper;

    // Integrate the system over time and output the results
    for (double t = t_start; t < t_end; t += dt) {
        // Print current state (time, r, theta, r_dot, theta_dot)
        // std::cout << "t = " << t
        //           << " r = " << x_polar[0]
        //           << " theta = " << x_polar[1]
        //           << " r_dot = " << x_polar[2]
        //           << " theta_dot = " << x_polar[3] << std::endl;

        // Perform the integration step
        stepper.do_step([&](const state_type &x, state_type &dxdt, double t) {
            dynamics_g(x, dxdt, t, params);
        }, x_polar, t, dt);
    }

    // stopwatch 
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time taken by function: "
              << duration.count() << " microseconds" << std::endl;

    return 0;
}