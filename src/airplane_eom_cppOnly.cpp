#include "nokolat2024_sim/airplane_eom.hpp"
#include <rclcpp/rclcpp.hpp>
#include <iostream>
#include <chrono>
#include <boost/numeric/odeint.hpp>
#include <cmath>

using namespace boost::numeric::odeint;
using namespace boost::numeric::ublas;

// --- パラメータ ---
const double g = 9.81;
const double U0 = 293.8;
const double W0 = 0;
const double theta0 = 0.0;

const double Xu = -0.0215;
const double Zu = -0.227;
const double Mu = 0;

const double Xa = 14.7;
const double Za = -236;
const double Ma = -3.78;
const double Madot = -0.28;

const double Xq = 0;
const double Zq = -5.76;
const double Mq = -0.992;

const double X_deltat = 0;
const double Z_deltae = -12.9;
const double Z_deltat = 0;
const double M_deltae = -2.48;
const double M_deltat = 0;

const double Yb = -45.4;
const double Lb = -1.71;
const double Nb = 0.986;

const double Yp = 0.716;
const double Lp = -0.962;
const double Np = -0.0632;

const double Yr = 2.66;
const double Lr = 0.271;
const double Nr = -0.215;

const double L_deltaa = 1.72;
const double N_deltaa = -0.0436;

const double Y_deltar = 9.17;
const double L_deltar = 0.244;
const double N_deltar = -0.666;

matrix<double> A_long(4, 4);
matrix<double> B_long(4, 2);
matrix<double> A_lat(5, 5);
matrix<double> B_lat(5, 2);

void initialize_matrices()
{
    A_long(0, 0) = Xu;
    A_long(0, 1) = Xa;
    A_long(0, 2) = -W0;
    A_long(0, 3) = -g * cos(theta0);
    A_long(1, 0) = Zu / U0;
    A_long(1, 1) = Za / U0;
    A_long(1, 2) = 1 + Zq / U0;
    A_long(1, 3) = (g / U0) * sin(theta0);
    A_long(2, 0) = Mu + Madot * (Zu / U0);
    A_long(2, 1) = Ma + Madot * (Za / U0);
    A_long(2, 2) = Mq + Madot * (1 + Zq / U0);
    A_long(2, 3) = (Madot * g / U0) * sin(theta0);
    A_long(3, 0) = 0;
    A_long(3, 1) = 0;
    A_long(3, 2) = 1;
    A_long(3, 3) = 0;

    B_long(0, 0) = 0;
    B_long(0, 1) = X_deltat;
    B_long(1, 0) = Z_deltae / U0;
    B_long(1, 1) = Z_deltat / U0;
    B_long(2, 0) = M_deltae + Madot * (Z_deltae / U0);
    B_long(2, 1) = M_deltat + Madot * (Z_deltat / U0);
    B_long(3, 0) = 0;
    B_long(3, 1) = 0;

    A_lat(0, 0) = Yb / U0;
    A_lat(0, 1) = (W0 + Yp) / U0;
    A_lat(0, 2) = (Yr / U0 - 1);
    A_lat(0, 3) = g * cos(theta0) / U0;
    A_lat(0, 4) = 0;
    A_lat(1, 0) = Lb;
    A_lat(1, 1) = Lp;
    A_lat(1, 2) = Lr;
    A_lat(1, 3) = 0;
    A_lat(1, 4) = 0;
    A_lat(2, 0) = Nb;
    A_lat(2, 1) = Np;
    A_lat(2, 2) = Nr;
    A_lat(2, 3) = 0;
    A_lat(2, 4) = 0;
    A_lat(3, 0) = 0;
    A_lat(3, 1) = 1;
    A_lat(3, 2) = tan(theta0);
    A_lat(3, 3) = 0;
    A_lat(3, 4) = 0;
    A_lat(4, 0) = 0;
    A_lat(4, 1) = 0;
    A_lat(4, 2) = 1 / cos(theta0);
    A_lat(4, 3) = 0;
    A_lat(4, 4) = 0;

    B_lat(0, 0) = 0;
    B_lat(0, 1) = Y_deltar / U0;
    B_lat(1, 0) = L_deltaa;
    B_lat(1, 1) = L_deltar;
    B_lat(2, 0) = N_deltaa;
    B_lat(2, 1) = N_deltar;
    B_lat(3, 0) = 0;
    B_lat(3, 1) = 0;
    B_lat(4, 0) = 0;
    B_lat(4, 1) = 0;
}

matrix<double> rotation_matrix(double psi, double theta, double phi)
{
    matrix<double> R_z(3, 3), R_y(3, 3), R_x(3, 3), R(3, 3);

    R_z(0, 0) = cos(psi);
    R_z(0, 1) = -sin(psi);
    R_z(0, 2) = 0;
    R_z(1, 0) = sin(psi);
    R_z(1, 1) = cos(psi);
    R_z(1, 2) = 0;
    R_z(2, 0) = 0;
    R_z(2, 1) = 0;
    R_z(2, 2) = 1;

    R_y(0, 0) = cos(theta);
    R_y(0, 1) = 0;
    R_y(0, 2) = sin(theta);
    R_y(1, 0) = 0;
    R_y(1, 1) = 1;
    R_y(1, 2) = 0;
    R_y(2, 0) = -sin(theta);
    R_y(2, 1) = 0;
    R_y(2, 2) = cos(theta);

    R_x(0, 0) = 1;
    R_x(0, 1) = 0;
    R_x(0, 2) = 0;
    R_x(1, 0) = 0;
    R_x(1, 1) = cos(phi);
    R_x(1, 2) = -sin(phi);
    R_x(2, 0) = 0;
    R_x(2, 1) = sin(phi);
    R_x(2, 2) = cos(phi);

    R = prod(R_z, R_y);
    R = prod(R, R_x);

    return R;
}

void aircraft_dynamics(const std::vector<double> &x, std::vector<double> &dxdt, double t)
{
    vector<double> x_long(4), x_lat(5), u_long(2), u_lat(2);
    for (int i = 0; i < 4; ++i)
        x_long(i) = x[i];
    for (int i = 0; i < 5; ++i)
        x_lat(i) = x[i + 4];
    u_long(0) = 0.01;
    u_long(1) = 0;
    u_lat(0) = 0.01;
    u_lat(1) = 0;

    vector<double> dxdt_long = prod(A_long, x_long) + prod(B_long, u_long);
    vector<double> dxdt_lat = prod(A_lat, x_lat) + prod(B_lat, u_lat);

    double u_b = U0 + x[0];
    double v_b = u_b * sin(x[4]);
    double w_b = u_b * tan(x[1]);
    vector<double> vel_b(3);
    vel_b(0) = u_b;
    vel_b(1) = v_b;
    vel_b(2) = w_b;

    double psi = x[8];
    double theta = x[3];
    double phi = x[7];
    matrix<double> R = rotation_matrix(psi, theta, phi);
    vector<double> vel_e = prod(R, vel_b);

    for (int i = 0; i < 4; ++i)
        dxdt[i] = dxdt_long(i);
    for (int i = 0; i < 5; ++i)
        dxdt[i + 4] = dxdt_lat(i);
    for (int i = 0; i < 3; ++i)
        dxdt[i + 9] = vel_e(i);
}

int main()
{
    using namespace std;
    chrono::system_clock::time_point start, end;

    start = chrono::system_clock::now();

    initialize_matrices();

    std::vector<double> x0(12, 0.0);
    std::vector<double> t(200);
    for (int i = 0; i < 200; ++i)
        t[i] = i;

    runge_kutta_dopri5<std::vector<double>> stepper;
    std::vector<std::vector<double>> x_history;

    for (int i = 0; i < 1000; ++i)
    {
        x_history.clear();
        integrate_times(stepper, aircraft_dynamics, x0, t.begin(), t.end(), 0.1, [&](const std::vector<double> &x, double t)
                        { x_history.push_back(x); });
        // std::cout << i << std::endl;
    }

    // for (const auto &x : x_history)
    // {
    //     std::cout << x[9] << " " << x[10] << " " << x[11] << std::endl;
    // }

    end = chrono::system_clock::now();

    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("time %lf[ms]\n", time);

    return 0;
}