#ifndef AIRPLANE_EOM_HPP
#define AIRPLANE_EOM_HPP

#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/point.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>

using namespace boost::numeric::odeint;
using namespace boost::numeric::ublas;

typedef vector<double> state_type;

class AirplaneEOM : public rclcpp::Node
{
public:
    AirplaneEOM();

private:
    void timer_callback();
    void aircraft_dynamics(const state_type &x, state_type &dxdt, double t, const vector<vector<double>> &A_long, const vector<vector<double>> &B_long, const vector<vector<double>> &A_lat, const vector<vector<double>> &B_lat, const vector<double> &u_long, const vector<double> &u_lat);
    void initialize_matrices();
    matrix<double> rotation_matrix(double psi, double theta, double phi);

    std::vector<double> state;
    rclcpp::Publisher<geometry_msgs::msg::Point>::SharedPtr publisher_;
    rclcpp::TimerBase::SharedPtr timer_;

    matrix<double> A_long;
    matrix<double> B_long;
    matrix<double> A_lat;
    matrix<double> B_lat;

    static constexpr double Xu = -0.02;
    static constexpr double Xa = 0.1;
    static constexpr double W0 = 10.0;
    static constexpr double g = 9.81;
    static constexpr double theta0 = 0.1;
    static constexpr double Zu = -0.1;
    static constexpr double U0 = 100.0;
    static constexpr double Za = -0.2;
    static constexpr double Zq = -0.3;
    static constexpr double Mu = 0.01;
    static constexpr double Madot = 0.02;
    static constexpr double Ma = 0.03;
    static constexpr double Mq = 0.04;
    static constexpr double X_deltat = 0.05;
    static constexpr double Z_deltae = 0.06;
    static constexpr double Z_deltat = 0.07;
    static constexpr double M_deltae = 0.08;
    static constexpr double M_deltat = 0.09;
    static constexpr double Yb = 0.1;
    static constexpr double Yp = 0.2;
    static constexpr double Yr = 0.3;
    static constexpr double Lb = 0.4;
    static constexpr double Lp = 0.5;
    static constexpr double Lr = 0.6;
    static constexpr double Nb = 0.7;
    static constexpr double Np = 0.8;
    static constexpr double Nr = 0.9;
    static constexpr double Y_deltar = 1.0;
    static constexpr double L_deltaa = 1.1;
    static constexpr double L_deltar = 1.2;
    static constexpr double N_deltaa = 1.3;
    static constexpr double N_deltar = 1.4;

    // 初期状態 (縦と横・方向、位置を結合)
    static constexpr state_type x0 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    static constexpr double t_start = 0.0;
    static constexpr double t_end = 1000.0;
    static constexpr double dt = 0.01;

    // 初期制御入力 (縦と横・方向)
    vector<double> u0_long; // 例: elevator=0, throttle=0
    vector<double> u0_lat;  // 例: aileron=0, rudder=0

    // シミュレーションの実行
    runge_kutta_dopri5<state_type> stepper;
    double t;
    state_type x;
};

#endif // AIRPLANE_EOM_HPP