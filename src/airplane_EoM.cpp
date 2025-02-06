#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <array>
#include <vector>

#include "rclcpp/rclcpp.hpp"
#include "geometry_msgs/msg/twist.hpp"
#include "geometry_msgs/msg/pose_stamped.hpp"

#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>

using namespace std::chrono_literals;
using boost::numeric::odeint::runge_kutta4;
using state_type = std::array<double, 12>; // 0-3: longitudinal, 4-8: lateral, 9-11: position

// 定数パラメータ
namespace Param
{
    // 基本定数
    constexpr double g = 9.81;
    constexpr double U0 = 4.0;
    constexpr double W0 = 0;
    constexpr double theta0 = 0.0;

    // 縦の安定微係数: From ChatGPT
    constexpr double Xu = -0.06;
    constexpr double Zu = -0.5;
    constexpr double Mu = -0.025;

    constexpr double Xa = -0.75;
    constexpr double Za = -3.5;
    constexpr double Ma = -1.5;
    constexpr double Madot = -2.0;

    constexpr double Xq = 0;
    constexpr double Zq = -3.5;
    constexpr double Mq = -3.5;

    constexpr double X_deltat = 0.75;
    constexpr double Z_deltae = -2.0;
    constexpr double Z_deltat = 0.0;
    constexpr double M_deltae = -3.5;
    constexpr double M_deltat = 0;

    constexpr double Yb = -1.25;
    constexpr double Lb = -0.6;
    constexpr double Nb = 0.25;

    constexpr double Yp = 0.0;
    constexpr double Lp = -1.0;
    constexpr double Np = 0.15;

    constexpr double Yr = 0.6;
    constexpr double Lr = 0.25;
    constexpr double Nr = -1.0;

    constexpr double L_deltaa = 2.0;
    constexpr double N_deltaa = 0.15;

    constexpr double Y_deltar = 1.5;
    constexpr double L_deltar = -1.0;
    constexpr double N_deltar = -2.0;
} // namespace Param

// 回転行列 R = Rz(psi)*Ry(theta)*Rx(phi)
Eigen::Matrix3d rotation_matrix(double psi, double theta, double phi)
{
    Eigen::Matrix3d Rz, Ry, Rx;
    Rz << std::cos(psi), -std::sin(psi), 0,
        std::sin(psi), std::cos(psi), 0,
        0, 0, 1;
    Ry << std::cos(theta), 0, std::sin(theta),
        0, 1, 0,
        -std::sin(theta), 0, std::cos(theta);
    Rx << 1, 0, 0,
        0, std::cos(phi), -std::sin(phi),
        0, std::sin(phi), std::cos(phi);
    return Rz * Ry * Rx;
}

// 運動方程式クラス（Boost Odeint 用）
class AircraftDynamics
{
public:
    AircraftDynamics()
    {
        // A_long (4x4)
        A_long_.resize(4, 4);
        A_long_ << Param::Xu, Param::Xa, -Param::W0, -Param::g * std::cos(Param::theta0),
            Param::Zu / Param::U0, Param::Za / Param::U0, 1 + Param::Zq / Param::U0, (Param::g / Param::U0) * std::sin(Param::theta0),
            Param::Mu + Param::Madot * (Param::Zu / Param::U0),
            Param::Ma + Param::Madot * (Param::Za / Param::U0),
            Param::Mq + Param::Madot * (1 + Param::Zq / Param::U0),
            (Param::Madot * Param::g / Param::U0) * std::sin(Param::theta0),
            0, 0, 1, 0;

        // B_long (4x2)
        B_long_.resize(4, 2);
        B_long_ << 0, Param::X_deltat,
            Param::Z_deltae / Param::U0, Param::Z_deltat / Param::U0,
            Param::M_deltae + Param::Madot * (Param::Z_deltae / Param::U0),
            Param::M_deltat + Param::Madot * (Param::Z_deltat / Param::U0),
            0, 0;

        // A_lat (5x5)
        A_lat_.resize(5, 5);
        A_lat_ << Param::Yb / Param::U0, (Param::W0 + Param::Yp) / Param::U0, (Param::Yr / Param::U0 - 1), Param::g * std::cos(Param::theta0) / Param::U0, 0,
            Param::Lb, Param::Lp, Param::Lr, 0, 0,
            Param::Nb, Param::Np, Param::Nr, 0, 0,
            0, 1, std::tan(Param::theta0), 0, 0,
            0, 0, 1 / std::cos(Param::theta0), 0, 0;

        // B_lat (5x2)
        B_lat_.resize(5, 2);
        B_lat_ << 0, Param::Y_deltar / Param::U0,
            Param::L_deltaa, Param::L_deltar,
            Param::N_deltaa, Param::N_deltar,
            0, 0,
            0, 0;
    }

    // ROS2 から更新される制御入力
    // u_long: [delta_e, delta_t], u_lat: [delta_a, delta_r]
    void set_control_inputs(double delta_e, double delta_t, double delta_a, double delta_r)
    {
        u_long_[0] = delta_e;
        u_long_[1] = delta_t;
        u_lat_[0] = delta_a;
        u_lat_[1] = delta_r;
        // std::cout << "Control inputs: " << u_long_[0] << ", " << u_long_[1] << ", " << u_lat_[0] << ", " << u_lat_[1] << std::endl;
    }

    // ODE用評価関数 operator():
    // x: state (長さ12)
    // dxdt: derivative
    // t: 時刻
    void operator()(const state_type &x, state_type &dxdt, const double /* t */)
    {
        // 縦運動部分
        Eigen::Vector4d x_long;
        for (int i = 0; i < 4; i++)
        {
            x_long(i) = x[i];
        }
        Eigen::Vector2d u_long;
        u_long << u_long_[0], u_long_[1];
        Eigen::Vector4d dx_long = A_long_ * x_long + B_long_ * u_long;

        // 横運動部分
        Eigen::VectorXd x_lat(5);
        for (int i = 0; i < 5; i++)
        {
            x_lat(i) = x[4 + i];
        }
        Eigen::Vector2d u_lat;
        u_lat << u_lat_[0], u_lat_[1];
        Eigen::VectorXd dx_lat = A_lat_ * x_lat + B_lat_ * u_lat;

        // 機体座標系での速度ベクトル
        double u_b = Param::U0 + x[0];
        double v_b = u_b * std::sin(x[4]); // x[4] を側方角度と仮定
        double w_b = u_b * std::tan(x[1]); // x[1] を迎角変動と仮定
        Eigen::Vector3d vel_b(u_b, v_b, w_b);

        // 全体座標系での速度ベクトル
        double psi = x[8];   // 横方向状態の最後が横ヨー角と仮定
        double theta = x[3]; // 縦のθ（ピッチ）
        double phi = x[7];   // 横ロール角
        Eigen::Vector3d vel_e = rotation_matrix(psi, theta, phi) * vel_b;

        // dxdt の組み立て
        for (int i = 0; i < 4; i++)
        {
            dxdt[i] = dx_long(i);
        }
        for (int i = 0; i < 5; i++)
        {
            dxdt[4 + i] = dx_lat(i);
        }
        for (int i = 0; i < 3; i++)
        {
            dxdt[9 + i] = vel_e(i);
        }
    }

private:
    // 行列メンバ
    Eigen::Matrix<double, 4, 4> A_long_;
    Eigen::Matrix<double, 4, 2> B_long_;
    Eigen::Matrix<double, 5, 5> A_lat_;
    Eigen::Matrix<double, 5, 2> B_lat_;

    // 制御入力（ロング：delta_e, delta_t; ラット：delta_a, delta_r）
    std::array<double, 2> u_long_{0.0, 0.0};
    std::array<double, 2> u_lat_{0.0, 0.0};
};

// ROS2ノード
class DynamicsSimulator : public rclcpp::Node
{
public:
    DynamicsSimulator()
        : Node("dynamics_simulator"), dynamics_()
    {
        // 初期状態 (全要素0とする)
        state_.fill(0.0);
        // サブスクライバー: 制御入力は geometry_msgs::msg::Twist で受信
        // マッピング例:
        // twist.angular.y -> elevator (delta_e)
        // twist.linear.x  -> throttle (delta_t)
        // twist.linear.y  -> aileron (delta_a)
        // twist.angular.z -> rudder (delta_r)
        twist_sub_ = this->create_subscription<geometry_msgs::msg::Twist>(
            "cmd_vel", 10,
            std::bind(&DynamicsSimulator::twist_callback, this, std::placeholders::_1));

        // タイマコールバックでODE１ステップ実行（dt秒）
        timer_ = this->create_wall_timer(
            10ms, std::bind(&DynamicsSimulator::timer_callback, this));

        // 機体位置と姿勢のパブリッシュ
        pose_pub_ = this->create_publisher<geometry_msgs::msg::PoseStamped>("pose", 10);
    }

private:
    void twist_callback(const geometry_msgs::msg::Twist::SharedPtr msg)
    {
        // 制御入力の更新
        double delta_e = msg->linear.x;  // elevator
        double delta_t = msg->linear.y;  // throttle
        double delta_a = msg->linear.y;  // aileron
        double delta_r = msg->angular.z; // rudder

        dynamics_.set_control_inputs(delta_e, delta_t, delta_a, delta_r);
    }

    void timer_callback()
    {
        // 時間刻み
        double dt = 0.01; // 10ms
        // Boost Odeint RK4 による１ステップ
        runge_kutta4<state_type> stepper;
        stepper.do_step(dynamics_, state_, t_, dt);
        t_ += dt;

        // 結果の出力（例として位置ベクトル）
        // RCLCPP_INFO(this->get_logger(), "Position: [%.2f, %.2f, %.2f] ",
        // state_[9], state_[10], state_[11]);

        // 姿勢をRoll, Pitch, YawからQuaternionに変換
        Eigen::Quaterniond q = Eigen::AngleAxisd(state_[7], Eigen::Vector3d::UnitX()) *
                               Eigen::AngleAxisd(state_[3], Eigen::Vector3d::UnitY()) *
                               Eigen::AngleAxisd(state_[8], Eigen::Vector3d::UnitZ());

        // PoseStamped メッセージの組み立て
        auto pose_msg = std::make_shared<geometry_msgs::msg::PoseStamped>();
        pose_msg->header.stamp = this->now();
        pose_msg->pose.position.x = state_[9];
        pose_msg->pose.position.y = state_[10];
        pose_msg->pose.position.z = state_[11];
        pose_msg->pose.orientation.w = q.w();
        pose_msg->pose.orientation.x = q.x();
        pose_msg->pose.orientation.y = q.y();
        pose_msg->pose.orientation.z = q.z();

        pose_pub_->publish(*pose_msg);
    }

    rclcpp::Subscription<geometry_msgs::msg::Twist>::SharedPtr twist_sub_;
    rclcpp::TimerBase::SharedPtr timer_;
    rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pose_pub_;
    AircraftDynamics dynamics_;
    state_type state_;
    double t_{0.0};
};

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<DynamicsSimulator>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}