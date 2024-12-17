#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/point.hpp>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

class AirplaneEOM : public rclcpp::Node
{
public:
    AirplaneEOM() : Node("airplane_eom")
    {
        publisher_ = this->create_publisher<geometry_msgs::msg::Point>("point", 10);
        timer_ = this->create_wall_timer(
            std::chrono::milliseconds(10),
            std::bind(&AirplaneEOM::timer_callback, this));
    }

private:
    void timer_callback()
    {
        integrate(ode_system, state, 0.0, 0.01, 0.001);
        auto message = geometry_msgs::msg::Point();
        message.x = state[0];
        message.y = state[1];
        message.z = state[2];
        publisher_->publish(message);
        // RCLCPP_INFO(this->get_logger(), "Publishing: '%f', '%f', '%f'", message.x, message.y, message.z);
    }

    static void ode_system(const std::vector<double> &x, std::vector<double> &dxdt, double t)
    {
        // ここに微分方程式を定義します
        dxdt[0] = x[1];
        dxdt[1] = -x[0];
        dxdt[2] = 0.0; // 例として単純な微分方程式
    }

    std::vector<double> state = {1.0, 0.0, 0.0}; // 初期状態

    rclcpp::Publisher<geometry_msgs::msg::Point>::SharedPtr publisher_;
    rclcpp::TimerBase::SharedPtr timer_;
};

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<AirplaneEOM>());
    rclcpp::shutdown();
    return 0;
}