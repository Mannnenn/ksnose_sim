#include <rclcpp/rclcpp.hpp>
#include <boost/version.hpp>

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    auto node = rclcpp::Node::make_shared("version_info_node");

    RCLCPP_INFO(node->get_logger(), "boostバージョン: %d", BOOST_VERSION);
    RCLCPP_INFO(node->get_logger(), "boostライブラリバージョン: %s", BOOST_LIB_VERSION);

    rclcpp::shutdown();
    return 0;
}
