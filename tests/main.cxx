
#include "EventShapes/EventShapes.h"

#include <Eigen/Dense>

#include <iostream>
#include <vector>

using namespace Eigen;

int test_example_eigen() {
    std::cout << "Running Eigen example" << std::endl;

    // Test eigen
    MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);
    std::cout << m << std::endl;

    std::cout << std::endl;
    return 0;
}

int test_example_eventshapes(unsigned int ndims = 2) {
    std::cout << "Running example Event Shapes tests" << std::endl;

    // Construct example input vectors
    std::vector<std::vector<float>> vs;

    EventShapes es;

    if (ndims == 3) {
        float epsilon = 0.05;
        vs.push_back({1., epsilon, 0.});
        vs.push_back({0., 1., epsilon});
        vs.push_back({epsilon, 0., 1.});

        es = EventShapes(vs, 3);

    } else {
        vs.push_back({0.5, 1.});
        vs.push_back({1., -2.5});
        vs.push_back({0., 1.5});
        vs.push_back({-1., 2.});

        es = EventShapes(vs, 2);

    }

    // Calculate
    es.calc_all();
    std::cout << "Thrust: " << es.get_thrust() << std::endl;
    std::cout << "Thrust major: " << es.get_thrust_major() << std::endl;
    std::cout << "Thrust minor: " << es.get_thrust_minor() << std::endl;
    std::cout << "Oblateness: " << es.get_oblateness() << std::endl;
    std::cout << "Broadening: " << es.get_broadening() << std::endl;
    std::cout << "S: " << es.get_lin_spher_S() << std::endl;
    std::cout << "A: " << es.get_lin_spher_A() << std::endl;
    std::cout << "C: " << es.get_lin_spher_C() << std::endl;
    std::cout << "D: " << es.get_lin_spher_D() << std::endl;


    std::cout << std::endl;
    return 0;
}


int test_lvs() {
    std::cout << "Running LVS tests" << std::endl;

    EventShapes es;

    std::vector<std::vector<float>> vs;
    vs.push_back({0., 0., 1.});
    vs.push_back({0., 0.5, 1.});
    vs.push_back({4., 1.2, 0.1});
    vs.push_back({2., 1., -3.});
    vs.push_back({-1, -2, 4});
    es = EventShapes(vs, 3);

    es.calc_all();
    std::cout << "Trad thrust: " << es.get_thrust() << std::endl;

    es.lvs_t();

    return 0;
}

int main() {

    std::cout << "Hello, World!" << std::endl;
    std::cout << "Running package event-shapes tests" << std::endl;

    test_example_eigen();

    test_example_eventshapes(3);
    test_example_eventshapes(2);

    test_lvs();

    return 0;
}
