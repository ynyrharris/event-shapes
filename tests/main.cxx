
#include "EventShapes/EventShapes.h"

#include <Eigen/Dense>

#include <iostream>
#include <vector>
#include <chrono>

using namespace Eigen;

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


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

    auto t1 = high_resolution_clock::now();
    es.calc_all();
    auto t2 = high_resolution_clock::now();
    std::cout << "Trad thrust: " << es.get_thrust() << std::endl;
    auto ms = duration<double, std::milli>(t2 - t1);
    std::cout << "Trad thrust took " << ms.count() << "ms" << std::endl;

    t1= high_resolution_clock::now();
    es.lvs_t();
    t2 = high_resolution_clock::now();
    ms = duration<double, std::milli>(t2 - t1);
    std::cout << "LVS thrust took " << ms.count() << "ms" << std::endl;



    float epsilon = 0.05;
    vs.clear();
    vs.push_back({1., epsilon, 0.});
    vs.push_back({0., 1., epsilon});
    vs.push_back({epsilon, 0., 1.});
    es = EventShapes(vs, 3);

    t1 = high_resolution_clock::now();
    es.calc_all();
    t2 = high_resolution_clock::now();
    std::cout << "Trad thrust: " << es.get_thrust() << std::endl;
    ms = duration<double, std::milli>(t2 - t1);
    std::cout << "Trad thrust took " << ms.count() << "ms" << std::endl;

    t1= high_resolution_clock::now();
    es.lvs_t();
    t2 = high_resolution_clock::now();
    ms = duration<double, std::milli>(t2 - t1);
    std::cout << "LVS thrust took " << ms.count() << "ms" << std::endl;

    vs.clear();
    for (unsigned int i = 0; i < 40; i++) {
        Eigen::Vector3f v = Eigen::Vector3f::Random();
        vs.push_back({v.x(), v.y(), v.z()});
    }

    es = EventShapes(vs, 3);

    t1 = high_resolution_clock::now();
    es.calc_all();
    t2 = high_resolution_clock::now();
    std::cout << "Trad thrust: " << es.get_thrust() << std::endl;
    ms = duration<double, std::milli>(t2 - t1);
    std::cout << "Trad thrust took " << ms.count() << "ms" << std::endl;

    t1= high_resolution_clock::now();
    es.lvs_t();
    t2 = high_resolution_clock::now();
    ms = duration<double, std::milli>(t2 - t1);
    std::cout << "LVS thrust took " << ms.count() << "ms" << std::endl;

    return 0;
}


int benchmark_lvs() {
    std::cout << "Benchmarking LVS agorithm" << std::endl;

    // Perform benchmarking
    std::vector<std::vector<float>> vs;

    std::vector<unsigned int> lenv;
    std::vector<double> lvs_times;
    std::vector<double> trad_times;
    std::vector<double> lvs_ts;
    std::vector<double> trad_ts;

    for (unsigned int i = 1; i < 30; i++) {
        unsigned int N = i * 2.;

        vs.clear();
        vs.reserve(N);

        for (unsigned int j = 0; j < N; j++) {
            Eigen::Vector3f v = Eigen::Vector3f::Random();
            vs.push_back({v.x(), v.y(), v.z()});
        }

        // Initialise EventShapes object
        EventShapes es = EventShapes(vs, 3);

        // Perform Trad calculation
        auto t_trad_0 = high_resolution_clock::now();
        es.calc_all();
        auto t_trad_1 = high_resolution_clock::now();
        auto t_trad_obj = duration<double, std::milli>(t_trad_1 - t_trad_0);
        double t_trad = t_trad_obj.count();

        // Perform LVS calculation
        auto t_lvs_0 = high_resolution_clock::now();
        es.lvs_t();
        auto t_lvs_1 = high_resolution_clock::now();
        auto t_lvs_obj = duration<double, std::milli>(t_lvs_1 - t_lvs_0);
        double t_lvs = t_lvs_obj.count();

        // Print benchmark parameters:
        // N input vectors  |  V.size()  |  t_lvs  |  t_trad  |  lvs_t  | trad_t
        std::cout << N;
        std::cout << "  |  " << es.m_lvs_lenV;
        std::cout << "  |  " << t_lvs;
        std::cout << "  |  " << t_trad;
        std::cout << "  |  " << es.m_lvs_t;
        std::cout << "  |  " << es.get_thrust();
        std::cout << std::endl;

        lenv.push_back(es.m_lvs_lenV);
        lvs_times.push_back(t_lvs);
        trad_times.push_back(t_trad);
        lvs_ts.push_back(es.m_lvs_t);
        trad_ts.push_back(es.get_thrust());
    }

    std::cout << "lenV: ";
    for (unsigned int a : lenv) {
        std::cout << a << ", ";
    } std::cout << std::endl;

    std::cout << "lvs times: ";
    for (double t : lvs_times) {
        std::cout << t << ", ";
    } std::cout << std::endl;

    std::cout << "trad times: ";
    for (double t : trad_times) {
        std::cout << t << ", ";
    } std::cout << std::endl;

    std::cout << "lvs ts: ";
    for (double t : lvs_ts) {
        std::cout << t << ", ";
    } std::cout << std::endl;

    std::cout << "trad ts: ";
    for (double t : trad_ts) {
        std::cout << t << ", ";
    } std::cout << std::endl;

    return 0;
}

int main() {

    std::cout << "Hello, World!" << std::endl;
    std::cout << "Running package event-shapes tests" << std::endl;

    test_example_eigen();

    test_example_eventshapes(3);
    test_example_eventshapes(2);

    // test_lvs();

    // benchmark_lvs();

    return 0;
}
