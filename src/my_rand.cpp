#include "my_rand.hpp"

using v3d = Eigen::Vector3d;

namespace my_rand {

double rand(double max, double min) {
    if (min > max) std::swap(min, max);
    double width = max - min;
    return (double) min + max * std::rand() / RAND_MAX;
};

int randint(int sup, int min) {
    if (min > sup) std::swap(min, sup);
    int width = sup - min;
    return min + std::rand() % width;
};

v3d randv3d(double max_r) {
    v3d r;
    for (int i = 0; i < 3; i++) {
        r(i) = rand(max_r, -max_r);
    }
    return r;
};

} // namespace rand