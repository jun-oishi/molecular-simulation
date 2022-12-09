#ifndef MY_RAND
#define MY_RAND

#include <bits/stdc++.h>
#include <eigen3/Eigen/Dense>

using v3d = Eigen::Vector3d;

namespace rand {

/**
 * @brief [min, max]なる実数をランダムに返す
 * @param max 最大値
 * @param min 最小値
*/
double rand(double max=1.0, double min=1.0) {
    if (min > max) std::swap(min, max);
    double width = max - min;
    return (double) min + max * std::rand() / RAND_MAX;
};

/**
 * @brief [min, max)なる整数をランダムに返す
 * @param sup 上界値
 * @param min 最小値
*/
int randint(int sup=8, int min=0) {
    if (min > sup) std::swap(min, sup);
    int width = sup - min;
    return min + std::rand() % width;
};

/**
 * @brief 各成分の絶対値がmax_r以下なる3次元ベクトルをランダムに返す
 * @param max_r 各成分の絶対値の最大値
*/
v3d randv3d(double max_r) {
    v3d r;
    for (int i = 0; i < 3; i++) {
        r(i) = rand(max_r, -max_r);
    }
    return r;
};

} // namespace rand

#endif // MY_RAND