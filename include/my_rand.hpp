#ifndef MY_RAND_HPP
#define MY_RAND_HPP

#include <eigen3/Eigen/Dense>

using v3d = Eigen::Vector3d;

namespace my_rand {

/**
 * @brief [min, max]なる実数をランダムに返す
 * @param max 最大値
 * @param min 最小値
*/
double rand(double max=1.0, double min=1.0);

/**
 * @brief [min, max)なる整数をランダムに返す
 * @param sup 上界値
 * @param min 最小値
*/
int randint(int sup=8, int min=0);

/**
 * @brief 各成分の絶対値がmax_r以下なる3次元ベクトルをランダムに返す
 * @param max_r 各成分の絶対値の最大値
*/
v3d randv3d(double max_r);

} // namespace rand

#endif // MY_RAND_HPP