#ifndef NVT_MC_CORE_HPP
#define NVT_MC_CORE_HPP

#include <vector>
#include <set>
#include <eigen3/Eigen/Dense>
#include "physical_constants.hpp"

using ui = unsigned int;
using v3d = Eigen::Vector3d;
using v3xd = Eigen::Matrix3Xd;

namespace NVT_MC {

class NVTState {
  public:
    v3xd r; // 各粒子の座標
    const ui N; // 粒子数
    double U; // 全系のポテンシャル
    const double T; // 温度
    std::vector<std::set<ui>> NeighborsList; // 隣接リスト
    // レナード・ジョーンズ相互作用の定数(Arの場合)
    const double EPSILON = physical_constants::k * 119.8; // Arのepsilon
    const double SIGMA = 3.405e-10; // Arのsigma
    double R_COFF = 2 * SIGMA; // カットオフ半径

    /**
     * @brief コンストラクタ
     * @param N 粒子数
     * @param T 温度
     * @param r 粒子の座標を表す3xN行列
    */
    NVTState(ui N, double T, v3xd r);

    /**
     * @brief 粒子を動かす試行を行い、結果を評価して受け入れる場合は座標とUを更新する
     * @param k 動かす粒子の番号
     * @param dr 移動のベクトル
    */
    void move(ui k, v3d dr);

    /**
     * @brief 隣接リストを更新する
     * @param P 隣接リストの更新頻度
     * @param dr_max 一ステップでの最大移動距離
    */
    void makeNeighborList(ui P, double dr_max);

    /**
     * @brief 全体のポテンシャルを計算する
    */
    void initU();

  private:
    /**
     * @brief レナード・ジョーンズ相互作用のポテンシャル
     * @param r 粒子間の相対位置ベクトル
    */
    double U_ij(v3d r);
};

class NVT_MC_Simulator {
  public:
    ui N; // 粒子数
    const ui M; // 1辺に並ぶ格子数
    const double N_DENSE; // 数密度
    double L; // 全体領域の一辺の長さ
    double a; // 格子点間の距離
    const ui REFRESH_NEIGHBOR_FREQ = 20; // 隣接リストの更新頻度
    double MAX_dR; // 一ステップでの最大移動距離
    const double T; // 温度
    double V; // 全体領域の体積

    std::vector<double> U; // エネルギーの履歴

    NVT_MC_Simulator(ui M, double n, double T);

    // 単純格子点上に初期配置して開始
    void run(ui n_steps=0xff);

  private:
    NVTState r_now; // 現在の状態
    // 単純格子点上に粒子を配置する
    void init();

    void shift();

};

}; // namespace nvt_mc

#endif // NVT_MC_HPP