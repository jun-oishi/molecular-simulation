#include <bits/stdc++.h>
#include <eigen3/Eigen/Dense>
#include "myrand.cpp"
#include "physical_constants.cpp"

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
    NVTState(ui N, double T, v3xd r) : N(N), T(T) {
        this->r = r;
    };

    /**
     * @brief 粒子を動かす試行
     * @param k 動かす粒子の番号
     * @param dr 移動のベクトル
    */
    void move(ui k, v3d dr) {
        double dU = 0;
        for (ui i = 0; i < N; i++) {
            // i番目とk番目の粒子の寄与によるエネルギーを更新する
            if (i == k || this->NeighborsList[k].count(i) == 0) continue;
            v3d r_ki = this->r.col(i) - this->r.col(k);
            dU = dU - this->U_ij(r_ki) + this->U_ij(r_ki + dr);
        }

        if (dU < 0
            || rand::rand() < exp(-dU / (physical_constants::k * T))
        ) {
            this->U += dU;
            this->r.col(k) += dr;
        }
    };

    /**
     * @brief 隣接リストを更新する
     * @param P 隣接リストの更新頻度
     * @param dr_max 一ステップでの最大移動距離
    */
    void makeNeighborList(ui P, double dr_max) {
        double r = R_COFF + 2*P*dr_max;
        for (int i = 0; i < this->N; i++) {
            for (int j = i+1; j < this->N; j++) {
                v3d r_ij = this->r.col(i) - this->r.col(j);
                if (r_ij.norm() < r) {
                    this->NeighborsList[i].insert(j);
                    this->NeighborsList[j].insert(i);
                }
            }
        }
    };

    /**
     * @brief 全体のポテンシャルを計算する
    */
    void initU() {
        this->U = 0;
        for (ui i = 0; i < N; i++) {
            for (ui j = i+1; j < N; j++) {
                if (this->NeighborsList[i].count(j) == 0) continue;
                v3d r_ij = this->r.col(i) - this->r.col(j);
                this->U += this->U_ij(r_ij);
            }
        }
    }

  private:
    /**
     * @brief レナード・ジョーンズ相互作用のポテンシャル
     * @param r 粒子間の相対位置ベクトル
    */
    double U_ij(v3d r) {
        double r_norm = r.norm();
        if (r_norm > R_COFF) return 0;
        double r6 = std::pow(r_norm/SIGMA, 6);
        return 4 * EPSILON * (std::pow(r6, -2) - std::pow(r6, -1));
    };
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

    NVT_MC_Simulator(ui M, double n, double T) : M(M), N_DENSE(n), T(T) {
        ui N = M*M*M;
        this->V = N / N_DENSE;
        this->L = std::cbrt(V);
        this->a = this->L / M;
        this->MAX_dR = this->a / 64;
    }

    // 単純格子点上に初期配置して開始
    void run(ui n_steps=0xff) {
        this->init();

        for (ui i = 0; i < n_steps; i++) {
            if (i % REFRESH_NEIGHBOR_FREQ == 0) {
                this->r_now.makeNeighborList(REFRESH_NEIGHBOR_FREQ, MAX_dR);
            }
            this->shift();
        }
    }

  private:
    NVTState r_now;
    // 単純格子点上に粒子を配置する
    NVTState init() {
        v3xd r = v3xd::Zero(3, N);
        for (ui i = 0; i < M; i++) {
            for (ui j = 0; j < M; j++) {
                for (ui k = 0; k < M; k++) {
                    r.col(i*M*M + j*M + k) = v3d(i, j, k) * a;
                }
            }
        }
        NVTState r_now = NVTState(N, T, r);
        r_now.makeNeighborList(REFRESH_NEIGHBOR_FREQ, MAX_dR);
        r_now.initU();
        this->U.push_back(r_now.U);
    }

    NVTState shift() {
        int k = rand::randint(0, N);
        Eigen::Vector3d dr = rand::randv3d(this->MAX_dR);
        NVTState r_new = r_now;
        r_new.move(k, dr);
        this->U.push_back(r_now.U);
    }

};

}; // namespace nvt_mc
