#include "NVT_MC_core.hpp"

#include "my_rand.hpp"

using ui = unsigned int;
using v3d = Eigen::Vector3d;
using v3xd = Eigen::Matrix3Xd;

namespace NVT_MC {

NVTState::NVTState(ui N, double T, v3xd r) : N(N), T(T) {
    this->r = r;
};

void NVTState::move(ui k, v3d dr) {
    double dU = 0;
    for (ui i = 0; i < N; i++) {
        // i番目とk番目の粒子の寄与によるエネルギーを更新する
        if (i == k || this->NeighborsList[k].count(i) == 0) continue;
        v3d r_ki = this->r.col(i) - this->r.col(k);
        dU = dU - this->U_ij(r_ki) + this->U_ij(r_ki + dr);
    }

    if (dU < 0
        || my_rand::rand() < exp(-dU / (physical_constants::k * T))
    ) {
        this->U += dU;
        this->r.col(k) += dr;
    }
};

void NVTState::makeNeighborList(ui P, double dr_max) {
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

void NVTState::initU() {
    this->U = 0;
    for (ui i = 0; i < N; i++) {
        for (ui j = i+1; j < N; j++) {
            if (this->NeighborsList[i].count(j) == 0) continue;
            v3d r_ij = this->r.col(i) - this->r.col(j);
            this->U += this->U_ij(r_ij);
        }
    }
}

double NVTState::U_ij(v3d r) {
    double r_norm = r.norm();
    if (r_norm > R_COFF) return 0;
    double r6 = std::pow(r_norm/SIGMA, 6);
    return 4 * EPSILON * (std::pow(r6, -2) - std::pow(r6, -1));
}

NVT_MC_Simulator::NVT_MC_Simulator(ui M, double n, double T) : M(M), N_DENSE(n), T(T) {
    ui N = M*M*M;
    this->V = N / N_DENSE;
    this->L = std::cbrt(V);
    this->a = this->L / M;
    this->MAX_dR = this->a / 64;
    this->init();
}

void NVT_MC_Simulator::run(ui n_steps) {
    for (ui i = 0; i < n_steps; i++) {
        if (i % this->REFRESH_NEIGHBOR_FREQ == 0) {
            this->r_now.makeNeighborList(this->REFRESH_NEIGHBOR_FREQ, MAX_dR);
        }
        this->shift();
    }
};

void NVT_MC_Simulator::init() {
    v3xd r = v3xd::Zero(3, N);
    for (ui i = 0; i < M; i++) {
        for (ui j = 0; j < M; j++) {
            for (ui k = 0; k < M; k++) {
                r.col(i*M*M + j*M + k) = v3d(i, j, k) * a;
            }
        }
    }
    this->r_now = NVTState(N, T, r);
    this->r_now.makeNeighborList(REFRESH_NEIGHBOR_FREQ, MAX_dR);
    this->r_now.initU();
    this->U.push_back(r_now.U);
};

void NVT_MC_Simulator::shift() {
    int k = my_rand::randint(0, N);
    Eigen::Vector3d dr = my_rand::randv3d(this->MAX_dR);
    this->r_now.move(k, dr);
    this->U.push_back(r_now.U);
};

}; // namespace nvt_mc
