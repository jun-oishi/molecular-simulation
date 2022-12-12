#include <bits/stdc++.h>
#include "NVT_MC_core.hpp"

int main() {
    int M, n_steps;
    double n, T;
    std::cout << "M: ";
    std::cin >> M;
    std::cout << "n: ";
    std::cin >> n;
    std::cout << "T: ";
    std::cin >> T;
    std::cout << "n_steps: ";
    std::cin >> n_steps;

    NVT_MC::NVT_MC_Simulator sim(M, n, T);
    sim.run(n_steps);

    std::string filename;
    std::cout << "filename: ";
    std::cin >> filename;
    std::ofstream ofs(filename);
    ofs << "M:" << M << ",n:" << n << ",T:" << T << ",n_steps:" << n_steps << std::endl;
    for (int i = 0; i < n_steps; i++) {
        ofs << i << "," << sim.U[i] << std::endl;
    }
    ofs.close();

    return 0;
}