
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <random>

#include "vec3d.h"

using namespace std;

const int NA = 30; // Aの粒子数
const int NB = 30; // Bの粒子数
const double k = 4; // Bの質量/Aの質量
const double T = 5; // 無次元温度
const double dt = 0.001; // 時間ステップ幅
const double n_dense = 0.3; // 粒子の数密度
const double L = cbrt((NA+NB)/n_dense);  // シミュレーション領域の1辺

const int n_steps = 10000;
const int save_interval = 10;
int current_step;

const string outdir = "out/";

// const uint32_t SEED = 10;
const uint32_t SEED = time(NULL);

struct Particle {
  char type;
  vec3d r, v, f;

  // レナードジョーンズポテンシャルに従う力(otherから受ける力)を計算する
  vec3d lj_force(Particle &other) {
    double dr = (this->r - other.r).norm();
    double r2 = dr*dr;
    double r6 = r2*r2*r2;
    double r12 = r6*r6;
    return (this->r - other.r) * 24 * (2 / r12 - 1 / r6) / r2;
  }
};

// min以上max以下のdoubleを返す
double frand(double min, double max);
// (3次元)maxwell分布に従う速度を返す
double maxwell(double T, double m);

// 粒子の初期位置と初期速度を設定する
void initialize(Particle *particles);
// 粒子に働く力を計算する
void compute_force(Particle *particles);
// 粒子を次の座標に移動させる
void move(Particle *particles);
// 粒子の位置を書き出す
void write(Particle *particles, double time, string outfile);

int main() {
  string outfile = "diffusion.xyz";
  // cout << "outfile: ";
  // cin >> outfile;
  outfile = outdir + outfile;

  Particle particles[NA+NB];
  for (int i=0; i<NA; i++) {
    particles[i].type = 'A';
  }
  for (int i=NA; i<NA+NB; i++) {
    particles[i].type = 'B';
  }
  initialize(particles);

  ofstream file(outfile, ios::out);
  file.close();

  write(particles, 0, outfile);

  double t = dt;
  for (current_step=1; current_step<=n_steps; current_step++) {
    compute_force(particles);
    move(particles);
    if (current_step % save_interval == 0) {
      write(particles, t, outfile);
    }
    t += dt;
  }
}

double drand(double min, double max) {
  static minstd_rand0 engine(SEED);
  static const uint32_t rand_max = engine.max();
  return min + (max - min) * ((double)engine() / rand_max);
}

// T, m は無次元化された値
double maxwell(double T, double m) {
  double ret = sqrt(-2*T*log(drand(1e-50,1))/m) * cos(2*M_PI*drand(0,1));
  if (isnan(ret)) {
    cout << "nan" << endl;
    exit(1);
  }
  return ret;
}

void initialize(Particle *particles) {
  // 初期位置の設定
  for (int i=0; i<NA+NB; i++) {
    particles[i].r.x = drand(0, L);
    particles[i].r.y = drand(0, L);
    particles[i].r.z = drand(0, L);

    // 近すぎる場合は置き直す
    bool overlap = false;
    for (int j=0; j<i; j++) {
      double dist = (particles[i].r - particles[j].r).norm();
      if (dist < 1) {
        overlap = true;
        i--;
        break;
      }
    }
    if (overlap) continue;
  }

  // 初期速度の設定
  // まずは軽いAから
  for (int i=0; i<NA; i++) {
    double m = 1;
    particles[i].v.x = maxwell(T, m);
    particles[i].v.x = maxwell(T, m);
    particles[i].v.y = maxwell(T, m);
  }
  // 合計運動量を0にする
  vec3d total_mom_A = {0, 0, 0};
  for (int i=0; i<NA; i++) {
    vec3d mom = particles[i].v * 1;
    total_mom_A += mom;
  }
  if (total_mom_A.norm() > 1e-10) {
    vec3d avg_mom_A = total_mom_A / NA;
    for (int i=0; i<NA; i++) {
      particles[i].v -= avg_mom_A / 1;
    }
  }

  // 次に重いB
  for (int i = NA; i < NA+NB; i++) {
    particles[i].v.x = maxwell(T, k);
    particles[i].v.x = maxwell(T, k);
    particles[i].v.y = maxwell(T, k);
  }
  // 合計運動量を0にする
  vec3d total_mom_B = {0, 0, 0};
  for (int i = 0; i < NA; i++) {
    vec3d mom = particles[i].v * k;
    total_mom_B += mom;
  }
  if (total_mom_B.norm() > 1e-10) {
    vec3d avg_mom_B = total_mom_B / NB;
    for (int i = NA; i < NA+NB; i++) {
      particles[i].v -= avg_mom_B / k;
    }
  }
}

void compute_force(Particle *particles) {
  // TODO 周期的境界条件
  // TODO 計算対象をVerlet neighborで絞りたい
  for (int i=0; i<NA+NB-1; i++) {
    for (int j=i+1; j<NA+NB; j++) {
      vec3d f = particles[i].lj_force(particles[j]);
      if (isnan(f.x) || isnan(f.y) || isnan(f.z)) {
        cout << "i:" << i << " j:" << j << endl;
        cout << "this : " << particles[i].r.x << " " << particles[i].r.y << " " << particles[i].r.z
             << endl;
        cout << "other: " << particles[j].r.x << " " << particles[j].r.y << " "
             << particles[j].r.z << endl;
        cout << "f: " << f.x << " " << f.y << " " << f.z << endl;
        exit(1);
      }
      particles[i].f = f;
      particles[j].f = -f;
    }
  }
}

void move(Particle *particles) {
  for (int i=0; i<NA+NB; i++) {
    double m = (particles[i].type == 'A') ? 1 : k;
    particles[i].r = particles[i].r + particles[i].v * dt + particles[i].f * (dt*dt / (2*m));
    particles[i].v += particles[i].f * (dt / m);

    if (particles[i].r.x < 0) {
      particles[i].r.x = -particles[i].r.x;
      particles[i].v.x = -particles[i].v.x;
    } else if (particles[i].r.x > L) {
      particles[i].r.x = 2*L-particles[i].r.x;
      particles[i].v.x = -particles[i].v.x;
    }
    if (particles[i].r.y < 0) {
      particles[i].r.y = -particles[i].r.y;
      particles[i].v.y = -particles[i].v.y;
    } else if (particles[i].r.y > L) {
      particles[i].r.y = 2*L-particles[i].r.y;
      particles[i].v.y = -particles[i].v.y;
    }
    if (particles[i].r.z < 0) {
      particles[i].r.z = -particles[i].r.z;
      particles[i].v.z = -particles[i].v.z;
    } else if (particles[i].r.z > L) {
      particles[i].r.z = 2*L-particles[i].r.z;
      particles[i].v.z = -particles[i].v.z;
    }
  }
}

void write(Particle *particles, double time, string outfile) {
  ofstream ofs(outfile, ios::app);
  ofs << NA+NB << endl;
  ofs << "Lattice \"" << L << " 0 0 0 " << L << " 0 0 0 " << L << "\"";
  ofs << " Properties=id:I:1:element:S:1:radius:R:1:pos:R:3:velo:R:3";
  ofs << " Time=" << time;
  ofs << endl;

  for (int i = 0; i < NA + NB; i++) {
    double r = (particles[i].type == 'A') ? 0.2 : 0.2*cbrt(k);
    ofs << i << " ";
    ofs << particles[i].type << " ";
    ofs << r << " ";
    ofs << particles[i].r.x << " ";
    ofs << particles[i].r.y << " ";
    ofs << particles[i].r.z << " ";
    ofs << particles[i].v.x << " ";
    ofs << particles[i].v.y << " ";
    ofs << particles[i].v.z << " ";
    ofs << endl;
  }
  ofs.close();
}