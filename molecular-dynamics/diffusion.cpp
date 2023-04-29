
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

const int NA = 20; // Aの粒子数
const int NB = 20; // Bの粒子数
const double k = 10; // Bの質量/Aの質量
const double T = 5; // 無次元温度
const double dt = 0.001; // 時間ステップ幅
const double n_dense = 0.1; // 粒子の数密度
const double L = cbrt((NA+NB)/n_dense);  // シミュレーション領域の1辺

const int n_steps = 10000;
const int save_interval = 10;

struct vec3d {
  double x, y, z;

  vec3d() {};

  vec3d(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  vec3d(const vec3d &other) {
    x = other.x;
    y = other.y;
    z = other.z;
  }

  vec3d& operator=(const vec3d &other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }

  vec3d& operator+=(const vec3d &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  vec3d& operator-=(const vec3d &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  vec3d& operator*=(double a) {
    x *= a;
    y *= a;
    z *= a;
    return *this;
  }

  vec3d& operator/=(double a) {
    x /= a;
    y /= a;
    z /= a;
    return *this;
  }

  double norm() {
    return sqrt(x*x + y*y + z*z);
  }
};

vec3d operator+(const vec3d &a, const vec3d &b) {return vec3d(a) += b;}
vec3d operator-(const vec3d &a, const vec3d &b) {return vec3d(a) -= b;}
vec3d operator*(const vec3d &a, double b) {return vec3d(a) *= b;}
vec3d operator*(double a, const vec3d &b) {return vec3d(b) *= a;}
vec3d operator/(const vec3d &a, double b) {return vec3d(a) /= b;}

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
  srand(time(NULL));
  cout << "RAND_MAX = " << RAND_MAX << endl;

  string outfile = "diffusion1.xyz";

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

  double t = 0;
  for (int i=0; i<n_steps; i++) {
    initialize(particles);
    compute_force(particles);
    move(particles);
    if (i % save_interval == 0) {
      write(particles, t, outfile);
    }
    t += dt;
  }
}

double frand(double min, double max) {
  double ret = min + (max-min)*rand()/RAND_MAX;
  if (isnan(ret)) {
    ret = frand(min, max);
  }
  cout << ret << " ";
  return ret;
}

// T, m は無次元化された値
double maxwell(double T, double m) {
  return sqrt(-2*T*log(frand(0,1))/m) * cos(2*M_PI*frand(0,1));
}

void initialize(Particle *particles) {
  // 初期位置の設定
  for (int i=0; i<NA+NB; i++) {
    particles[i].r.x = frand(0, L);
    particles[i].r.y = frand(0, L);
    particles[i].r.z = frand(0, L);

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
  vec3d avg_mom_A = total_mom_A / NA;
  for (int i=0; i<NA; i++) {
    particles[i].v -= avg_mom_A / 1;
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
  vec3d avg_mom_B = total_mom_B / NA;
  for (int i = 0; i < NA; i++) {
    particles[i].v -= avg_mom_B / k;
  }

}

void compute_force(Particle *particles) {
  // TODO 計算対象をVerlet neighborで絞りたい
  // TODO 周期的境界条件...
  for (int i=0; i<NA+NB; i++) {
    for (int j=i+1; j<NA+NB; j++) {
      vec3d fij = particles[i].lj_force(particles[j]);
      particles[i].f = fij;
      vec3d zero = {0, 0, 0};
      vec3d fji = zero - fij;
      particles[j].f = fji;
    }
  }
}

void move(Particle *particles) {
  for (int i=0; i<NA+NB; i++) {
    double m = (particles[i].type == 'A') ? 1 : k;
    particles[i].r = particles[i].r + particles[i].v * dt + particles[i].f * (dt*dt / (2*m));
    particles[i].v += particles[i].f * (dt / m);

    if (particles[i].r.x < 0) {
      particles[i].r.x += L;
    } else if (particles[i].r.x > L) {
      particles[i].r.x -= L;
    }
    if (particles[i].r.y < 0) {
      particles[i].r.y += L;
    } else if (particles[i].r.y > L) {
      particles[i].r.y -= L;
    }
    if (particles[i].r.z < 0) {
      particles[i].r.z += L;
    } else if (particles[i].r.z > L) {
      particles[i].r.z -= L;
    }
  }
}

void write(Particle *particles, double time, string outfile) {
  ofstream ofs(outfile, ios::app);
  ofs << NA+NB << endl;
  ofs << "Lattice \"" << L << " 0 0 0 " << L << " 0 0 0 " << L << "\"";
  ofs << " Properties=id:I:1:element:S:1:radius:R:1:pos:R:3";
  ofs << " Time=" << time;
  ofs << endl;

  for (int i = 0; i < NA + NB; i++) {
    double r = (particles[i].type == 'A') ? 0.2 : 0.02*cbrt(k);
    ofs << i << " ";
    ofs << particles[i].type << " ";
    ofs << r << " ";
    ofs << particles[i].r.x << " ";
    ofs << particles[i].r.y << " ";
    ofs << particles[i].r.z << endl;
  }
  ofs.close();
}