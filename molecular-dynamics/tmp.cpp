#include <cmath>
#include <iostream>
#include <random>
#include <string>

using namespace std;

const uint32_t SEED = 6;

double drand(double min, double max) {
  static minstd_rand0 engine(SEED);
  static const uint32_t rand_max = engine.max();
  return min + (max-min)*((double)engine()/rand_max);
}

double maxwell(double T, double m) {
  double r1 = drand(0, 1);
  double r2 = drand(0, 1);
  double ret =
      sqrt(-2 * T * log(r1) / m) * cos(2 * M_PI * r2);
  if (isnan(ret)) {
    cout << "nan occurred" << endl;
    cout << "r1: " << r1 << endl;
    cout << "r2: " << r2 << endl;
    cout << "T: " << T << endl;
    cout << "m: " << m << endl;
    exit(1);
  }
  return ret;
}

int main() {
  for (int i = 0; i < 300; i+=10) {
    cout << "log(1e-"<<i<<") = " << log(pow(10, -i)) << endl;
  }
  cout << "no nan occurred" << endl;
}