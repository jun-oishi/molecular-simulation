#include <stdio.h>
#include <stdlib.h>
#include <time.h>

const double LATTICE_SIZE = 10.0;
double max_velocity = LATTICE_SIZE / 20;

double frand(double min, double max) {
  return min + (max - min) * (double)rand() / RAND_MAX;
  return (double)rand() / RAND_MAX;
}

void move1d(double *x, double *v) {
  *x += *v;
  if (*x < 0) {
    *x = -*x;
    *v = -*v;
  } else if (*x > LATTICE_SIZE) {
    *x = 2 * LATTICE_SIZE - *x;
    *v = -*v;
  }
}

void move(double *x, double *y, double *z, double *u, double *v, double *w) {
  move1d(x, u);
  move1d(y, v);
  move1d(z, w);
}

int main() {
  const int NUM_PARTICLES = 30;
  const int NUM_STEPS = 1000;

  double x[NUM_PARTICLES], y[NUM_PARTICLES], z[NUM_PARTICLES];
  double u[NUM_PARTICLES], v[NUM_PARTICLES], w[NUM_PARTICLES];

  srand(time(NULL));
  for (int i=0; i<NUM_PARTICLES; i++) {
    x[i] = frand(0, LATTICE_SIZE);
    y[i] = frand(0, LATTICE_SIZE);
    z[i] = frand(0, LATTICE_SIZE);
    u[i] = frand(-max_velocity, max_velocity);
    v[i] = frand(-max_velocity, max_velocity);
    w[i] = frand(-max_velocity, max_velocity);
  }

  for (int i=0; i<NUM_STEPS; i++) {
    printf("%d\n", NUM_PARTICLES);
    printf("Lattice=\"%lf 0 0 0 %lf 0 0 0 %lf\"", LATTICE_SIZE, LATTICE_SIZE, LATTICE_SIZE);
    printf(" Properties=id:I:1:element:S:1:radius:R:1:pos:R:3");
    printf(" Time=%lf", (double)i);
    printf("\n");
    for (int j=0; j<NUM_PARTICLES; j++) {
      char element = 'A' + (j % 3);
      double radius = 0.9 + 0.1 * (j % 3);
      printf("%d %c %lf %lf %lf %lf\n", j, element, radius, x[j], y[j], z[j]);
      move(&x[j], &y[j], &z[j], &u[j], &v[j], &w[j]);
    }
  }
  return 0;
}