/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * jello, struct point t[8][8][8]);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);


// Implemented Functions for calculations:
double calcDist(struct point t, struct point b);
point applyHooke(double k, double rest, struct point t, struct point b);
point applyDampen(double d, struct point t, struct point b, struct point t1_vec, struct point t2_vec);

#endif

