
#include "jello.h"
#include "physics.h"
#include <string>

// Calculate Distance.
double calcDist(struct point a, struct point b)
{
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

// Apply Hooke's Law.
point applyHooke(double k, double rest_len, struct point a, struct point b)
{
	point new_d;
	double dist, T;

	dist = calcDist(a, b);
	T = -k * (dist - rest_len);

	new_d.x = T * (a.x - b.x) / dist;
	new_d.y = T * (a.y - b.y) / dist;
	new_d.z = T * (a.z - b.z) / dist;

	return new_d;
}

// Apply Dampening.
point applyDampen(double kd, struct point a, struct point b, struct point t1_vec, struct point t2_vec)
{
	point new_d, L;
	double c, dist, T;

	pDIFFERENCE(t1_vec, t2_vec, L);

	dist = calcDist(a, b);

	c = (L.x * (a.x - b.x) + L.y * (a.y - b.y) + L.z * (a.z - b.z)) / dist;

	T = -(kd * c);

	new_d.x = T * (a.x - b.x) / dist;
	new_d.y = T * (a.y - b.y) / dist;
	new_d.z = T * (a.z - b.z) / dist;

	return new_d;
}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 't'. */
void computeAcceleration(struct world * jello, struct point t[8][8][8])
{
	int a, b, c; // (i,j,k)
	double rest_length = 0.142857, rest_shear, rest_bend, rest_diag;

	for (a = 0; a <= 7; a++)
		for (b = 0; b <= 7; b++)
			for (c = 0; c <= 7; c++)
			{
				point N;	// normal
				point collide_F = {}, F = {};	// (force)

				rest_shear = rest_length * sqrt(2);
				rest_diag = rest_length * sqrt(3);
				rest_bend = rest_length * 2;

				// Bend Spring Setup:

				if (a + 2 > -1 && a + 2 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a + 2][b][c], jello->v[a][b][c], jello->v[a + 2][b][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_bend, jello->p[a][b][c], jello->p[a + 2][b][c]), F);
				}

				if (a - 2 > -1 && a - 2 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a - 1][b][c], jello->v[a][b][c], jello->v[a - 1][b][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_bend, jello->p[a][b][c], jello->p[a - 2][b][c]), F);
				}

				if (b + 2 > -1 && b + 2 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b + 2][c], jello->v[a][b][c], jello->v[a][b + 2][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_bend, jello->p[a][b][c], jello->p[a][b + 2][c]), F);
				}

				if (b - 2 > -1 && b - 2 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b - 2][c], jello->v[a][b][c], jello->v[a][b - 2][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_bend, jello->p[a][b][c], jello->p[a][b - 2][c]), F);
				}

				if (c + 2 > -1 && c + 2 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b][c + 2], jello->v[a][b][c], jello->v[a][b][c + 2]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_bend, jello->p[a][b][c], jello->p[a][b][c + 2]), F);
				}

				if (c - 2 > -1 && c - 2 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b][c - 2], jello->v[a][b][c], jello->v[a][b][c - 2]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_bend, jello->p[a][b][c], jello->p[a][b][c - 2]), F);
				}

				// Structural Spring Setup:

				if (a + 1 > -1 && a + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a + 1][b][c], jello->v[a][b][c], jello->v[a + 1][b][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_length, jello->p[a][b][c], jello->p[a + 1][b][c]), F);
				}

				if (a - 1 > -1 && a - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a - 1][b][c], jello->v[a][b][c], jello->v[a - 1][b][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_length, jello->p[a][b][c], jello->p[a - 1][b][c]), F);
				}

				if (b + 1 > -1 && b + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b + 1][c], jello->v[a][b][c], jello->v[a][b + 1][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_length, jello->p[a][b][c], jello->p[a][b + 1][c]), F);
				}

				if (b - 1 > -1 && b - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b - 1][c], jello->v[a][b][c], jello->v[a][b - 1][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_length, jello->p[a][b][c], jello->p[a][b - 1][c]), F);
				}

				if (c + 1 > -1 && c + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b][c + 1], jello->v[a][b][c], jello->v[a][b][c + 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_length, jello->p[a][b][c], jello->p[a][b][c + 1]), F);
				}

				if (c - 1 > -1 && c - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b][c - 1], jello->v[a][b][c], jello->v[a][b][c - 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_length, jello->p[a][b][c], jello->p[a][b][c - 1]), F);
				}

				// Shear Spring Setup:

				if (a + 1 > -1 && a + 1 < 8 && b + 1 > -1 && b + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a + 1][b + 1][c], jello->v[a][b][c], jello->v[a + 1][b + 1][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a + 1][b + 1][c]), F);
				}

				if (a - 1 > -1 && a - 1 < 8 && b + 1 > -1 && b + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a - 1][b + 1][c], jello->v[a][b][c], jello->v[a - 1][b + 1][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a - 1][b + 1][c]), F);
				}

				if (a - 1 > -1 && a - 1 < 8 && b - 1 > -1 && b - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a - 1][b - 1][c], jello->v[a][b][c], jello->v[a - 1][b - 1][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a - 1][b - 1][c]), F);
				}

				if (a + 1 > -1 && a + 1 < 8 && b - 1 > -1 && b - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a + 1][b - 1][c], jello->v[a][b][c], jello->v[a + 1][b - 1][c]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a + 1][b - 1][c]), F);
				}

				if (b + 1 > -1 && b + 1 < 8 && c + 1 > -1 && c + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b + 1][c + 1], jello->v[a][b][c], jello->v[a][b + 1][c + 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a][b + 1][c + 1]), F);
				}

				if (b - 1 > -1 && b - 1 < 8 && c + 1 > -1 && c + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b - 1][c + 1], jello->v[a][b][c], jello->v[a][b - 1][c + 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a][b - 1][c + 1]), F);
				}

				if (b - 1 > -1 && b - 1 < 8 && c - 1 > -1 && c - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b - 1][c - 1], jello->v[a][b][c], jello->v[a][b - 1][c - 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a][b - 1][c - 1]), F);
				}

				if (b + 1 > -1 && b + 1 < 8 && c - 1 > -1 && c - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a][b + 1][c - 1], jello->v[a][b][c], jello->v[a][b + 1][c - 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a][b + 1][c - 1]), F);
				}

				if (a + 1 > -1 && a + 1 < 8 && c + 1 > -1 && c + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a + 1][b][c + 1], jello->v[a][b][c], jello->v[a + 1][b][c + 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a + 1][b][c + 1]), F);
				}

				if (a - 1 > -1 && a - 1 < 8 && c + 1 > -1 && c + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a - 1][b][c + 1], jello->v[a][b][c], jello->v[a - 1][b][c + 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a - 1][b][c + 1]), F);
				}

				if (a - 1 > -1 && a - 1 < 8 && c - 1 > -1 && c - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a - 1][b][c - 1], jello->v[a][b][c], jello->v[a - 1][b][c - 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a - 1][b][c - 1]), F);
				}

				if (a + 1 > -1 && a + 1 < 8 && c - 1 > -1 && c - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a + 1][b][c - 1], jello->v[a][b][c], jello->v[a + 1][b][c - 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_shear, jello->p[a][b][c], jello->p[a + 1][b][c - 1]), F);
				}

				// For Diagonals.

				if (a + 1 > -1 && a + 1 < 8 && b + 1 > -1 && b + 1 < 8 && c + 1 > -1 && c + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a + 1][b + 1][c + 1], jello->v[a][b][c], jello->v[a + 1][b + 1][c + 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_diag, jello->p[a][b][c], jello->p[a + 1][b + 1][c + 1]), F);
				}

				if (a - 1 > -1 && a - 1 < 8 && b + 1 > -1 && b + 1 < 8 && c + 1 > -1 && c + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a - 1][b + 1][c + 1], jello->v[a][b][c], jello->v[a - 1][b + 1][c + 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_diag, jello->p[a][b][c], jello->p[a - 1][b + 1][c + 1]), F);
				}

				if (a - 1 > -1 && a - 1 < 8 && b - 1 > -1 && b - 1 < 8 && c + 1 > -1 && c + 1 < 8) 
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a - 1][b - 1][c + 1], jello->v[a][b][c], jello->v[a - 1][b - 1][c + 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_diag, jello->p[a][b][c], jello->p[a - 1][b - 1][c + 1]), F);
				}

				if (a + 1 > -1 && a + 1 < 8 && b - 1 > -1 && b - 1 < 8 && c + 1 > -1 && c + 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a + 1][b - 1][c + 1], jello->v[a][b][c], jello->v[a + 1][b - 1][c + 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_diag, jello->p[a][b][c], jello->p[a + 1][b - 1][c + 1]), F);
				}

				if (a + 1 > -1 && a + 1 < 8 && b + 1 > -1 && b + 1 < 8 && c - 1 > -1 && c - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a + 1][b + 1][c - 1], jello->v[a][b][c], jello->v[a + 1][b + 1][c - 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_diag, jello->p[a][b][c], jello->p[a + 1][b + 1][c - 1]), F);
				}

				if (a - 1 > -1 && a - 1 < 8 && b + 1 > -1 && b + 1 < 8 && c - 1 > -1 && c - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a - 1][b + 1][c - 1], jello->v[a][b][c], jello->v[a - 1][b + 1][c - 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_diag, jello->p[a][b][c], jello->p[a - 1][b + 1][c - 1]), F);
				}

				if (a - 1 > -1 && a - 1 < 8 && b - 1 > -1 && b - 1 < 8 && c - 1 > -1 && c - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a - 1][b - 1][c - 1], jello->v[a][b][c], jello->v[a - 1][b - 1][c - 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_diag, jello->p[a][b][c], jello->p[a - 1][b - 1][c - 1]), F);
				}

				if (a + 1 > -1 && a + 1 < 8 && b - 1 > -1 && b - 1 < 8 && c - 1 > -1 && c - 1 < 8)
				{
					pSUM(F, applyDampen(jello->dElastic, jello->p[a][b][c], jello->p[a + 1][b - 1][c - 1], jello->v[a][b][c], jello->v[a + 1][b - 1][c - 1]), F);
					pSUM(F, applyHooke(jello->kElastic, rest_diag, jello->p[a][b][c], jello->p[a + 1][b - 1][c - 1]), F);
				}

				// Calculate for Collisions:

				if (jello->incPlanePresent == 1)		// If an Inclined Plane is present.
				{
					point dir;
					double T1, above, below;

					// Surfaces:

					below = jello->a * jello->a + jello->b * jello->b + jello->c * jello->c;
					above = jello->a * jello->p[a][b][c].x + jello->b * jello->p[a][b][c].y + jello->c * jello->p[a][b][c].z + jello->d;

					T1 = -(above / below);

					dir.x = jello->p[a][b][c].x - (jello->a * T1 + jello->p[a][b][c].x);
					dir.y = jello->p[a][b][c].y - (jello->b * T1 + jello->p[a][b][c].y);
					dir.z = jello->p[a][b][c].z - (jello->c * T1 + jello->p[a][b][c].z);

					// Find colliding force:

					if (dir.x * Normal.x + dir.y * Normal.y + dir.z * Normal.z < 0)
					{
						collide_F.x += jello->kCollision * fabs(dir.x) * Normal.x - 1.0 * jello->dCollision * jello->v[a][b][c].x;
						collide_F.y += jello->kCollision * fabs(dir.y) * Normal.y - 1.0 * jello->dCollision * jello->v[a][b][c].y;
						collide_F.z += jello->kCollision * fabs(dir.z) * Normal.z - 1.0 * jello->dCollision * jello->v[a][b][c].z;
					}
				}

				if (jello->p[a][b][c].x > 2)
				{
					N.x = -1;
					collide_F.x += -jello->dCollision * jello->v[a][b][c].x + jello->kCollision * fabs(jello->p[a][b][c].x - 2) * N.x;
					N.z = N.y = 0;
				}

				if (jello->p[a][b][c].x < -2)
				{
					N.x = 1;
					collide_F.x += -jello->dCollision * jello->v[a][b][c].x + jello->kCollision * fabs(jello->p[a][b][c].x + 2) * N.x;
					N.z = N.y = 0;
				}

				if (jello->p[a][b][c].y > 2)
				{
					N.y = -1;
					collide_F.y += -jello->dCollision * jello->v[a][b][c].y + jello->kCollision * fabs(jello->p[a][b][c].y - 2) * N.y;
					N.z = N.x = 0;
				}

				if (jello->p[a][b][c].y < -2)
				{
					N.y = 1;
					collide_F.y += -jello->dCollision * jello->v[a][b][c].y + jello->kCollision * fabs(jello->p[a][b][c].y + 2) * N.y;
					N.z = N.x = 0;
				}

				if (jello->p[a][b][c].z > 2)
				{
					N.z = -1;
					collide_F.z += -jello->dCollision * jello->v[a][b][c].z + jello->kCollision * fabs(jello->p[a][b][c].z - 2) * N.z;
					N.y = N.x = 0;
				}

				if (jello->p[a][b][c].z < -2)
				{
					N.z = 1;
					collide_F.z += -jello->dCollision * jello->v[a][b][c].z + jello->kCollision * fabs(jello->p[a][b][c].z + 2) * N.z;
					N.y = N.x = 0;
				}

				pSUM(F, collide_F, F);		// Combine.

				// Force Field Calculation:

				double size;
				int x=0, y=0, z=0, temp, new_x, new_y, new_z;

				if (jello->resolution != 0)
				{

					size = 4.0 / jello->resolution;

					if (jello->p[a][b][c].x >= 2)
						x = (4 / size) - 1;
					else
						x = (jello->p[a][b][c].x + 2) / size;

					if (jello->p[a][b][c].y >= 2)
						y = (4 / size) - 1;
					else
						y = (jello->p[a][b][c].y + 2) / size;

					if (jello->p[a][b][c].z >= 2)
						z = (4 / size) - 1;
					else
						z = (jello->p[a][b][c].z + 2) / size;

					temp = x * jello->resolution * jello->resolution + y * jello->resolution + z;

					pSUM(F, jello->forceField[temp], F);

				}

				new_x = F.x * 100;
				F.x = new_x / 100.0;

				new_y = F.y * 100;
				F.y = new_y / 100.0;

				new_z = F.z * 100;
				F.z = new_z / 100.0;

				pMULTIPLY(F, 1/jello->mass, t[a][b][c]);
			}
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point t[8][8][8];

  computeAcceleration(jello, t);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * t[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * t[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * t[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point t[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, t);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(t[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, t);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(t[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, t);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(t[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, t);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(t[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
