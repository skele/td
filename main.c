/* Compute time-variable spectrum from the tidal disruption echo in a
circum-binary gas disc as seen from a specific observation angle
Patrick Brem, AEI Potsdam, 07/2012 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define VSCALE 613476.295
#define TIMEINDAYS 23286.24 

void loopthroughgas(double* bhpos, double clight, double phi, double theta)
{
  int id;
  double mass,rij,baselambda,lambda,time,rij_eye;
  FILE *f, *fout;
  char line[160];
  double tpos[3],tv[3];
  double rotate[3],v,strength;
  baselambda = 656.0;
  f = fopen("gas.dat", "rt");
  fout = fopen("out.dat", "w");
  //define rotation matrix; need only to obtain the final z component so only need last line of the matrix, e.g. a vector;
  theta = theta/180.0*3.14159;
  phi = phi/180.0*3.14159;
  rotate[0] = sin(theta)*sin(phi);
  rotate[1] = sin(theta)*cos(phi);
  rotate[2] = cos(theta);
  while (fgets(line, 160, f) != NULL)
    {
      sscanf(line, "%d %lg %lg %lg %lg %lg %lg %lg\n",&id,&mass,&tpos[0],&tpos[1],&tpos[2],&tv[0],&tv[1],&tv[2]);
      rij = sqrt((tpos[0]-bhpos[0])*(tpos[0]-bhpos[0])+(tpos[1]-bhpos[1])*(tpos[1]-bhpos[1])+(tpos[2]-bhpos[2])*(tpos[2]-bhpos[2]));
      //need also distance to the observer, who is located at the direction (phi,theta) in degrees
      //rotation so that z points into observer direction (confirm definitions/notations.....)
      rij_eye = rotate[0]*tpos[0]+rotate[1]*tpos[1]+rotate[2]*tpos[2];
      //also project velocity to the line of sight
      v = rotate[0]*tv[0]+rotate[1]*tv[1]+rotate[2]*tv[2];
      time = (rij+rij_eye)/clight*TIMEINDAYS; //when the TD reaches the gas particle and then the radiation reaches the eye
      lambda = baselambda/(1.0-v/clight);
      strength = 1.0/(rij*rij);
      //strength is because of the weaker TD impact over distance^2
      fprintf(fout,"%f\t%f\t%f\n",time,lambda,strength);
    }

}


int main (int argc, char** argv)
{

  double bhpos[3], clight, phi, theta;
  bhpos[0] = 0.012866;
  bhpos[1] = 0.3103;
  bhpos[2] = -0.02399;
  clight = 3.0e8/VSCALE;
  phi = 0.0;
  theta = 45.0;

  loopthroughgas(bhpos,clight,phi,theta);


  return 0;

}
