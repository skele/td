/* Compute time-variable spectrum from the tidal disruption echo in a
circum-binary gas disc as seen from a specific observation angle
Patrick Brem, AEI Potsdam, 07/2012 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct gasparticle {
  double x[3],v[3];
  struct gasparticle * next;
};

typedef struct gasparticle gp;

void loopthroughgas(double* bhpos, double clight)
{
  int i,id;
  double mass,rij,basefreq,freq,time;
  FILE *f, *fout;
  char line[160];
  double tpos[3],tv[3];
  basefreq = 1000.0;
  f = fopen("gas.dat", "rt");
  fout = fopen("out.dat", "w");
  while (fgets(line, 160, f) != NULL)
    {
      sscanf(line, "%d %lg %lg %lg %lg %lg %lg %lg\n",&id,&mass,&tpos[0],&tpos[1],&tpos[2],&tv[0],&tv[1],&tv[2]);
      rij = sqrt((tpos[0]-bhpos[0])*(tpos[0]-bhpos[0])+(tpos[1]-bhpos[1])*(tpos[1]-bhpos[1])+(tpos[2]-bhpos[2])*(tpos[2]-bhpos[2]));
      //first assumption: see everything instantly after lighting up. looking from left so only x-component of velocity makes doppler-shift.
      time = rij/clight; //when the TD reaches the gas particle
      freq = (1.0-tv[0]/clight)*basefreq;
      fprintf(fout,"%f\t%f\n",time,freq);
    }

}


int main (int argc, char** argv)
{

  double bhpos[3], clight;
  //  gp** radial;
  bhpos[0] = 0.012866;
  bhpos[1] = 0.3103;
  bhpos[2] = -0.02399;
  clight = 10.0;

  loopthroughgas(bhpos,clight);


  return 0;

}
