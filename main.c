/* Compute time-variable spectrum from the tidal disruption echo in a
circum-binary gas disc as seen from a specific observation angle
Patrick Brem, AEI Potsdam, 07/2012 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define XRES 200
#define YRES 200
#define ZRES 200

#define XMIN (-200)
#define XMAX 200
#define YMIN (-200)
#define YMAX 200
#define ZMIN (-200)
#define ZMAX 50


#define VSCALE 613476.295
#define TIMEINDAYS 23286.24 

void fillarrays(float ***dens, float ***velocity, char* infile, float phi, float theta)
{
  FILE *f;
  int i,j,k,id,accept,discard;
  float tempx,tempy,tempz;
  char line[160];
  float tpos[3],tv[3],mass,rotate[3],volume,columndens;

  theta = theta/180.0*3.14159;
  phi = phi/180.0*3.14159;
  rotate[0] = sin(theta)*sin(phi);
  rotate[1] = sin(theta)*cos(phi);
  rotate[2] = cos(theta);
  discard = 0;
  accept = 0;
  f = fopen(infile, "rt");
  while (fgets(line, 160, f) != NULL)
    {
      sscanf(line, "%d %g %g %g %g %g %g %g\n",&id,&mass,&tpos[0],&tpos[1],&tpos[2],&tv[0],&tv[1],&tv[2]);
      //convert x y and z coordinates into an integer between 0 and X/Y/Z RES
      //first convert into a float between 0 and 1
      tempx = (tpos[0] - XMIN)/(XMAX - XMIN);
      tempy = (tpos[1] - YMIN)/(YMAX - YMIN);
      tempz = (tpos[2] - ZMIN)/(ZMAX - ZMIN);
      if ((tempx < 0.0) || (tempx >= 1.0) || (tempy < 0.0) || (tempy >= 1.0) || (tempz < 0.0) || (tempz >= 1.0)) 
	{
	  discard++;
	  continue;
	}
      //convert it into the array index
      i = floor(tempx*XRES);
      j = floor(tempy*YRES);
      k = floor(tempz*ZRES);
      dens[i][j][k] += mass;
      //ATTENTION: need to divide by volume later to get a real density
      //save memory: need only the line of sight component of the velocity!
      velocity[i][j][k] += mass*(rotate[0]*tv[0]+rotate[1]*tv[1]+rotate[2]*tv[2]);
      //ATTENTION: need to divide by total mass later to get c.o.m. l.o.s. velocity
      accept++;
    }
  //prepare arrays for physical application
  volume = (XMAX - XMIN)*(YMAX - YMIN)*(ZMAX - ZMIN)/(XRES*YRES*ZRES);
  for (i = 0; i < XRES; i++)
    for (j = 0; j < YRES; j++)
      {
	columndens = 0.0;
      for (k = 0; k < ZRES; k++)
	{
	  velocity[i][j][k] /= dens[i][j][k];
	  dens[i][j][k] /= volume;
	  columndens += dens[i][j][k];
	}
      printf("%d\t%d\t%e\n",i,j,columndens);
      }
  //  printf("Accepted %d\tdiscarded %d\n",accept,discard);
}

void loopthroughgas(char* infile, double clight, double phi, double theta)
{
  int id;
  double mass,rij,baselambda,lambda,time,rij_eye;
  FILE *f, *fout, *fout2;
  char line[160];
  double tpos[3],tv[3],bh1pos[3],bh1v[3],bh2pos[3],bh2v[3];
  double rotate[3],v,strength,strength2,rij2,time2;
  baselambda = 656.0;
  f = fopen(infile, "rt");
  fout = fopen("out.dat", "w");
  fout2 = fopen("out2.dat", "w");
  //define rotation matrix; need only to obtain the final z component so only need last line of the matrix, e.g. a vector;
  theta = theta/180.0*3.14159;
  phi = phi/180.0*3.14159;
  rotate[0] = sin(theta)*sin(phi);
  rotate[1] = sin(theta)*cos(phi);
  rotate[2] = cos(theta);
  //read the 2 BHs
  fgets(line,160,f);
  sscanf(line, "%d %lg %lg %lg %lg %lg %lg %lg\n",&id,&mass,&bh1pos[0],&bh1pos[1],&bh1pos[2],&bh1v[0],&bh1v[1],&bh1v[2]);
  fgets(line,160,f);
  sscanf(line, "%d %lg %lg %lg %lg %lg %lg %lg\n",&id,&mass,&bh2pos[0],&bh2pos[1],&bh2pos[2],&bh2v[0],&bh2v[1],&bh2v[2]);
  while (fgets(line, 160, f) != NULL)
    {
      sscanf(line, "%d %lg %lg %lg %lg %lg %lg %lg\n",&id,&mass,&tpos[0],&tpos[1],&tpos[2],&tv[0],&tv[1],&tv[2]);
      rij = sqrt((tpos[0]-bh1pos[0])*(tpos[0]-bh1pos[0])+(tpos[1]-bh1pos[1])*(tpos[1]-bh1pos[1])+(tpos[2]-bh1pos[2])*(tpos[2]-bh1pos[2]));
      rij2 = sqrt((tpos[0]-bh2pos[0])*(tpos[0]-bh2pos[0])+(tpos[1]-bh2pos[1])*(tpos[1]-bh2pos[1])+(tpos[2]-bh2pos[2])*(tpos[2]-bh2pos[2]));
      //need also distance to the observer, who is located at the direction (phi,theta) in degrees
      //rotation so that z points into observer direction (confirm definitions/notations.....)
      rij_eye = rotate[0]*tpos[0]+rotate[1]*tpos[1]+rotate[2]*tpos[2];
      //also project velocity to the line of sight
      v = rotate[0]*tv[0]+rotate[1]*tv[1]+rotate[2]*tv[2];
      time = (rij+rij_eye)/clight*TIMEINDAYS; //when the TD reaches the gas particle and then the radiation reaches the eye
      lambda = baselambda/(1.0-v/clight);
      strength = 1.0/(rij*rij);
      time2 = (rij2+rij_eye)/clight*TIMEINDAYS; //when the TD reaches the gas particle and then the radiation reaches the eye
      strength2 = 1.0/(rij2*rij2);
      //strength is because of the weaker TD impact over distance^2
      fprintf(fout,"%f\t%f\t%f\n",time,lambda,strength);
      fprintf(fout2,"%f\t%f\t%f\n",time2,lambda,strength2);
    }

}


int main (int argc, char** argv)
{
  //this is the 3D density and velocity grids
  float ***dens, ***velocity;
  double clight, phi, theta;
  char* infile;
  int i,j,k;

  infile = (char*) malloc(50*sizeof(char));

  dens = (float ***) malloc(XRES*sizeof(float **));
  for (i = 0; i < XRES; i++)
    {
      dens[i] = (float **) malloc(YRES*sizeof(float*));
      for (j = 0; j < YRES; j++)
	{
	  dens[i][j] = (float *) malloc (ZRES*sizeof(float));
	  for (k = 0; k < ZRES; k++)
	    dens[i][j][k] = 0.0;
	}
    }
  velocity = (float ***) malloc(XRES*sizeof(float **));
  for (i = 0; i < XRES; i++)
    {
      velocity[i] = (float **) malloc(YRES*sizeof(float*));
      for (j = 0; j < YRES; j++)
	{
	  velocity[i][j] = (float *) malloc (ZRES*sizeof(float));
	  for (k = 0; k < ZRES; k++)
	    velocity[i][j][k] = 0.0;
	}
    }
  clight = 3.0e8/VSCALE;
  phi = 0.0;
  theta = 45.0;
  infile = argv[1];
  fillarrays(dens,velocity,infile,phi,theta);
  //  loopthroughgas(infile,clight,phi,theta);


  return 0;

}
