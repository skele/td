/* Compute time-variable spectrum from the tidal disruption echo in a
circum-binary gas disc as seen from a specific observation angle
Patrick Brem, AEI Potsdam, 07/2012 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define XRES 200.0
#define YRES 200.0
#define ZRES 200.0

#define XMIN (-10.0)
#define XMAX 14.0
#define YMIN (-9.0)
#define YMAX 11.0
#define ZMIN (-6.0)
#define ZMAX 10.0

#define NT 10
#define L_0 1
#define DURATION 300
#define SRES 30
#define SMIN 655.3
#define SMAX 657.3

#define OFFSET 10
//additional eye distance so that all the arrival times are > 0

#define VSCALE 613476.295
#define TIMEINDAYS 23286.24 


void fillarrays(double ***dens, double ***velocity, char* infile, double phi, double theta)
{
  FILE *f;
  int i,j,k,id,accept,discard;
  double tempx,tempy,tempz;
  char line[160];
  double tpos[3],tv[3],mass,rotate[3],volume,columndens;

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
      sscanf(line, "%d %lg %lg %lg %lg %lg %lg %lg\n",&id,&mass,&tpos[0],&tpos[1],&tpos[2],&tv[0],&tv[1],&tv[2]);
      //convert x y and z coordinates into an integer between 0 and X/Y/Z RES
      //first convert into a double between 0 and 1
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
  printf("volume %e\n",volume);
  for (i = 0; i < XRES; i++)
    for (j = 0; j < YRES; j++)
      {
	columndens = 0.0;
      for (k = 0; k < ZRES; k++)
	{
	  if (dens[i][j][k] > 0.0)
	    {
	      velocity[i][j][k] /= dens[i][j][k];
	      dens[i][j][k] /= volume;
	      columndens += dens[i][j][k];
	    }
	}
      }
  //  printf("Accepted %d\tdiscarded %d\n",accept,discard);
}

double luminosity(double time)
{
  double temp;
  temp = L_0*pow((time+10.0)/TIMEINDAYS,-5./3.);
  return temp;
}

void loopthroughgrid(double ***dens, double ***velocity, double clight, double phi, double theta, char* infile)
{
  int id;
  double mass,rij,baselambda,lambda,time,rij_eye;
  FILE *f, *fout, *fout2;
  char line[160];
  double tpos[3],bh1pos[3],bh1v[3],bh2pos[3],bh2v[3];
  double rotate[3],rij2,time2;
  baselambda = 656.3;
  int time0_1,time0_2,accepted,rejected;
  int i,j,k,l,freqbin;
  double ts1,ts2,strength0_1,strength0_2;
  double **spectrum,**spectrum2;
  spectrum = (double **) malloc(DURATION * sizeof(double *));
  for (i = 0; i < DURATION; i++)
    spectrum[i] = (double *) malloc(SRES * sizeof(double));
  spectrum2 = (double **) malloc(DURATION * sizeof(double *));
  for (i = 0; i < DURATION; i++)
    spectrum2[i] = (double *) malloc(SRES * sizeof(double));
  for (i = 0; i < DURATION; i++)
    for (j = 0; j < SRES; j++)
      {
	spectrum[i][j] = 0.0;
	spectrum2[i][j] = 0.0;
      }
  f = fopen(infile, "rt");
  accepted = 0;
  rejected = 0;
  fout = fopen("out1.dat", "w");
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

  for (i = 0; i < XRES; i++)
    for (j = 0; j < YRES; j++)
      for (k = 0; k < ZRES; k++)
	{
	  if (dens[i][j][k] > 0.0)
	    {
	      //calculate 3d position of the cell, first make number between 0 and 1
	      tpos[0] = i/((double) XRES)*(XMAX - XMIN) + XMIN;
	      tpos[1] = j/((double) YRES)*(YMAX - YMIN) + YMIN;
	      tpos[2] = k/((double) ZRES)*(ZMAX - ZMIN) + ZMIN;
	      rij = sqrt((tpos[0]-bh1pos[0])*(tpos[0]-bh1pos[0])+(tpos[1]-bh1pos[1])*(tpos[1]-bh1pos[1])+(tpos[2]-bh1pos[2])*(tpos[2]-bh1pos[2]));
	      rij2 = sqrt((tpos[0]-bh2pos[0])*(tpos[0]-bh2pos[0])+(tpos[1]-bh2pos[1])*(tpos[1]-bh2pos[1])+(tpos[2]-bh2pos[2])*(tpos[2]-bh2pos[2]));
	      //printf("%f\t%f\n",rij,rij2);
	      //already know the l.o.s. velocity
	      //need also distance to the observer, who is located at the direction (phi,theta) in degrees
	      //rotation so that z points into observer direction (confirm definitions/notations.....)
	      rij_eye = rotate[0]*tpos[0]+rotate[1]*tpos[1]+rotate[2]*tpos[2];
	      //calculate the light travel time for both BHs
	      time = (rij+rij_eye)/clight*TIMEINDAYS + OFFSET; //when the TD reaches the gas particle and then the radiation reaches the eye
	      time2 = (rij2+rij_eye)/clight*TIMEINDAYS + OFFSET; //when the TD reaches the gas particle and then the radiation reaches the eye
	      //	  printf("%f\t%f\n",time,time2);
	      lambda = baselambda/(1.0-velocity[i][j][k]/clight);
	      //find the time bin where the first entry belongs (t=0 TD emission)
	      time0_1 = floor(time);
	      time0_2 = floor(time2);
	      //find the frequency bin
	      freqbin = floor((lambda-SMIN)/(SMAX - SMIN) * SRES);
	      if ((freqbin >= 0) && (freqbin < SRES))
		accepted++;
	      else
		{
		  //	      printf("lambda %f time %f time2 %f\n",lambda,time,time2);
		  rejected++;
		  continue;
		}
	      strength0_1 = 1.0/(rij*rij);
	      strength0_2 = 1.0/(rij2*rij2);
	      for (l = 0; l < NT; l++)
		{
		  //now calculate the flux for NT different times
		  ts1 = strength0_1*luminosity((double) l)*(dens[i][j][k]);//*(dens[i][j][k]);
		  ts2 = strength0_2*ts1/strength0_1;
		  //printf("%e\t%e\t%e\t%e\t%e\n",ts1,ts2,strength0_1,luminosity((double) l),dens[i][j][k]);
		  if ((time0_1 + l) < DURATION)
		    spectrum[time0_1 + l][freqbin] += ts1;
		  if ((time0_2 + l) < DURATION)
		    spectrum2[time0_2 + l][freqbin] += ts2;
		}
	    }
	}
  for (i = 0; i < DURATION; i++)
    for (j = 0; j < SRES; j++)
      {
	//recover frequency value
	lambda = (double ) j/SRES*(SMAX - SMIN) + SMIN;
	fprintf(fout,"%d\t%f\t%f\n",i,lambda,spectrum[i][j]);
	fprintf(fout2,"%d\t%f\t%f\n",i,lambda,spectrum2[i][j]);
      }
  printf("Accepted %d\tRejected %d\n",accepted,rejected);
}


int main (int argc, char** argv)
{
  //this is the 3D density and velocity grids
  double ***dens, ***velocity;
  double clight, phi, theta;
  char *infile,*infileBH;
  int i,j,k;


  infile = (char *) malloc(100*sizeof(char));
  infileBH = (char *) malloc(100*sizeof(char));
  dens = (double ***) malloc(XRES*sizeof(double **));
  for (i = 0; i < XRES; i++)
    {
      dens[i] = (double **) malloc(YRES*sizeof(double*));
      for (j = 0; j < YRES; j++)
	{
	  dens[i][j] = (double *) malloc (ZRES*sizeof(double));
	  for (k = 0; k < ZRES; k++)
	    dens[i][j][k] = 0.0;
	}
    }
  velocity = (double ***) malloc(XRES*sizeof(double **));
  for (i = 0; i < XRES; i++)
    {
      velocity[i] = (double **) malloc(YRES*sizeof(double*));
      for (j = 0; j < YRES; j++)
	{
	  velocity[i][j] = (double *) malloc (ZRES*sizeof(double));
	  for (k = 0; k < ZRES; k++)
	    velocity[i][j][k] = 0.0;
	}
    }
  clight = 3.0e8/VSCALE;
  phi = 0.0;
  theta = 45.0;
  infile = argv[1];
  infileBH = argv[2];
  fillarrays(dens,velocity,infile,phi,theta);
  loopthroughgrid(dens,velocity,clight,phi,theta,infileBH);


  return 0;

}
