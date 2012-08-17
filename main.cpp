/* Compute time-variable spectrum from the tidal disruption echo in a
circum-binary gas disc as seen from a specific observation angle
Patrick Brem, AEI Potsdam, 07/2012 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define XRES 200.0
#define YRES 200.0
#define ZRES 100.0

#define XMIN (-20.0)
#define XMAX 20.0
#define YMIN (-20.0)
#define YMAX 20.0
#define ZMIN (-10.0)
#define ZMAX 10.0

#define NT 20
#define L_0 1.0
#define DURATION 300
#define SRES 30
#define SMIN 655.3
#define SMAX 657.3

#define OFFSET 20
//additional eye distance so that all the arrival times are > 0

// define what 1 code unit is in length and mass
const double PC_SCALE = (0.04);
const double MSUN_SCALE = (3.5E+06);

const double M_SCALE = (PC_SCALE*3.086E+016);
const double KG_SCALE = (1.99E+30*MSUN_SCALE);
const double G_SI = (6.67E-11);
//derive velocity scale

const double TEMPCOEFF = 0.000001;

const double VSCALE = sqrt(G_SI*KG_SCALE/M_SCALE);
const double TIMEINSECONDS = sqrt(M_SCALE*M_SCALE*M_SCALE/G_SI/KG_SCALE);
const double TIMEINDAYS = TIMEINSECONDS/(24.0*3600.0);

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
  printf("Boxed %2.2f%% of the particles\n",(double) accept*100.0/(accept+discard));
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
  double mass,rij,baselambda,lambda,dlambda1,dlambda2,time,rij_eye,temperature1,temperature2,base_ts1,base_ts2;
  FILE *f, *fout, *fout2;
  char line[160];
  double tpos[3],bh1pos[3],bh1v[3],bh2pos[3],bh2v[3];
  double rotate[3],rij2,time2,difflambda;
  baselambda = 656.3;
  int time0_1,time0_2,accepted,rejected;
  int i,j,k,l,ll,freqbin;
  double ts1,ts2,strength0_1,strength0_2;
  double **spectrum,**spectrum2,tcoeff1,tcoeff2;
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
	      //calculate 3d position of the cell, first make number between 0 and 1 then also add the half length so that it is in the center of the cell
	      tpos[0] = i/((double) XRES)*(XMAX - XMIN) + XMIN + (XMAX-XMIN)/XRES/2.0;
	      tpos[1] = j/((double) YRES)*(YMAX - YMIN) + YMIN + (YMAX-YMIN)/YRES/2.0;
	      tpos[2] = k/((double) ZRES)*(ZMAX - ZMIN) + ZMIN + (ZMAX-ZMIN)/ZRES/2.0;
	      rij = sqrt((tpos[0]-bh1pos[0])*(tpos[0]-bh1pos[0])+(tpos[1]-bh1pos[1])*(tpos[1]-bh1pos[1])+(tpos[2]-bh1pos[2])*(tpos[2]-bh1pos[2]));
	      rij2 = sqrt((tpos[0]-bh2pos[0])*(tpos[0]-bh2pos[0])+(tpos[1]-bh2pos[1])*(tpos[1]-bh2pos[1])+(tpos[2]-bh2pos[2])*(tpos[2]-bh2pos[2]));
	      


	      //already know the l.o.s. velocity
	      //need also distance to the observer, who is located at the direction (phi,theta) in degrees
	      //rotation so that z points into observer direction (confirm definitions/notations.....)
	      rij_eye = rotate[0]*tpos[0]+rotate[1]*tpos[1]+rotate[2]*tpos[2];
	      //calculate the light travel time for both BHs
	      time = (rij+rij_eye)/clight*TIMEINDAYS + OFFSET; //when the TD reaches the gas particle and then the radiation reaches the eye
	      time2 = (rij2+rij_eye)/clight*TIMEINDAYS + OFFSET; //when the TD reaches the gas particle and then the radiation reaches the eye
	      //	  printf("%f\t%f\n",time,time2);
	      lambda = baselambda/(1.0-velocity[i][j][k]/clight);

	      //rij2 will now also be used to assign a certain temperature for doppler broadening
	      temperature1 = TEMPCOEFF/rij/rij;
	      dlambda1 = lambda*sqrt(temperature1); // the broadening
	      temperature2 = TEMPCOEFF/rij2/rij2;
	      dlambda2 = lambda*sqrt(temperature2);// the broadening

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
		  base_ts1 = strength0_1*luminosity((double) l)*(dens[i][j][k]);//*(dens[i][j][k]);
		  base_ts2 = strength0_2*ts1/strength0_1;

		  //loop over whole frequency and calculate corresponding doppler broadened flux value
		  for (ll = 0; ll < SRES; ll++)
		    {
		      //convert the difference in frequency bins into a difference in lambda
		      difflambda = (freqbin-ll)/((float) SRES)*(SMAX - SMIN);
		      tcoeff1 = 1.0/dlambda1*exp(-(difflambda*difflambda)/(dlambda1*dlambda1));
		      tcoeff2 = 1.0/dlambda2*exp(-(difflambda*difflambda)/(dlambda2*dlambda2));
		      ts1 = base_ts1*tcoeff1;
		      ts2 = base_ts2*tcoeff2;
		      // printf("diff lambda: %e strength coeff: %e temp: %e  dlambda: %e\n ",difflambda,tcoeff1,temperature1,dlambda1);
		      if (((time0_1 + l) < DURATION) && ((time0_1 + l) >= 0))
			spectrum[time0_1 + l][ll] += ts1;
		      if (((time0_2 + l) < DURATION) && ((time0_2 + l) >= 0))
			spectrum2[time0_2 + l][ll] += ts2;
		      
		      
		    }
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
  printf("Accepted %2.2f%% of time and spectrum data\n",(double) accepted*100.0/(accepted+rejected));
}


int main (int argc, char** argv)
{
  //this is the 3D density and velocity grids
  double ***dens, ***velocity;
  double clight, phi, theta;
  char *infile,*infileBH;
  int i,j,k;

  if (argc != 3)
    {
      printf("Usage: \n");
      printf("./td <inputfile_gas> <inputfile_BHs>\n");
      return 0;
    }

  printf("Scales:\n1 CODE MASS =\t%e kg = %e Msun\n",KG_SCALE,MSUN_SCALE);
  printf("1 CODE LENGTH =\t%e m = %e pc\n",M_SCALE,PC_SCALE);
  printf("1 CODE VEL =\t%e m/s\n",VSCALE);
  printf("1 CODE TIME =\t%e sec = %e days\n",TIMEINSECONDS,TIMEINDAYS);

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
