/* Compute time-variable spectrum from the tidal disruption echo in a
circum-binary gas disc as seen from a specific observation angle
Patrick Brem, AEI Potsdam, 07/2012 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PHIRES 200.0
#define SINTHETARES 200.0
#define RRES 100.0

#define PI 3.14159

#define PHIMIN (0.0)
#define PHIMAX (2.0*PI)
#define SINTHETAMIN (-1.0)
#define SINTHETAMAX (1.0)
#define RMIN (1.0)
#define RMAX 200.0

#define NT 50

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

const double TEMPCOEFF = 0.000000001;

const double VSCALE = sqrt(G_SI*KG_SCALE/M_SCALE);
const double TIMEINSECONDS = sqrt(M_SCALE*M_SCALE*M_SCALE/G_SI/KG_SCALE);
const double TIMEINDAYS = TIMEINSECONDS/(24.0*3600.0);

const double PROTONMASS = 1.67E-27 / KG_SCALE;

//define ionization rate value in SI units

const double I_0_SI = 3.8E+31;
const double I_0 = I_0_SI/(M_SCALE*M_SCALE/TIMEINSECONDS);
const double HV0 = 2.17E-11; //photon energy in erg

const double DL = 3.1E+24;

int getsinthetaindex(double theta)
{
  /* Determine an index between 0 and SINTHETARES so that sin(theta) is equally spaced wiht resolution SINTHETARES*/
  double temp;
  int index;
  temp = sin(theta);
  index = floor((temp+1.0)/2.0*SINTHETARES);
  return index;
}

void fillarrays(double ***dens, double ***velocity, char* infile, char* bhfile, double phi, double theta)
{
  FILE *f,*fbh;
  int i,j,k,id,accept,discard;
  double tempx,tempy,tempz;
  char line[160];
  double tpos[3],tv[3],mass,rotate[3],volume,columndens;
  double bhpos[3],bhv[3],tempr,tempphi,temptheta;
  theta = theta/180.0*3.14159;
  phi = phi/180.0*3.14159;
  rotate[0] = sin(theta)*sin(phi);
  rotate[1] = sin(theta)*cos(phi);
  rotate[2] = cos(theta);
  discard = 0;
  accept = 0;
  f = fopen(infile, "rt");
  fbh = fopen(bhfile, "rt");
  //read BHs and take first BH as center for spherical coordinates
  fgets(line,160,fbh);
  sscanf(line, "%d %lg %lg %lg %lg %lg %lg %lg\n",&id,&mass,&bhpos[0],&bhpos[1],&bhpos[2],&bhv[0],&bhv[1],&bhv[2]);

  while (fgets(line, 160, f) != NULL)
    {
      sscanf(line, "%d %lg %lg %lg %lg %lg %lg %lg\n",&id,&mass,&tpos[0],&tpos[1],&tpos[2],&tv[0],&tv[1],&tv[2]);
      //compute relative distance and velocity to BH
      tempx = (tpos[0]-bhpos[0]);
      tempy = (tpos[1]-bhpos[1]);
      tempz = (tpos[2]-bhpos[2]);
      //transform to spherical coordinates phi,theta,r
      tempr = sqrt(tempx*tempx+tempy*tempy+tempz*tempz);
      tempphi = atan2(tempy,tempx)+PI;
      temptheta = asin(tempz/tempr);
      //check whether r is in the range
      if ((tempr > RMAX) || (tempr < RMIN))
	{
	  discard++;
	  continue;
	}
      //convert it into array index
      i = floor((tempphi/(2.0*PI))*PHIRES);
      j = getsinthetaindex(temptheta);
      k = floor((tempr-RMIN)/(RMAX-RMIN)*RRES);
      printf("tempphi %f\ttemptheta %f\n",tempphi,temptheta);
      printf("Indeces: %d\t%d\t%d\n",i,j,k);
      dens[i][j][k] += mass;
      //ATTENTION: need to divide by volume later to get a real density
      //save memory: need only the line of sight component of the velocity!
      velocity[i][j][k] += mass*(rotate[0]*tv[0]+rotate[1]*tv[1]+rotate[2]*tv[2]);
      //ATTENTION: need to divide by total mass later to get c.o.m. l.o.s. velocity
      accept++;
    }
  //prepare arrays for physical application
  
  /* volume = (XMAX - XMIN)*(YMAX - YMIN)*(ZMAX - ZMIN)/(PHIRES*YRES*RRES); */
  /* printf("volume %e\n",volume); */
  for (i = 0; i < PHIRES; i++)
    for (j = 0; j < SINTHETARES; j++)
      {
	columndens = 0.0;
      for (k = 0; k < RRES; k++)
	{
	  if (dens[i][j][k] > 0.0)
	    {
	      velocity[i][j][k] /= dens[i][j][k];
	      //dens[i][j][k] /= volume;
	      //without the /= volume this is the absolute number of mass in a cell
	      columndens += dens[i][j][k];
	    }
	}
      }
      printf("Boxed %2.2f%% of the particles\n",(double) accept*100.0/(accept+discard));
}

double luminosity(double time, double t_min)
{
  double temp;
  if (time > t_min)
    {
      temp = I_0*pow((time)/t_min,-5./3.); 
    }
  else
    {
      temp = 0.0;
    }
  return temp;
}

void loopthroughgrid(double ***dens, double ***velocity, double clight, double phi, double theta, char* infile)
{
  int id;
  double mass1,mass2,rij,baselambda,lambda,dlambda1,dlambda2,time,rij_eye,temperature1,temperature2,base_ts1,base_ts2;
  FILE *f, *fout, *fout2;
  char line[160];
  double tpos[3],bh1pos[3],bh1v[3],bh2pos[3],bh2v[3];
  double rotate[3],rij2,time2,difflambda;
  baselambda = 656.3;
  int time0_1,time0_2,accepted,rejected,intstarttime;
  int i,j,k,l,ll,freqbin;
  double ts1,ts2,strength0_1,strength0_2;
  double **spectrum,**spectrum2,tcoeff1,tcoeff2;
  double t_min1,t_min2;
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

  //derive the shape of the tidal flare depending on BH mass (Ulmer 1997)
  //order of magnitude estimate in DAYS
  t_min1 = 0.11*sqrt(mass1*MSUN_SCALE/(1.0E+06))*365.0;
  t_min2 = 0.11*sqrt(mass2*MSUN_SCALE/(1.0E+06))*365.0;
  //earliest time when the tidal flare leaves the BH;
  intstarttime = floor(t_min1);
  if (intstarttime > floor(t_min2))
    {
      intstarttime = floor(t_min2);
    }
  printf("T_min1 = %e days\nT_min2 = %e days\n",t_min1,t_min2);

  for (i = 0; i < PHIRES; i++)
    for (j = 0; j < SINTHETARES; j++)
      for (k = 0; k < RRES; k++)
	{
	  if (dens[i][j][k] > 0.0)
	    {
	      //calculate 3d position of the cell, first make number between 0 and 1 then also add the half length so that it is in the center of the cell
	      /*
	      tpos[0] = i/((double) PHIRES)*(XMAX - XMIN) + XMIN + (XMAX-XMIN)/PHIRES/2.0;
	      tpos[1] = j/((double) YRES)*(YMAX - YMIN) + YMIN + (YMAX-YMIN)/YRES/2.0;
	      tpos[2] = k/((double) RRES)*(ZMAX - ZMIN) + ZMIN + (ZMAX-ZMIN)/RRES/2.0;
	      rij = sqrt((tpos[0]-bh1pos[0])*(tpos[0]-bh1pos[0])+(tpos[1]-bh1pos[1])*(tpos[1]-bh1pos[1])+(tpos[2]-bh1pos[2])*(tpos[2]-bh1pos[2]));
	      rij2 = sqrt((tpos[0]-bh2pos[0])*(tpos[0]-bh2pos[0])+(tpos[1]-bh2pos[1])*(tpos[1]-bh2pos[1])+(tpos[2]-bh2pos[2])*(tpos[2]-bh2pos[2]));
	      */	      


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
	      for (l = intstarttime; l < (intstarttime + NT); l++)
		{
		  //now calculate the flux for NT different times
		  //this is ergs emitted per unit time
		  base_ts1 = strength0_1*luminosity((double) l,t_min1)*(dens[i][j][k])/PROTONMASS*HV0;//*(dens[i][j][k]);
		  base_ts2 = strength0_2*luminosity((double) l,t_min2)*(dens[i][j][k])/PROTONMASS*HV0;
		  //calculate the flux at some distance given in cm now all in physical units
		  /*		  base_ts1 = base_ts1/((16*PI*PI)*DL*DL);
		  base_ts2 = base_ts2/((16*PI*PI)*DL*DL);
		  */
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
  printf("1 CODE Ionization rate =\t%e = %e m^2/s\n",I_0,I_0_SI);

  infile = (char *) malloc(100*sizeof(char));
  infileBH = (char *) malloc(100*sizeof(char));
  dens = (double ***) malloc(PHIRES*sizeof(double **));
  for (i = 0; i < PHIRES; i++)
    {
      dens[i] = (double **) malloc(SINTHETARES*sizeof(double*));
      for (j = 0; j < SINTHETARES; j++)
	{
	  dens[i][j] = (double *) malloc (RRES*sizeof(double));
	  for (k = 0; k < RRES; k++)
	    dens[i][j][k] = 0.0;
	}
    }
  velocity = (double ***) malloc(PHIRES*sizeof(double **));
  for (i = 0; i < PHIRES; i++)
    {
      velocity[i] = (double **) malloc(SINTHETARES*sizeof(double*));
      for (j = 0; j < SINTHETARES; j++)
	{
	  velocity[i][j] = (double *) malloc (RRES*sizeof(double));
	  for (k = 0; k < RRES; k++)
	    velocity[i][j][k] = 0.0;
	}
    }
  clight = 3.0e8/VSCALE;
  phi = 0.0;
  theta = 45.0;
  infile = argv[1];
  infileBH = argv[2];
  fillarrays(dens,velocity,infile,infileBH,phi,theta);
  loopthroughgrid(dens,velocity,clight,phi,theta,infileBH);


  return 0;

}
