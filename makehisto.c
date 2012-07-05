#include <stdio.h>
#include <math.h>

int main(int argc, char** argv)
{
  FILE *fin,*fout;
  double histo[500][40];
  char line[100];
  double time,lambda,baselambda,strength;
  int timebin,lambdabin,i,j;
  int accept,dismissal;

  accept=0;
  dismissal=0;

  baselambda = 656.0;

  fin = fopen(argv[1],"r");
  fout = fopen(argv[2],"w");
  
  for (i = 0; i < 500; i++)
    for (j = 0; j < 40; j++)
      histo[i][j] = 0.0;

  while (fgets(line,100,fin) != NULL)
    {
      sscanf(line, "%lf %lf %lf\n",&time,&lambda,&strength);
      //which time bin is it?
      timebin = floor(time);
      lambdabin = floor((lambda-baselambda)*20.0+20.0);//THINK HERE
      if ((timebin >= 0) && (timebin < 500) && (lambdabin >= 0) && (lambdabin < 40))
	{
	  histo[timebin][lambdabin] += strength;
	  accept++;
	}
      else
	{
	  dismissal++;
	}
    }
  for (i = 0; i < 500; i++)
    {
      time = i;
    for (j = 0; j < 40; j++)
      {
	lambda = (j-20.0)/20.0+baselambda;
	fprintf(fout,"%d\t%f\t%f\n",i,lambda,histo[i][j]);
      }
    }
  printf("Accepted %d Dismissed %d\n",accept,dismissal);

  return 0;
}
