#include <stdio.h>
#include <math.h>

int main(int argc, char** argv)
{
  FILE *fin,*fout;
  int histo[100][40];
  char line[100];
  double time,lambda,baselambda;
  int timebin,lambdabin,i,j;
  int accept,dismissal;

  accept=0;
  dismissal=0;

  baselambda = 656.0;

  fin = fopen(argv[1],"r");
  fout = fopen(argv[2],"w");
  
  for (i = 0; i < 100; i++)
    for (j = 0; j < 40; j++)
      histo[i][j] = 0;

  while (fgets(line,100,fin) != NULL)
    {
      sscanf(line, "%lf %lf\n",&time,&lambda);
      //which time bin is it?
      timebin = floor(time*10.0);
      lambdabin = floor((lambda-baselambda)*10.0+20.0);//THINK HERE
      if ((timebin >= 0) && (timebin < 100) && (lambdabin >= 0) && (lambdabin < 40))
	{
	  histo[timebin][lambdabin]++;
	  accept++;
	}
      else
	{
	  dismissal++;
	}
    }
  for (i = 0; i < 100; i++)
    {
      time = i*10.0;
    for (j = 0; j < 40; j++)
      {
	lambda = j+baselambda-20.0;
	fprintf(fout,"%d\t%f\t%d\n",i,lambda,histo[i][j]);
      }
    }
  printf("Accepted %d Dismissed %d\n",accept,dismissal);

  return 0;
}
