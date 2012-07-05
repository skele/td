#include <stdio.h>
#include <math.h>

int main(int argc, char** argv)
{
  FILE *fin,*fout;
  int histo[100][40];
  char line[100];
  double time,freq,basefreq;
  int timebin,freqbin,i,j;
  int accept,dismissal;

  accept=0;
  dismissal=0;

  basefreq = 1000.0;

  fin = fopen(argv[1],"r");
  fout = fopen(argv[2],"w");
  
  for (i = 0; i < 100; i++)
    for (j = 0; j < 40; j++)
      histo[i][j] = 0;

  while (fgets(line,100,fin) != NULL)
    {
      sscanf(line, "%lf %lf\n",&time,&freq);
      //which time bin is it?
      timebin = floor(time);
      freqbin = floor((freq-basefreq)/10.0+20.0);
      if ((timebin >= 0) && (timebin < 100) && (freqbin >= 0) && (freqbin < 40))
	{
	  histo[timebin][freqbin]++;
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
	freq = j*10.0+basefreq-200.0;
	fprintf(fout,"%d\t%f\t%d\n",i,freq,histo[i][j]);
      }
    }
  printf("Accepted %d Dismissed %d\n",accept,dismissal);

  return 0;
}
