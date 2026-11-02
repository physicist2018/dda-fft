#include "randomize.hpp"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>


void
random_init()
{
  FILE *f=NULL;
  time_t seed = time(NULL);
  if((f=fopen("/dev/random","rb"))!=NULL)
    {
      fread(&seed, sizeof(seed), 1, f);
      fclose(f);
    }
  else
    {
      seed=time(NULL);
    }
  srand(seed);    
}

