#include "randomize.hpp"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>


static unsigned long long g_base_seed = 0;

/* Per-thread state: xorshift64* returns the same 31-bit range as rand(). */

static thread_local unsigned long long t_state = 0;


void
random_init()
{
  FILE *f=NULL;
  unsigned long long seed = (unsigned long long)time(NULL);
  if((f=fopen("/dev/random","rb"))!=NULL)
    {
      fread(&seed, sizeof(seed), 1, f);
      fclose(f);
    }
  g_base_seed = seed;
  dda_srand(seed);
}


unsigned long long
random_base_seed()
{
  return g_base_seed;
}


void
dda_srand(unsigned long long seed)
{
  t_state = seed ? seed : 0x9E3779B97F4A7C15ULL;
}


int
dda_rand()
{
  t_state ^= t_state >> 12;
  t_state ^= t_state << 25;
  t_state ^= t_state >> 27;
  return (int)((t_state * 2685821657736338717ULL) >> 33);
}
