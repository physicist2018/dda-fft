#ifndef RANDOMIZE_H
#define RANDOMIZE_H

/* The C library rand() keeps its state in a global variable, so it cannot be
   used once several threads generate samples at the same time. Every worker
   thread gets its own generator instead; dda_rand() returns values in exactly
   the same range as rand(), so all the "(real)dda_rand()/DDA_RAND_MAX" and
   "(real)dda_rand()*2./DDA_RAND_MAX-1." expressions are unchanged. */

#define DDA_RAND_MAX 2147483647

void
random_init(void);

/* Same as random_init(), but with a caller supplied base seed, so that a run
   can be repeated exactly (--seed). */

void
random_init_seed(unsigned long long seed);

unsigned long long
random_base_seed(void);

void
dda_srand(unsigned long long seed);

int
dda_rand(void);

#endif
