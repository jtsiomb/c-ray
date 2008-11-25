#include <stdlib.h>

/** Generates a random number in [0, range) */
static inline scalar_t frand(scalar_t range)
{
	return range * (float)rand() / (float)RAND_MAX;
}


/** linear interpolation */
static inline scalar_t lerp(scalar_t a, scalar_t b, scalar_t t)
{
	return a + (b - a) * t;
}
