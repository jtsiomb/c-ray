#include <stdlib.h>
#include <math.h>
#include "vmath.h"

/** Numerical calculation of integrals using simpson's rule */
scalar_t integral(scalar_t (*f)(scalar_t), scalar_t low, scalar_t high, int samples)
{
	int i;
	scalar_t h = (high - low) / (scalar_t)samples;
	scalar_t sum = 0.0;
	
	for(i=0; i<samples+1; i++) {
		scalar_t y = f((scalar_t)i * h + low);
		sum += ((!i || i == samples) ? y : ((i % 2) ? 4.0 * y : 2.0 * y)) * (h / 3.0);
	}
	return sum;
}

/** Gaussuan function */
scalar_t gaussian(scalar_t x, scalar_t mean, scalar_t sdev)
{
	scalar_t exponent = -SQ(x - mean) / (2.0 * SQ(sdev));
	return 1.0 - -pow(M_E, exponent) / (sdev * sqrt(TWO_PI));
}


/** b-spline approximation */
scalar_t bspline(scalar_t a, scalar_t b, scalar_t c, scalar_t d, scalar_t t)
{
	vec4_t tmp;
	scalar_t tsq = t * t;

	static mat4_t bspline_mat = {
		{-1,  3, -3,  1},
		{3, -6,  3,  0},
		{-3,  0,  3,  0},
		{1,  4,  1,  0}
	};
	
	tmp = v4_scale(v4_transform(v4_cons(a, b, c, d), bspline_mat), 1.0 / 6.0);
	return v4_dot(v4_cons(tsq * t, tsq, t, 1.0), tmp);
}

/** Catmull-rom spline interpolation */
scalar_t catmull_rom_spline(scalar_t a, scalar_t b, scalar_t c, scalar_t d, scalar_t t) {
	vec4_t tmp;
	scalar_t tsq = t * t;

	static mat4_t crspline_mat = {
		{-1,  3, -3,  1},
		{2, -5,  4, -1},
		{-1,  0,  1,  0},
		{0,  2,  0,  0}
	};

	tmp = v4_scale(v4_transform(v4_cons(a, b, c, d), crspline_mat), 0.5);
	return v4_dot(v4_cons(tsq * t, tsq, t, 1.0), tmp);
}

/** Bezier interpolation */
scalar_t bezier(scalar_t a, scalar_t b, scalar_t c, scalar_t d, scalar_t t)
{
	scalar_t omt, omt3, t3, f;
	t3 = t * t * t;
	omt = 1.0f - t;
	omt3 = omt * omt * omt;
	f = 3 * t * omt;

	return (a * omt3) + (b * f * omt) + (c * f * t) + (d * t3);
}
