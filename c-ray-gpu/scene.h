#ifndef _SCENE_H_
#define _SCENE_H_

#include "vmath.h"

#define MAX_LIGHTS		16				/* maximum number of lights */

struct ray {
	vec3_t orig, dir;
};

struct material {
	vec3_t col;		/* color */
	float spow;		/* specular power */
	float refl;		/* reflection intensity */
};

struct sphere {
	vec3_t pos;
	float rad;
	struct material mat;
	struct sphere *next;
};

struct spoint {
	vec3_t pos, normal, vref;	/* position, normal and view reflection */
	float dist;		/* parametric distance of intersection along the ray */
};

struct camera {
	vec3_t pos, targ;
	float fov;
};

struct sphere *obj_list;
int obj_num;
struct camera cam;

void load_scene(FILE *fp);

#endif	/* _SCENE_H_ */
