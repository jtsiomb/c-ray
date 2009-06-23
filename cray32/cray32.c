#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <float.h>
#include "gfx.h"
#include "inttypes.h"

struct vec3 {
	float x, y, z;
};

struct ray {
	struct vec3 orig, dir;
};

struct material {
	struct vec3 col;	/* color */
	float spow;		/* specular power */
	float refl;		/* reflection intensity */
};

struct sphere {
	struct vec3 pos;
	float rad;
};

struct plane {
	struct vec3 pos;
	struct vec3 norm;
};

#define OBJ_SPHERE		1
#define OBJ_PLANE		2

struct object {
	int type;
	union {
		struct sphere sph;
		struct plane pln;
	} o;
	struct material mat;

	struct object *next;
};

struct spoint {
	struct vec3 pos, normal, vref;	/* position, normal and view reflection */
	float dist;		/* parametric distance of intersection along the ray */
};

struct camera {
	struct vec3 pos, targ;
	float fov;
};

void render(int xsz, int ysz, uint32_t *fb, int samples);
struct vec3 trace(struct ray ray, int depth);
struct vec3 shade(struct object *obj, struct spoint *sp, int depth);
struct vec3 reflect(struct vec3 v, struct vec3 n);
struct vec3 cross_product(struct vec3 v1, struct vec3 v2);
struct ray get_primary_ray(int x, int y, int sample);
struct vec3 get_sample_pos(int x, int y, int sample);
struct vec3 jitter(int x, int y, int s);
int ray_sphere(const struct sphere *sph, struct ray ray, struct spoint *sp);
void load_scene(FILE *fp);

#define MAX_LIGHTS		16				/* maximum number of lights */
#define RAY_MAG			1000.0			/* trace rays of this magnitude */
#define MAX_RAY_DEPTH	4				/* raytrace recursion limit */
#define FOV				0.78539816		/* field of view in rads (pi/4) */
#define HALF_FOV		(FOV * 0.5)
#define ERR_MARGIN		1e-6			/* an arbitrary error margin to avoid surface acne */

/* some helpful macros... */
#define SQ(x)		((x) * (x))
#define MAX(a, b)	((a) > (b) ? (a) : (b))
#define MIN(a, b)	((a) < (b) ? (a) : (b))
#define DOT(a, b)	((a).x * (b).x + (a).y * (b).y + (a).z * (b).z)
#define NORMALIZE(a)  do {\
	float len = sqrt(DOT(a, a));\
	(a).x /= len; (a).y /= len; (a).z /= len;\
} while(0);

/* global state */
int xres = 640;
int yres = 480;
int bpp = 32;
int rshift, gshift, bshift;
float aspect;
struct object *obj_list;
struct vec3 lights[MAX_LIGHTS];
int lnum = 0;
struct camera cam;

#define NRAN	1024
#define MASK	(NRAN - 1)
struct vec3 urand[NRAN];
int irand[NRAN];


const char *usage = {
	"Usage: cray32 [options]\n"
	"Options:\n"
	"  -s WxH     where W is the width and H the height of the image\n"
	"  -r <rays>  shoot <rays> rays per pixel (antialiasing)\n"
	"  -i <file>  read from <file> instead of stdin\n"
	"  -h         this help screen\n\n"
	"Note: at the moment the only valid resolutions are: 320x240 and 320x200\n\n"
};


int main(int argc, char **argv) {
	int i;
	uint32_t *pixels;
	int rays_per_pixel = 1;
	FILE *infile = stdin;

	for(i=1; i<argc; i++) {
		if(argv[i][0] == '-' && argv[i][2] == 0) {
			char *sep;
			switch(argv[i][1]) {
			case 's':
				if(!isdigit(argv[++i][0]) || !(sep = strchr(argv[i], 'x')) || !isdigit(*(sep + 1))) {
					fputs("-s must be followed by something like \"640x480\"\n", stderr);
					return 1;
				}
				xres = atoi(argv[i]);
				yres = atoi(sep + 1);
				break;

			case 'i':
				if(!(infile = fopen(argv[++i], "rb"))) {
					fprintf(stderr, "failed to open input file %s: %s\n", argv[i], strerror(errno));
					return 1;
				}
				break;

			case 'r':
				if(!isdigit(argv[++i][0])) {
					fputs("-r must be followed by a number (rays per pixel)\n", stderr);
					return EXIT_FAILURE;
				}
				rays_per_pixel = atoi(argv[i]);
				break;


			case 'h':
				fputs(usage, stdout);
				return 0;
				
			default:
				fprintf(stderr, "unrecognized argument: %s\n", argv[i]);
				fputs(usage, stderr);
				return 1;
			}
		} else {
			fprintf(stderr, "unrecognized argument: %s\n", argv[i]);
			fputs(usage, stderr);
			return 1;
		}
	}

	aspect = (float)xres / (float)yres;

	/*
	if(!(pixels = malloc(xres * yres * sizeof *pixels))) {
		perror("pixel buffer allocation failed");
		return 1;
	}
	*/

	load_scene(infile);
	if(infile != stdin) fclose(infile);

	/* initialize the random number tables for the jitter */
	for(i=0; i<NRAN; i++) urand[i].x = (float)rand() / RAND_MAX - 0.5;
	for(i=0; i<NRAN; i++) urand[i].y = (float)rand() / RAND_MAX - 0.5;
	for(i=0; i<NRAN; i++) irand[i] = (int)(NRAN * ((float)rand() / RAND_MAX));


	printf("initializing graphics: %dx%d %dbpp\n", xres, yres, bpp);
	if(!(pixels = set_video_mode(xres, yres, bpp))) {
		return 1;
	}
	if((bpp = get_color_depth()) == 32) {
	    get_color_shift(&rshift, &gshift, &bshift);
		printf("got 32bpp (r << %d, g << %d, b << %d)\n", rshift, gshift, bshift);
	} else {
		printf("got 24bpp\n");
	}

	render(xres, yres, pixels, rays_per_pixel);

	for(;;) {
		char c = getch();

		if(c == 'q' || c == 27) break;
	}

	set_text_mode();
	return 0;
}

#define SCAN_SKIP	8
void render(int xsz, int ysz, uint32_t *fb, int samples) {
	int i, x, y, s;
	uint32_t *p32;
	uint8_t *p24;
	float rcp_samples = 1.0 / (float)samples;

	for(i=0; i<SCAN_SKIP; i++) {
		for(y=i; y<ysz; y+=SCAN_SKIP) {
			if(bpp == 32) {
			    p32 = fb + y * xsz;
			} else {
			    p24 = (unsigned char*)fb + y * xsz * 3;
			}

			for(x=0; x<xsz; x++) {
				unsigned int cr, cg, cb;
				float r, g, b;
				r = g = b = 0.0;
			
				for(s=0; s<samples; s++) {
					struct vec3 col = trace(get_primary_ray(x, y, s), 0);
					r += col.x;
					g += col.y;
					b += col.z;
				}

				cr = (unsigned int)(r * rcp_samples * 255.0);
				cg = (unsigned int)(g * rcp_samples * 255.0);
				cb = (unsigned int)(b * rcp_samples * 255.0);

				if(cr > 255) cr = 255;
				if(cg > 255) cg = 255;
				if(cb > 255) cb = 255;

				if(bpp == 32) {
				    *p32++ = (cr << rshift) | (cg << gshift) | (cb << bshift);
				} else {
				    *p24++ = cr;
				    *p24++ = cg;
				    *p24++ = cb;
				}
			}

			/* TODO use interrupts to avoid polling */
			while(kbhit()) {
				if(getch() == 27) {
					return;
				}
			}
		}
	}
}

/* trace a ray throught the scene recursively (the recursion happens through
 * shade() to calculate reflection rays if necessary).
 */
struct vec3 trace(struct ray ray, int depth) {
	struct vec3 col;
	struct spoint sp, nearest_sp;
	struct object *nearest_obj = 0;
	struct object *iter = obj_list->next;

	/* if we reached the recursion limit, bail out */
	if(depth >= MAX_RAY_DEPTH) {
		col.x = col.y = col.z = 0.0;
		return col;
	}
	
	/* find the nearest intersection ... */
	while(iter) {
		if(ray_object(iter, ray, &sp)) {
			if(!nearest_obj || sp.dist < nearest_sp.dist) {
				nearest_obj = iter;
				nearest_sp = sp;
			}
		}
		iter = iter->next;
	}

	/* and perform shading calculations as needed by calling shade() */
	if(nearest_obj) {
		col = shade(nearest_obj, &nearest_sp, depth);
	} else {
		col.x = col.y = col.z = 0.0;
	}

	return col;
}

/* Calculates direct illumination with the phong reflectance model.
 * Also handles reflections by calling trace again, if necessary.
 */
struct vec3 shade(struct object *obj, struct spoint *sp, int depth) {
	int i;
	struct vec3 col = {0.05, 0.05, 0.05};

	/* for all lights ... */
	for(i=0; i<lnum; i++) {
		float ispec, idiff;
		struct vec3 ldir;
		struct ray shadow_ray;
		struct object *iter = obj_list->next;
		int in_shadow = 0;

		ldir.x = lights[i].x - sp->pos.x;
		ldir.y = lights[i].y - sp->pos.y;
		ldir.z = lights[i].z - sp->pos.z;

		shadow_ray.orig = sp->pos;
		shadow_ray.dir = ldir;

		/* shoot shadow rays to determine if we have a line of sight with the light */
		while(iter) {
			if(ray_object(iter, shadow_ray, 0)) {
				in_shadow = 1;
				break;
			}
			iter = iter->next;
		}

		/* and if we're not in shadow, calculate direct illumination with the phong model. */
		if(!in_shadow) {
			NORMALIZE(ldir);

			idiff = MAX(DOT(sp->normal, ldir), 0.0);
			ispec = obj->mat.spow > 0.0 ? pow(MAX(DOT(sp->vref, ldir), 0.0), obj->mat.spow) : 0.0;

			col.x += idiff * obj->mat.col.x + ispec;
			col.y += idiff * obj->mat.col.y + ispec;
			col.z += idiff * obj->mat.col.z + ispec;
		}
	}

	/* Also, if the object is reflective, spawn a reflection ray, and call trace()
	 * to calculate the light arriving from the mirror direction.
	 */
	if(obj->mat.refl > 0.0) {
		struct ray ray;
		struct vec3 rcol;

		ray.orig = sp->pos;
		ray.dir = sp->vref;
		ray.dir.x *= RAY_MAG;
		ray.dir.y *= RAY_MAG;
		ray.dir.z *= RAY_MAG;

		rcol = trace(ray, depth + 1);
		col.x += rcol.x * obj->mat.refl;
		col.y += rcol.y * obj->mat.refl;
		col.z += rcol.z * obj->mat.refl;
	}

	return col;
}

/* calculate reflection vector */
struct vec3 reflect(struct vec3 v, struct vec3 n) {
	struct vec3 res;
	float dot = v.x * n.x + v.y * n.y + v.z * n.z;
	res.x = -(2.0 * dot * n.x - v.x);
	res.y = -(2.0 * dot * n.y - v.y);
	res.z = -(2.0 * dot * n.z - v.z);
	return res;
}

struct vec3 cross_product(struct vec3 v1, struct vec3 v2) {
	struct vec3 res;
	res.x = v1.y * v2.z - v1.z * v2.y;
	res.y = v1.z * v2.x - v1.x * v2.z;
	res.z = v1.x * v2.y - v1.y * v2.x;
	return res;
}

/* determine the primary ray corresponding to the specified pixel (x, y) */
struct ray get_primary_ray(int x, int y, int sample) {
	struct ray ray;
	float m[3][3];
	struct vec3 i, j = {0, 1, 0}, k, dir, orig = {0, 0, 0};

	k.x = cam.targ.x - cam.pos.x;
	k.y = cam.targ.y - cam.pos.y;
	k.z = cam.targ.z - cam.pos.z;
	NORMALIZE(k);

	i = cross_product(j, k);
	j = cross_product(k, i);
	m[0][0] = i.x; m[0][1] = j.x; m[0][2] = k.x;
	m[1][0] = i.y; m[1][1] = j.y; m[1][2] = k.y;
	m[2][0] = i.z; m[2][1] = j.z; m[2][2] = k.z;
	
	ray.dir = get_sample_pos(x, y, sample);
	ray.dir.z = 1.0 / HALF_FOV;
	ray.dir.x *= RAY_MAG;
	ray.dir.y *= RAY_MAG;
	ray.dir.z *= RAY_MAG;
	
	dir = ray.dir;
	ray.dir.x = dir.x * m[0][0] + dir.y * m[0][1] + dir.z * m[0][2];
	ray.dir.y = dir.x * m[1][0] + dir.y * m[1][1] + dir.z * m[1][2];
	ray.dir.z = dir.x * m[2][0] + dir.y * m[2][1] + dir.z * m[2][2];

	ray.orig.x = orig.x * m[0][0] + orig.y * m[0][1] + orig.z * m[0][2] + cam.pos.x;
	ray.orig.y = orig.x * m[1][0] + orig.y * m[1][1] + orig.z * m[1][2] + cam.pos.y;
	ray.orig.z = orig.x * m[2][0] + orig.y * m[2][1] + orig.z * m[2][2] + cam.pos.z;

	ray.dir.x += orig.x;
	ray.dir.y += orig.y;
	ray.dir.z += orig.z;
	
	return ray;
}


struct vec3 get_sample_pos(int x, int y, int sample) {
	struct vec3 pt;
	static float sf = 0.0;

	if(sf == 0.0) {
		sf = 1.5 / (float)xres;
	}

	pt.x = ((float)x / (float)xres - 0.5) * 2.0 * aspect;
	pt.y = (0.5 - (float)y / (float)yres) * 2.0;

	if(sample) {
		struct vec3 jt = jitter(x, y, sample);
		pt.x += jt.x * sf;
		pt.y += jt.y * sf / aspect;
	}
	return pt;
}

/* jitter function taken from Graphics Gems I. */
struct vec3 jitter(int x, int y, int s) {
	struct vec3 pt;
	pt.x = urand[(x + (y << 2) + irand[(x + s) & MASK]) & MASK].x;
	pt.y = urand[(y + (x << 2) + irand[(y + s) & MASK]) & MASK].y;
	return pt;
}

int ray_object(const struct object *obj, struct ray ray, struct spoint *sp) {
	switch(obj->type) {
	case OBJ_SPHERE:
		return ray_sphere(&obj->o.sph, ray, sp);

	case OBJ_PLANE:
		return ray_plane(&obj->o.pln, ray, sp);

	default:
		break;
	}
	return 0;
}

/* Calculate ray-sphere intersection, and return {1, 0} to signify hit or no hit.
 * Also the surface point parameters like position, normal, etc are returned through
 * the sp pointer if it is not NULL.
 */
int ray_sphere(const struct sphere *sph, struct ray ray, struct spoint *sp) {
	float a, b, c, d, sqrt_d, t1, t2;
	
	a = SQ(ray.dir.x) + SQ(ray.dir.y) + SQ(ray.dir.z);
	b = 2.0 * ray.dir.x * (ray.orig.x - sph->pos.x) +
				2.0 * ray.dir.y * (ray.orig.y - sph->pos.y) +
				2.0 * ray.dir.z * (ray.orig.z - sph->pos.z);
	c = SQ(sph->pos.x) + SQ(sph->pos.y) + SQ(sph->pos.z) +
				SQ(ray.orig.x) + SQ(ray.orig.y) + SQ(ray.orig.z) +
				2.0 * (-sph->pos.x * ray.orig.x - sph->pos.y * ray.orig.y - sph->pos.z * ray.orig.z) - SQ(sph->rad);
	
	if((d = SQ(b) - 4.0 * a * c) < 0.0) return 0;

	sqrt_d = sqrt(d);
	t1 = (-b + sqrt_d) / (2.0 * a);
	t2 = (-b - sqrt_d) / (2.0 * a);

	if((t1 < ERR_MARGIN && t2 < ERR_MARGIN) || (t1 > 1.0 && t2 > 1.0)) {
		return 0;
	}

	if(sp) {
		if(t1 < ERR_MARGIN) t1 = t2;
		if(t2 < ERR_MARGIN) t2 = t1;
		sp->dist = t1 < t2 ? t1 : t2;
		
		sp->pos.x = ray.orig.x + ray.dir.x * sp->dist;
		sp->pos.y = ray.orig.y + ray.dir.y * sp->dist;
		sp->pos.z = ray.orig.z + ray.dir.z * sp->dist;
		
		sp->normal.x = (sp->pos.x - sph->pos.x) / sph->rad;
		sp->normal.y = (sp->pos.y - sph->pos.y) / sph->rad;
		sp->normal.z = (sp->pos.z - sph->pos.z) / sph->rad;

		sp->vref = reflect(ray.dir, sp->normal);
		NORMALIZE(sp->vref);
	}
	return 1;
}

/* calculate ray-plane intersection... etc */
int ray_plane(const struct plane *pln, struct ray ray, struct spoint *sp) {
	struct vec3 to_orig;
	float t;

	if(fabs(DOT(pln->norm, ray.dir)) < ERR_MARGIN) {
		return 0;
	}

	to_orig.x = ray.orig.x - pln->pos.x;
	to_orig.y = ray.orig.y - pln->pos.y;
	to_orig.z = ray.orig.z - pln->pos.z;
	t = -DOT(pln->norm, to_orig) / DOT(pln->norm, ray.dir);

	if(t < ERR_MARGIN || t > 1.0) {
		return 0;
	}

	if(sp) {
		sp->dist = t;

		sp->pos.x = ray.orig.x + ray.dir.x * t;
		sp->pos.y = ray.orig.y + ray.dir.y * t;
		sp->pos.z = ray.orig.z + ray.dir.z * t;

		sp->normal.x = pln->norm.x;
		sp->normal.y = pln->norm.y;
		sp->normal.z = pln->norm.z;

		sp->vref = reflect(ray.dir, sp->normal);
		NORMALIZE(sp->vref);
	}
	return 1;
}

/* Load the scene from an extremely simple scene description file */
#define DELIM	" \t\n"
void load_scene(FILE *fp) {
	char line[256], *ptr, type;

	obj_list = malloc(sizeof *obj_list);
	obj_list->next = 0;
	
	while((ptr = fgets(line, 256, fp))) {
		int i;
		struct vec3 pos, norm, col;
		float rad, spow, refl;
		struct object *obj;
		
		while(*ptr == ' ' || *ptr == '\t') ptr++;
		if(*ptr == '#' || *ptr == '\n') continue;

		if(!(ptr = strtok(line, DELIM))) continue;
		type = *ptr;
		
		for(i=0; i<3; i++) {
			if(!(ptr = strtok(0, DELIM))) break;
			*((float*)&pos.x + i) = atof(ptr);
		}

		if(type == 'l') {
			lights[lnum++] = pos;
			continue;
		}

		if(!(ptr = strtok(0, DELIM))) continue;
		rad = atof(ptr);

		for(i=0; i<3; i++) {
			if(!(ptr = strtok(0, DELIM))) break;
			*((float*)&col.x + i) = atof(ptr);
		}

		if(type == 'c') {
			cam.pos = pos;
			cam.targ = col;
			cam.fov = rad;
			continue;
		}

		if(!(ptr = strtok(0, DELIM))) continue;
		spow = atof(ptr);

		if(!(ptr = strtok(0, DELIM))) continue;
		refl = atof(ptr);

		obj = malloc(sizeof *obj);

		obj->mat.col = col;
		obj->mat.spow = spow;
		obj->mat.refl = refl;

		if(type == 's') {
			obj->type = OBJ_SPHERE;
			obj->o.sph.pos = pos;
			obj->o.sph.rad = rad;

			obj->next = obj_list->next;
			obj_list->next = obj;
			continue;
		}

		for(i=0; i<3; i++) {
			if(!(ptr = strtok(0, DELIM))) break;
			*((float*)&norm.x + i) = atof(ptr);
		}
		NORMALIZE(norm);

		if(type == 'p') {
			obj->type = OBJ_PLANE;
			obj->o.pln.pos = pos;
			obj->o.pln.norm = norm;

			obj->next = obj_list->next;
			obj_list->next = obj;
		} else {
			fprintf(stderr, "unknown type: %c\n", type);
			free(obj);
		}
	}
}
