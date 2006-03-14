#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include "libgba.h"

/*int __gba_multiboot;*/

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
	struct material mat;
	struct sphere *next;
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
struct vec3 shade(struct sphere *obj, struct spoint *sp, int depth);
struct vec3 reflect(struct vec3 v, struct vec3 n);
struct vec3 cross_product(struct vec3 v1, struct vec3 v2);
struct ray get_primary_ray(int x, int y, int sample);
struct vec3 get_sample_pos(int x, int y, int sample);
struct vec3 jitter(int x, int y, int s);
int ray_sphere(const struct sphere *sph, struct ray ray, struct spoint *sp);
void create_scene(void);
unsigned long get_msec(void);

#define MAX_LIGHTS		16				/* maximum number of lights */
#define RAY_MAG			1000.0			/* trace rays of this magnitude */
#define MAX_RAY_DEPTH	5				/* raytrace recursion limit */
#define FOV				0.78539816		/* field of view in rads (pi/4) */
#define HALF_FOV		(FOV * 0.5)
#define ERR_MARGIN		1e-5			/* an arbitrary error margin to avoid surface acne */

/* bit-shift ammount for packing each color into a 32bit uint */
#ifdef LITTLE_ENDIAN
#define RSHIFT	16
#define BSHIFT	0
#else	/* big endian */
#define RSHIFT	0
#define BSHIFT	16
#endif	/* endianess */
#define GSHIFT	8	/* this is the same in both byte orders */

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
int xres = 800;
int yres = 600;
float aspect = 1.333333;
struct sphere *obj_list;
struct vec3 lights[MAX_LIGHTS];
int lnum = 0;
struct camera cam;
unsigned int rays = 0, prim_rays = 0;


int main(int argc, char **argv) {
	unsigned long rend_time, start_time;
	static char str[256];

	xres = 240;
	yres = 160;
	aspect = (float)xres / (float)yres;

	gba_init();
	set_video_mode(VMODE_LFB_240x160_16, 1);
	clear_buffer(front_buffer, RGB(160, 200, 255));

	draw_string("rendering...", 0, 8, front_buffer);
	create_scene();

	start_time = get_millisec();
	render(xres, yres, front_buffer->pixels, 1);
	rend_time = get_millisec() - start_time;

	memcpy(back_buffer->pixels, front_buffer->pixels, xres * yres * sizeof(uint16_t));

	{	/* stats */
		int line = 0;
		sprintf(str, "rendering time: %lu sec", rend_time / 1000);
		draw_string(str, 0, line++ * 16, front_buffer);
		sprintf(str, "    (%lu msec)", rend_time);
		draw_string(str, 0, line++ * 16, front_buffer);
		line++;

		sprintf(str, "total rays cast: %u", rays);
		draw_string(str, 0, line++ * 16, front_buffer);
		sprintf(str, "   primary rays: %u", prim_rays);
		draw_string(str, 0, line++ * 16, front_buffer);
		line++;

		draw_string("Press any key to clear text", 0, line++ * 16, front_buffer);
		getchar();
	}
	flip();
	
	clr_int();
	halt();
	return 0;
}

#define START_PSIZE		16

/* render a frame of xsz/ysz dimensions into the provided framebuffer */
void render(int xsz, int ysz, uint32_t *fb, int samples) {
	int i, j, k, x, y;
	uint16_t *fb16;
	char *wbuf;
	char *tmp = malloc(38400);

	for(k=START_PSIZE; k>0; k>>=1) {
		int scanoffs = (k - 1) * xsz;
		
		fb16 = (uint16_t*)fb;
		wbuf = tmp;
		for(j=0; j<ysz; j+=k) {
			for(i=0; i<xsz; i+=k) {
				uint16_t pixel;
				int calced = *wbuf == 1;

				if(!calced) {
					struct vec3 col = trace(get_primary_ray(i, j, 0), 0);
					prim_rays++;
				
					pixel = RGB(MIN((uint16_t)(col.x * 255.0), 255),
								MIN((uint16_t)(col.y * 255.0), 255),
								MIN((uint16_t)(col.z * 255.0), 255));
					*wbuf = 1;
				} else {
					pixel = *fb16;
				}
				
				uint16_t *sub = fb16;
				for(y=0; y<k; y++) {
					for(x=0; x<k; x++) {
						if(x || y || !calced) {
							*sub = pixel;
						}
						sub++;
					}
					sub += xsz - k;
				}

				fb16 += k;
				wbuf += k;
			}
			fb16 += scanoffs;
			wbuf += scanoffs;
		}
	}

	free(tmp);
}

/* trace a ray throught the scene recursively (the recursion happens through
 * shade() to calculate reflection rays if necessary).
 */
struct vec3 trace(struct ray ray, int depth) {
	struct vec3 col;
	struct spoint sp, nearest_sp;
	struct sphere *nearest_obj = 0;
	struct sphere *iter = obj_list->next;

	rays++;

	/* if we reached the recursion limit, bail out */
	if(depth >= MAX_RAY_DEPTH) {
		col.x = col.y = col.z = 0.0;
		return col;
	}
	
	/* find the nearest intersection ... */
	while(iter) {
		if(ray_sphere(iter, ray, &sp)) {
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
struct vec3 shade(struct sphere *obj, struct spoint *sp, int depth) {
	int i;
	struct vec3 col = {0, 0, 0};

	/* for all lights ... */
	for(i=0; i<lnum; i++) {
		float ispec, idiff;
		struct vec3 ldir;
		struct ray shadow_ray;
		struct sphere *iter = obj_list->next;
		int in_shadow = 0;

		ldir.x = lights[i].x - sp->pos.x;
		ldir.y = lights[i].y - sp->pos.y;
		ldir.z = lights[i].z - sp->pos.z;

		shadow_ray.orig = sp->pos;
		shadow_ray.dir = ldir;

		/* shoot shadow rays to determine if we have a line of sight with the light */
		while(iter) {
			if(ray_sphere(iter, shadow_ray, 0)) {
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
	struct vec3 i, j = {0, 1, 0}, k, dir, orig, foo;

	k.x = cam.targ.x - cam.pos.x;
	k.y = cam.targ.y - cam.pos.y;
	k.z = cam.targ.z - cam.pos.z;
	NORMALIZE(k);

	i = cross_product(j, k);
	j = cross_product(k, i);
	m[0][0] = i.x; m[0][1] = j.x; m[0][2] = k.x;
	m[1][0] = i.y; m[1][1] = j.y; m[1][2] = k.y;
	m[2][0] = i.z; m[2][1] = j.z; m[2][2] = k.z;
	
	ray.orig.x = ray.orig.y = ray.orig.z = 0.0;
	ray.dir = get_sample_pos(x, y, sample);
	ray.dir.z = 1.0 / HALF_FOV;
	ray.dir.x *= RAY_MAG;
	ray.dir.y *= RAY_MAG;
	ray.dir.z *= RAY_MAG;
	
	dir.x = ray.dir.x + ray.orig.x;
	dir.y = ray.dir.y + ray.orig.y;
	dir.z = ray.dir.z + ray.orig.z;
	foo.x = dir.x * m[0][0] + dir.y * m[0][1] + dir.z * m[0][2];
	foo.y = dir.x * m[1][0] + dir.y * m[1][1] + dir.z * m[1][2];
	foo.z = dir.x * m[2][0] + dir.y * m[2][1] + dir.z * m[2][2];

	orig.x = ray.orig.x * m[0][0] + ray.orig.y * m[0][1] + ray.orig.z * m[0][2] + cam.pos.x;
	orig.y = ray.orig.x * m[1][0] + ray.orig.y * m[1][1] + ray.orig.z * m[1][2] + cam.pos.y;
	orig.z = ray.orig.x * m[2][0] + ray.orig.y * m[2][1] + ray.orig.z * m[2][2] + cam.pos.z;

	ray.orig = orig;
	ray.dir.x = foo.x + orig.x;
	ray.dir.y = foo.y + orig.y;
	ray.dir.z = foo.z + orig.z;
	
	return ray;
}


struct vec3 get_sample_pos(int x, int y, int sample) {
	struct vec3 pt;
	/*float xsz = 2.0, ysz = xres / aspect;*/
	static float sf = 0.0;

	if(sf == 0.0) {
		sf = 2.0 / (float)xres;
	}

	pt.x = ((float)x / (float)xres) - 0.5;
	pt.y = -(((float)y / (float)yres) - 0.65) / aspect;

	/*
	if(sample) {
		struct vec3 jt = jitter(x, y, sample);
		pt.x += jt.x * sf;
		pt.y += jt.y * sf / aspect;
	}
	*/
	return pt;
}

/* jitter function taken from Graphics Gems I. */
/*
struct vec3 jitter(int x, int y, int s) {
	struct vec3 pt;
	pt.x = urand[(x + (y << 2) + irand[(x + s) & MASK]) & MASK].x;
	pt.y = urand[(y + (x << 2) + irand[(y + s) & MASK]) & MASK].y;
	return pt;
}
*/

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

	if((t1 < ERR_MARGIN && t2 < ERR_MARGIN) || (t1 > 1.0 && t2 > 1.0)) return 0;

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

void create_scene(void) {
	struct sphere *sph;

	/* add lights */
	lnum = 2;
	lights[0] = (struct vec3){-50.0, 100.0, -50.0};
	lights[1] = (struct vec3){40.0, 40.0, 150.0};

	/* add camera */
	cam.pos = (struct vec3){0.0, 6.0, -17.0};
	cam.targ = (struct vec3){0.0, -1.0, 0.0};
	cam.fov = 45.0;

	/* add spheres */
	obj_list = malloc(sizeof(struct sphere));
	obj_list->next = 0;
	
	sph = malloc(sizeof *sph);
	sph->next = obj_list->next;
	obj_list->next = sph;
	sph->pos = (struct vec3){-1.5, -0.3, -1.0};
	sph->rad = 0.7;
	sph->mat.col = (struct vec3){1.0, 0.2, 0.05};
	sph->mat.spow = 50.0;
	sph->mat.refl = 0.3;

	sph = malloc(sizeof *sph);
	sph->next = obj_list->next;
	obj_list->next = sph;
	sph->pos = (struct vec3){1.5, -0.4, 0.0};
	sph->rad = 0.6;
	sph->mat.col = (struct vec3){0.1, 0.85, 1.0};
	sph->mat.spow = 50.0;
	sph->mat.refl = 0.4;

	sph = malloc(sizeof *sph);
	sph->next = obj_list->next;
	obj_list->next = sph;
	sph->pos = (struct vec3){0.0, -1000.0, 2.0};
	sph->rad = 999.0;
	sph->mat.col = (struct vec3){0.1, 0.2, 0.6};
	sph->mat.spow = 80.0;
	sph->mat.refl = 0.8;

	sph = malloc(sizeof *sph);
	sph->next = obj_list->next;
	obj_list->next = sph;
	sph->pos = (struct vec3){0.0, 0.0, 2.0};
	sph->rad = 1.0;
	sph->mat.col = (struct vec3){1.0, 0.5, 0.1};
	sph->mat.spow = 60.0;
	sph->mat.refl = 0.7;
}
