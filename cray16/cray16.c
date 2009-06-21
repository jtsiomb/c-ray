#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <float.h>
#include <alloc.h>
#include "gfx.h"
#include "color.h"
#include "log.h"

typedef unsigned long uint32_t;
typedef unsigned short uint16_t;
typedef unsigned char uint8_t;

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

void render(int xsz, int ysz, uint8_t far *fb);
void dither(void);
struct vec3 trace(struct ray ray, int depth);
struct vec3 shade(struct object *obj, struct spoint *sp, int depth);
struct vec3 reflect(struct vec3 v, struct vec3 n);
struct vec3 cross_product(struct vec3 v1, struct vec3 v2);
struct ray get_primary_ray(int x, int y);
struct vec3 get_sample_pos(int x, int y);
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
int xres = 320;
int yres = 240;
float aspect;
struct object *obj_list;
struct vec3 lights[MAX_LIGHTS];
int lnum = 0;
struct camera cam;

uint16_t far * far *scanline;


const char *usage = {
	"Usage: cray16 [options]\n"
	"Options:\n"
	"  -s WxH     where W is the width and H the height of the image\n"
	"  -i <file>  read from <file> instead of stdin\n"
	"  -h         this help screen\n\n"
	"Note: at the moment the only valid resolutions are: 320x240 and 320x200\n\n"
};


int main(int argc, char **argv) {
	int i;
	unsigned long rend_time, start_time;
	uint8_t far *pixels;
	FILE *infile = stdin;

	logfoo("starting c-ray16\n");

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

	if(!(scanline = farmalloc(yres * sizeof *scanline))) {
		perror("malloc failed");
		return 1;
	}
	for(i=0; i<yres; i++) {
		if(!(scanline[i] = farmalloc(xres * sizeof *scanline[i]))) {
			perror("malloc failed");
			return 1;
		}
	}

	/*
	if(!(pixels = malloc(xres * yres * sizeof *pixels))) {
		perror("pixel buffer allocation failed");
		return 1;
	}
	*/

	load_scene(infile);
	if(infile != stdin) fclose(infile);

	logfoo("initializing graphics: %dx%d 8bpp\n", xres, yres);
	if(set_video_mode(xres, yres, 8) == -1) {
		return 1;
	}

	init_palette();

	logfoo("starting rendering\n");
	render(xres, yres, vidmem);
	logfoo("rendering done\n");

	for(;;) {
		char c = getch();

		if(c == 'q' || c == 27) break;
		if(c == 'd') dither();
	}

	logfoo("restoring text mode\n");
	restore_vga();
	return 0;
}

#define SCAN_SKIP	8
void render(int xsz, int ysz, uint8_t far *fb) {
	int i, j, x, y;
	uint8_t far *ptr;

	int pitch = xsz;
	int xinc = 1;
	int planes = 1;

	if(xsz == 320 && ysz == 240) {	/* mode x */
		pitch >>= 2;
		xinc = 4;
		planes = 4;
	}

	for(i=0; i<SCAN_SKIP; i++) {
		for(j=0; j<planes; j++) {
			if(planes > 1) {
				set_plane_mask(1 << j);
			}

			for(y=i; y<ysz; y+=SCAN_SKIP) {
				ptr = fb + y * pitch;

				for(x=j; x<xsz; x+=xinc) {
					unsigned int cr, cg, cb;
					struct vec3 col = trace(get_primary_ray(x, y), 0);

					cr = (unsigned int)(col.x * 255.0);
					cg = (unsigned int)(col.y * 255.0);
					cb = (unsigned int)(col.z * 255.0);

					if(cr > 255) cr = 255;
					if(cg > 255) cg = 255;
					if(cb > 255) cb = 255;

					scanline[y][x] = ((cr << 8) & 0xf800) | ((cg << 3) & 0x7e0) |
						((cb >> 3) & 0x1f);
					/**ptr++ = alloc_color(cr, cg, cb);*/
					/**ptr++ = (cr & 0xe0) | ((cg >> 3) & 0x1c) | ((cb >> 6) & 3);*/
					*ptr++ = (((cr >> 5) & 7) << 5) |
						(((cg >> 5) & 7) << 2) |
						((cb >> 6) & 3);
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

#define PACK_RGB16(r, g, b) \
	((((r) << 8) & 0xf800) | (((g) << 3) & 0x7e0) | (((b) >> 3) & 0x1f))

#define UNPACK_RED16(x)		(((x) & 0xf800) >> 8)
#define UNPACK_GREEN16(x)	(((x) & 0x7e0) >> 3)
#define UNPACK_BLUE16(x)	(((x) & 0x1f) << 3)

void add_to_pixel(uint16_t far *pix, int err_r, int err_g, int err_b, float fact);

void dither(void)
{
	int x, y;
	unsigned char far *fb = vidmem;

	int pitch = xres;
	int planes = 1;

	if(xres == 320 && yres == 240) {	/* mode x */
		pitch >>= 2;
		planes = 4;
	}


	for(y=0; y<yres; y++) {
		for(x=0; x<xres; x++) {
			int ref_r, ref_g, ref_b, r, g, b, idx;
			int err_r, err_g, err_b, p;

			ref_r = UNPACK_RED16(scanline[y][x]);
			ref_g = UNPACK_GREEN16(scanline[y][x]);
			ref_b = UNPACK_BLUE16(scanline[y][x]);

			idx = alloc_color(ref_r, ref_g, ref_b);
			lookup_color(idx, &r, &g, &b);

			err_r = ref_r - r;
			err_g = ref_g - g;
			err_b = ref_b - b;

			if(x < xres - 1) {
				add_to_pixel(&scanline[y][x + 1], err_r, err_g, err_b, 0.4375);
			}
			if(y < yres - 1) {
				if(x > 0) {
					add_to_pixel(&scanline[y + 1][x - 1], err_r, err_g, err_b, 0.1875);
				}
				add_to_pixel(&scanline[y + 1][x], err_r, err_g, err_b, 0.3125);
				if(x < xres - 1) {
					add_to_pixel(&scanline[y + 1][x + 1], err_r, err_g, err_b, 0.0625);
				}
			}

			if(planes > 1) {
				p = x & 3;
				set_plane_mask(1 << p);

				*fb = idx;
				if(p == 3) {
					fb++;
				}
			} else {
				*fb++ = idx;
			}
		}
	}
}

void add_to_pixel(uint16_t far *pix, int err_r, int err_g, int err_b, float fact)
{
	int r = UNPACK_RED16(*pix);
	int g = UNPACK_GREEN16(*pix);
	int b = UNPACK_BLUE16(*pix);

	r += (int)(fact * err_r);
	g += (int)(fact * err_g);
	b += (int)(fact * err_b);

	if(r > 255) r = 255;
	if(g > 255) g = 255;
	if(b > 255) b = 255;

	*pix = PACK_RGB16(r, g, b);
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
struct ray get_primary_ray(int x, int y) {
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
	
	ray.dir = get_sample_pos(x, y);
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


struct vec3 get_sample_pos(int x, int y) {
	struct vec3 pt;
	static float sf = 0.0;

	if(sf == 0.0) {
		sf = 1.5 / (float)xres;
	}

	pt.x = ((float)x / (float)xres - 0.5) * 2.0 * aspect;
	pt.y = (0.5 - (float)y / (float)yres) * 2.0;
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

	logfoo("loading scene file\n");

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
