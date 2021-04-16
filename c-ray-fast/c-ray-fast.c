/* c-ray-fast - a simple raytracing filter based on c-ray-f and c-ray-mt.
 * Copyright (C) 2006-2021 John Tsiombikas <nuclear@member.fsf.org>
 *
 * You are free to use, modify and redistribute this program under the terms
 * of the GNU General Public License v3 or later. See COPYING for details.
 * ---------------------------------------------------------------------
 * Usage:
 *   compile:  cc -o c-ray-fast c-ray-fast.c -lm -lpthread
 *   run:      cat scene | ./c-ray-fast >foo.ppm
 *   enjoy:    display foo.ppm (with imagemagick)
 *      or:    imgview foo.ppm (on IRIX)
 * ---------------------------------------------------------------------
 * Scene file format:
 *   # options (samples/gamma/ambient/bgcolor)
 *   ? option  value
 *   # material
 *   M "name"  dr dg db  sr sg sb  shininess  reflect  transmit  ior
 *   # sphere
 *   s  x y z  rad   "material name"
 *   p  x y z  nx ny nz  "material name"
 *   # light
 *   l  x y z  r g b
 *   # camera
 *   c  x y z  fov   tx ty tz
 * ---------------------------------------------------------------------
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <errno.h>
#include <pthread.h>

#define VER_STR		"c-ray-fast v1.0"
#define COMMENT		"# rendered with " VER_STR

#if defined(__APPLE__) && defined(__MACH__)
#ifndef __unix__
#define __unix__	1
#endif	/* unix */
#ifndef __bsd__
#define __bsd__	1
#endif	/* bsd */
#endif	/* apple */

/* find the appropriate way to define explicitly sized types */
/* for C99 or GNU libc (also mach's libc) we can use stdint.h */
#if (__STDC_VERSION__ >= 199900) || defined(__GLIBC__) || defined(__MACH__)
#include <stdint.h>
#elif defined(unix) || defined(__unix__)	/* some UNIX systems have them in sys/types.h */
#include <sys/types.h>
#elif defined(__WIN32__) || defined(WIN32)	/* the nameless one */
typedef unsigned __int8 uint8_t;
typedef unsigned __int32 uint32_t;
#endif	/* sized type detection */

struct vec3 {
	double x, y, z;
};

struct color {
	float r, g, b, a;
};

struct ray {
	struct vec3 orig, dir, invdir;
};

struct material {
	char *name;
	struct color kd, ks; /* color (diffuse/specular) */
	double spow;		/* specular power */
	double refl;		/* reflection intensity */
	double refr;		/* refraction (transmission) intensity */
	double ior;			/* index of refraction */
};
struct material def_mtl = {0, {1, 1, 1}, {1, 1, 1}, 50, 0, 0, 0};

struct sphere {
	struct vec3 pos;
	double rad;
};

struct plane {
	struct vec3 pos;
	struct vec3 norm;
};

enum {
	OBJ_SPHERE = 1,
	OBJ_PLANE = 2
};

struct object {
	int type;
	union {
		struct sphere sph;
		struct plane pln;
	} o;
	struct material mat;
};

struct spoint {
	struct vec3 pos, normal, vref;	/* position, normal and view reflection */
	double dist;		/* parametric distance of intersection along the ray */
	struct object *obj;
};

struct camera {
	struct vec3 pos, targ;
	double fov;
};

struct light {
	struct vec3 pos;
	struct color col;
};


struct aabox {
	struct vec3 vmin, vmax;
};

struct octnode {
	struct aabox aabb;
	struct object **obj;
	int num_obj, max_obj;
	struct octnode *sub[8];
};

struct block {
	struct color *start;
	int x, y, width, height;
};

void render(void);
void *thread_func(void *arg);
void trace(struct ray *ray, int depth, struct color *col);
void shade(struct ray *ray, struct spoint *sp, int depth, struct color *col);
int ray_scene(struct ray *ray, struct spoint *sp);
void reflect(struct vec3 *res, struct vec3 *v, struct vec3 *n);
void refract(struct vec3 *res, struct vec3 *v, struct vec3 *n, double from_ior, double to_ior);
void cross_product(struct vec3 *res, struct vec3 *v1, struct vec3 *v2);
void get_primary_ray(struct ray *ray, int x, int y, int sample);
void get_sample_pos(struct vec3 *res, int x, int y, int sample);
void jitter(struct vec3 *res, int x, int y, int s);
int ray_object(struct object *obj, struct ray *ray, struct spoint *sp);
int ray_sphere(struct sphere *sph, struct ray *ray, struct spoint *sp);
int ray_plane(struct plane *pln, struct ray *ray, struct spoint *sp);
int ray_aabox(struct aabox *aabox, struct ray *ray);
int box_sphere(struct aabox *box, struct sphere *sph);
int load_scene(FILE *fp);
int build_octree(struct octnode *node, int depth);
int populate_octnode(struct octnode *node, struct object **obj, int num_obj);
int isect_octnode(struct octnode *node, struct ray *ray, struct spoint *spret);
void print_obj(void);
void print_progr(int x, int max);
unsigned long get_msec(void);
int detect_nproc(void);

#define RAY_MAG			1000.0			/* trace rays of this magnitude */
#define MAX_RAY_DEPTH	5				/* raytrace recursion limit */
#define ERR_MARGIN		1e-6			/* an arbitrary error margin to avoid surface acne */

#define OCT_MAX_OBJ		64
#define OCT_MAX_DEPTH	9

/* some helpful macros... */
#define SQ(x)		((x) * (x))
#define MAX(a, b)	((a) > (b) ? (a) : (b))
#define MIN(a, b)	((a) < (b) ? (a) : (b))
#define CLAMP(x, a, b)	MIN(MAX(x, a), b)
#define DOT(a, b)	((a).x * (b).x + (a).y * (b).y + (a).z * (b).z)
#define NORMALIZE(a)  do {\
	double len = sqrt(DOT(a, a));\
	(a).x /= len; (a).y /= len; (a).z /= len;\
} while(0);
#define LERP(a, b, t)	((a) + ((b) - (a)) * (t))

/* global state */
int xres = 800;
int yres = 600;
double aspect = 1.333333;
struct color *pixels;

struct object *objects;
struct object *infobj;	/* infinite objects (planes) */
struct material *materials;
struct light *lights;
int num_obj, num_inf, num_mtl, num_lt;
struct camera cam;
double inv_gamma = 1.0 / 2.2;
int verbose, progress;
struct color bgcol;
struct color ambcol = {0.01, 0.01, 0.01, 0};
int rays_per_pixel = 1;
double rcp_samples;

enum { OPT_GAMMA = 1, OPT_SAMPLES = 2, OPT_BGCOL = 4 };
unsigned int override = 0;

struct octnode *octroot;

int blk_width, blk_height;
struct block *blocks;
int num_blk, cur_blk;
pthread_mutex_t job_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t job_cond = PTHREAD_COND_INITIALIZER;
pthread_t *threads;
int num_thr;

#define NRAN	1024
#define JITTER_MASK	(NRAN - 1)
struct vec3 urand[NRAN];
int irand[NRAN];

const char *usage = {
	"Usage: c-ray-f [options]\n"
	"  Reads a scene file from stdin, writes the image to stdout, and stats to stderr.\n\n"
	"Options:\n"
	"  -s WxH     where W is the width and H the height of the image\n"
	"  -r <rays>  shoot <rays> rays per pixel (antialiasing)\n"
	"  -i <file>  read from <file> instead of stdin\n"
	"  -o <file>  write to <file> instead of stdout\n"
	"  -a <file>  write alpha mask to <file>\n"
	"  -g <gamma> set gamma correction value (default: 2.2)\n"
	"  -p         print progress to stderr while rendering\n"
	"  -b <color> set background color (example: #ff0000 or 255,0,0 w/o spaces)\n"
	"  -v         verbose output (to stderr)\n"
	"  -h         this help screen\n\n"
};

int main(int argc, char **argv)
{
	int i, r, g, b, a;
	unsigned long rend_time, start_time;
	FILE *infile = stdin, *outfile = stdout, *outfile_alpha = 0;
	char *endp;
	float gamma;
	struct object *obj;
	struct aabox box;

	for(i=1; i<argc; i++) {
		if(argv[i][0] == '-') {
			if(argv[i][2] != 0) {
invalopt:		fprintf(stderr, "invalid option: %s\n", argv[i]);
				fputs(usage, stderr);
				return 1;
			}

			switch(argv[i][1]) {
			case 's':
				if(!argv[++i] || sscanf(argv[i], "%dx%d", &xres, &yres) != 2) {
					fputs("-s must be followed by something like \"640x480\"\n", stderr);
					return 1;
				}
				aspect = (double)xres / (double)yres;
				break;

			case 'i':
				if(!argv[++i] || !(infile = fopen(argv[i], "rb"))) {
					fprintf(stderr, "failed to open input file %s: %s\n", argv[i], strerror(errno));
					return 1;
				}
				break;

			case 'o':
				if(!argv[++i] || !(outfile = fopen(argv[i], "wb"))) {
					fprintf(stderr, "failed to open output file %s: %s\n", argv[i], strerror(errno));
					return 1;
				}
				break;

			case 'r':
				if(!argv[++i] || (rays_per_pixel = atoi(argv[i])) <= 0) {
					fputs("-r must be followed by a number (rays per pixel)\n", stderr);
					return 1;
				}
				override |= OPT_SAMPLES;
				break;

			case 'v':
				verbose = 1;
				break;

			case 'g':
				if(!argv[++i] || (gamma = atof(argv[i])) == 0.0) {
					fputs("-g must be followed by a (non-zero) gamma value\n", stderr);
					return 1;
				}
				inv_gamma = 1.0 / gamma;
				override |= OPT_GAMMA;
				break;

			case 'p':
				progress = 1;
				break;

			case 'a':
				if(!argv[++i] || !(outfile_alpha = fopen(argv[i], "wb"))) {
					fprintf(stderr, "failed to open alpha output file %s: %s\n", argv[i], strerror(errno));
					return 1;
				}
				break;

			case 'b':
				if(!argv[++i]) {
					fprintf(stderr, "-b must be followed by a background color\n");
					return 1;
				}
				if(argv[i][0] == '#') {
					long col = strtol(argv[i] + 1, &endp, 16);
					if(endp == argv[i] + 1) {
						fprintf(stderr, "-b followed by invalid HTML color format\n");
						return 1;
					}
					bgcol.r = ((col >> 16) & 0xff) / 255.0f;
					bgcol.g = ((col >> 8) & 0xff) / 255.0f;
					bgcol.b = (col & 0xff) / 255.0f;
				} else {
					int r, g, b;
					if(sscanf(argv[i], "%d,%d,%d", &r, &g, &b) != 3 || r < 0 || r >= 256 ||
							g < 0 || g >= 256 || b < 0 || b >= 256) {
						fprintf(stderr, "-b followed by invalid R,G,B color value\n");
						return 1;
					}
					bgcol.r = r / 255.0f;
					bgcol.g = g / 255.0f;
					bgcol.b = b / 255.0f;
				}
				bgcol.a = 0.0f;
				override |= OPT_BGCOL;
				break;

			case 'B':
				if(!argv[++i]) {
					fprintf(stderr, "-B must be followed by the rendering block size (WxH)\n");
					return 1;
				}
				if(sscanf(argv[i], "%dx%d", &blk_width, &blk_height) < 1 || blk_width <= 0 || blk_height < 0) {
					fprintf(stderr, "-B followed by invalid block size: %s\n", argv[i]);
					return 1;
				}
				if(!blk_height) blk_height = blk_width;
				break;

			case 't':
				if(!argv[++i]) {
					fprintf(stderr, "-t must be followed by the number of threads\n");
					return 1;
				}
				num_thr = atoi(argv[i]);
				break;

			case 'h':
				fputs(usage, stdout);
				return 0;

			default:
				goto invalopt;
			}
		} else {
			if(!(infile = fopen(argv[i], "rb"))) {
				fprintf(stderr, "failed to open input file %s: %s\n", argv[i], strerror(errno));
				return 1;
			}
		}
	}

	if(num_thr <= 0) {
		num_thr = detect_nproc();
	}
	if(!(threads = malloc(num_thr * sizeof *threads))) {
		perror("failed to allocate threads");
		return 1;
	}

	for(i=0; i<num_thr; i++) {
		if(pthread_create(threads + i, 0, thread_func, 0) == -1) {
			fprintf(stderr, "failed to spawn thread %d\n", i);
			return 1;
		}
	}

	if(verbose) {
		fputs(VER_STR "\n", stderr);
		fprintf(stderr, "Using %d threads\n", num_thr);
	}

	/* initialize the random number tables for the jitter */
	for(i=0; i<NRAN; i++) urand[i].x = (double)rand() / RAND_MAX - 0.5;
	for(i=0; i<NRAN; i++) urand[i].y = (double)rand() / RAND_MAX - 0.5;
	for(i=0; i<NRAN; i++) irand[i] = (int)(NRAN * ((double)rand() / RAND_MAX));

	if(!(pixels = malloc(xres * yres * sizeof *pixels))) {
		perror("pixel buffer allocation failed");
		return 1;
	}
	if(load_scene(infile) == -1) {
		return 1;
	}

	if(!blk_width) {
		blk_width = blk_height = 32;
	}
	if(blk_width > xres) blk_width = xres;
	if(blk_height > yres) blk_height = yres;

	/* prepare acceleration structure */
	if(!(octroot = malloc(sizeof *octroot)) || !(octroot->obj = malloc(num_obj * sizeof *octroot->obj))) {
		perror("failed to allocate octree root node");
		return 1;
	}
	octroot->num_obj = 0;

	box.vmin.x = box.vmin.y = box.vmin.z = DBL_MAX;
	box.vmax.x = box.vmax.y = box.vmax.z = -DBL_MAX;

	obj = objects;
	for(i=0; i<num_obj; i++) {
		if(obj->type == OBJ_SPHERE) {
			struct sphere *sph = &obj->o.sph;
			double d;

			d = sph->pos.x - sph->rad;
			if(d < box.vmin.x) box.vmin.x = d;
			d = sph->pos.x + sph->rad;
			if(d > box.vmax.x) box.vmax.x = d;
			d = sph->pos.y - sph->rad;
			if(d < box.vmin.y) box.vmin.y = d;
			d = sph->pos.y + sph->rad;
			if(d > box.vmax.y) box.vmax.y = d;
			d = sph->pos.z - sph->rad;
			if(d < box.vmin.z) box.vmin.z = d;
			d = sph->pos.z + sph->rad;
			if(d > box.vmax.z) box.vmax.z = d;

			octroot->obj[octroot->num_obj++] = obj;
		}
		obj++;
	}
	octroot->aabb = box;

	if(verbose) fprintf(stderr, "Building acceleration structure... ");
	start_time = get_msec();
	if(build_octree(octroot, 0) == -1) {
		return 1;
	}
	if(verbose) {
		rend_time = get_msec() - start_time;
		fprintf(stderr, "%lu seconds (%lu milliseconds)\n", rend_time / 1000, rend_time);
	}

	/* start rendering */
	start_time = get_msec();
	render();
	rend_time = get_msec() - start_time;

	/* output statistics to stderr */
	if(progress) fputc('\n', stderr);
	fprintf(stderr, "Rendering took: %lu seconds (%lu milliseconds)\n", rend_time / 1000, rend_time);

	/* output the image */
	fprintf(outfile, "P6\n" COMMENT "\n%d %d\n255\n", xres, yres);
	if(outfile_alpha) {
		fprintf(outfile_alpha, "P5\n" COMMENT "\n%d %d\n255\n", xres, yres);
	}
	for(i=0; i<xres * yres; i++) {
		r = pow(pixels->r, inv_gamma) * 255.0;
		g = pow(pixels->g, inv_gamma) * 255.0;
		b = pow(pixels->b, inv_gamma) * 255.0;

		fputc(CLAMP(r, 0, 255), outfile);
		fputc(CLAMP(g, 0, 255), outfile);
		fputc(CLAMP(b, 0, 255), outfile);

		if(outfile_alpha) {
			a = pixels->a * 255.0;
			fputc(CLAMP(a, 0, 255), outfile_alpha);
		}
		pixels++;
	}

	if(infile != stdin) fclose(infile);
	if(outfile != stdout) {
		fclose(outfile);
	} else {
		fflush(outfile);
	}
	if(outfile_alpha) fclose(outfile_alpha);
	return 0;
}

void render(void)
{
	int i, j, bx, by, bheight, nxblk, nyblk;
	struct color *fb;
	struct block *blk;

	nxblk = (xres + blk_width - 1) / blk_width;
	nyblk = (yres + blk_height - 1) / blk_height;
	if(!(blocks = malloc(nxblk * nyblk * sizeof *blocks))) {
		perror("failed to allocate render blocks");
		return;
	}
	blk = blocks;

	fb = pixels;
	rcp_samples = 1.0 / (double)rays_per_pixel;

	pthread_mutex_lock(&job_mutex);
	for(j=0; j<nyblk; j++) {
		by = j * blk_height;
		bheight = yres - by < blk_height ? yres - by : blk_height;

		for(i=0; i<nxblk; i++) {
			bx = i * blk_width;

			blk->start = fb;
			blk->x = bx;
			blk->y = by;
			blk->width = xres - bx < blk_width ? xres - bx : blk_width;
			blk->height = bheight;

			fb += blk->width;
			blk++;
		}

		fb += xres * (blk_height - 1);
	}

	num_blk = nxblk * nyblk;
	cur_blk = 0;
	pthread_mutex_unlock(&job_mutex);
	pthread_cond_broadcast(&job_cond);

	for(i=0; i<num_thr; i++) {
		pthread_join(threads[i], 0);
	}
}

void *thread_func(void *arg)
{
	int i, j, x, y, s;
	float r, g, b, a;
	struct ray ray;
	struct color col, *pix;
	struct block *blk;

	pthread_mutex_lock(&job_mutex);
	while(!num_blk) {
		pthread_cond_wait(&job_cond, &job_mutex);
	}

	while(cur_blk < num_blk) {
		blk = blocks + cur_blk++;
		pthread_mutex_unlock(&job_mutex);

		pix = blk->start;
		y = blk->y;
		for(i=0; i<blk->height; i++) {
			x = blk->x;
			for(j=0; j<blk->width; j++) {
				r = g = b = a = 0.0;

				for(s=0; s<rays_per_pixel; s++) {
					get_primary_ray(&ray, x, y, s);
					trace(&ray, 0, &col);
					r += col.r;
					g += col.g;
					b += col.b;
					a += col.a;
				}

				pix->r = r * rcp_samples;
				pix->g = g * rcp_samples;
				pix->b = b * rcp_samples;
				pix->a = a * rcp_samples;
				pix++;
				x++;
			}
			pix += xres - blk->width;
			y++;
		}

		pthread_mutex_lock(&job_mutex);
		if(progress) {
			print_progr(cur_blk, num_blk);
		}
	}
	pthread_mutex_unlock(&job_mutex);
	return 0;
}

#define USE_OCTREE

int ray_scene(struct ray *ray, struct spoint *retsp)
{
	int i;
	struct spoint sp, sp0;

#ifdef USE_OCTREE
	ray->invdir.x = 1.0 / ray->dir.x;
	ray->invdir.y = 1.0 / ray->dir.y;
	ray->invdir.z = 1.0 / ray->dir.z;

	if(!retsp) {
		for(i=0; i<num_inf; i++) {
			if(ray_object(infobj + i, ray, 0)) {
				return 1;
			}
		}
		return isect_octnode(octroot, ray, 0);
	}

	sp0.dist = DBL_MAX;
	sp0.obj = 0;

	for(i=0; i<num_inf; i++) {
		if(ray_object(infobj + i, ray, &sp) && sp.dist < sp0.dist) {
			sp0 = sp;
		}
	}
	if(isect_octnode(octroot, ray, &sp) && sp.dist < sp0.dist) {
		sp0 = sp;
	}

#else
	if(!retsp) {
		for(i=0; i<num_obj; i++) {
			if(ray_object(objects + i, ray, 0)) {
				return 1;
			}
		}
		return 0;
	}

	sp0.dist = DBL_MAX;
	sp0.obj = 0;

	/* find the nearest intersection ... */
	for(i=0; i<num_obj; i++) {
		if(ray_object(objects + i, ray, &sp)) {
			if(sp.dist < sp0.dist) {
				sp0 = sp;
			}
		}
	}
#endif

	if(sp0.obj) {
		if(retsp) *retsp = sp0;
		return 1;
	}
	return 0;
}

/* trace a ray throught the scene recursively (the recursion happens through
 * shade() to calculate reflection rays if necessary).
 */
void trace(struct ray *ray, int depth, struct color *col)
{
	struct spoint sp;

	/* if we reached the recursion limit, bail out */
	if(depth >= MAX_RAY_DEPTH) {
		*col = bgcol;
		return;
	}

	if(ray_scene(ray, &sp)) {
		shade(ray, &sp, depth, col);
		col->a = 1.0f;
	} else {
		*col = bgcol;
	}
}

/* Calculates direct illumination with the phong reflectance model.
 * Also handles reflections by calling trace again, if necessary.
 */
void shade(struct ray *ray, struct spoint *sp, int depth, struct color *col)
{
	int i, entering;
	float ispec, idiff;
	struct vec3 ldir, rdir;
	struct ray rray, shadow_ray;
	struct color rcol;
	double from_ior, to_ior;
	struct material *mtl = &sp->obj->mat;

	if(DOT(ray->dir, sp->normal) > 0.0) {
		sp->normal.x = -sp->normal.x;
		sp->normal.y = -sp->normal.y;
		sp->normal.z = -sp->normal.z;
		entering = 0;
	} else {
		entering = 1;
	}

	*col = ambcol;

	/* for all lights ... */
	for(i=0; i<num_lt; i++) {
		ldir.x = lights[i].pos.x - sp->pos.x;
		ldir.y = lights[i].pos.y - sp->pos.y;
		ldir.z = lights[i].pos.z - sp->pos.z;

		shadow_ray.orig = sp->pos;
		shadow_ray.dir = ldir;

		/* shoot shadow rays to determine if we have a line of sight with the light */
		if(!ray_scene(&shadow_ray, 0)) {
			NORMALIZE(ldir);

			idiff = MAX(DOT(sp->normal, ldir), 0.0);
			ispec = mtl->spow > 0.0 ? pow(MAX(DOT(sp->vref, ldir), 0.0), mtl->spow) : 0.0;

			col->r += (idiff * mtl->kd.r + ispec * mtl->ks.r) * lights[i].col.r;
			col->g += (idiff * mtl->kd.g + ispec * mtl->ks.g) * lights[i].col.g;
			col->b += (idiff * mtl->kd.b + ispec * mtl->ks.b) * lights[i].col.b;
		}
	}

	/* Also, if the object is reflective, spawn a reflection ray, and call trace()
	 * to calculate the light arriving from the mirror direction.
	 */
	if(mtl->refl > 0.0) {
		rray.orig = sp->pos;
		rray.dir = sp->vref;
		rray.dir.x *= RAY_MAG;
		rray.dir.y *= RAY_MAG;
		rray.dir.z *= RAY_MAG;

		trace(&rray, depth + 1, &rcol);
		col->r += rcol.r * mtl->refl;
		col->g += rcol.g * mtl->refl;
		col->b += rcol.b * mtl->refl;
	}

	if(mtl->refr > 0.0) {
		from_ior = entering ? 1.0 : mtl->ior;
		to_ior = entering ? mtl->ior : 1.0;

		rdir = ray->dir;
		NORMALIZE(rdir);

		rray.orig = sp->pos;
		refract(&rray.dir, &rdir, &sp->normal, from_ior, to_ior);
		rray.dir.x *= RAY_MAG;
		rray.dir.y *= RAY_MAG;
		rray.dir.z *= RAY_MAG;

		trace(&rray, depth + 1, &rcol);
		col->r += rcol.r * mtl->refr;
		col->g += rcol.g * mtl->refr;
		col->b += rcol.b * mtl->refr;
	}
}

/* calculate reflection vector */
void reflect(struct vec3 *res, struct vec3 *v, struct vec3 *n)
{
	double dot = v->x * n->x + v->y * n->y + v->z * n->z;
	res->x = -(2.0 * dot * n->x - v->x);
	res->y = -(2.0 * dot * n->y - v->y);
	res->z = -(2.0 * dot * n->z - v->z);
}

void refract(struct vec3 *res, struct vec3 *v, struct vec3 *n, double from_ior, double to_ior)
{
	double cos_inc, ior, radical, beta;
	struct vec3 neg_n;

	neg_n.x = -n->x;
	neg_n.y = -n->y;
	neg_n.z = -n->z;

	cos_inc = DOT(*v, neg_n);
	ior = from_ior / to_ior;

	radical = 1.0 + SQ(ior) * (SQ(cos_inc) - 1.0);
	if(radical < 0.0) {
		reflect(res, v, n);		/* total internal reflection */
		return;
	}

	beta = ior * cos_inc - sqrt(radical);

	res->x = v->x * ior + n->x * beta;
	res->y = v->y * ior + n->y * beta;
	res->z = v->z * ior + n->z * beta;
}

void cross_product(struct vec3 *res, struct vec3 *v1, struct vec3 *v2)
{
	struct vec3 tmp;
	tmp.x = v1->y * v2->z - v1->z * v2->y;
	tmp.y = v1->z * v2->x - v1->x * v2->z;
	tmp.z = v1->x * v2->y - v1->y * v2->x;
	*res = tmp;
}

/* determine the primary ray corresponding to the specified pixel (x, y) */
void get_primary_ray(struct ray *ray, int x, int y, int sample)
{
	float m[3][3];
	struct vec3 i, j = {0, 1, 0}, k, dir;

	k.x = cam.targ.x - cam.pos.x;
	k.y = cam.targ.y - cam.pos.y;
	k.z = cam.targ.z - cam.pos.z;
	NORMALIZE(k);

	cross_product(&i, &j, &k);
	cross_product(&j, &k, &i);
	m[0][0] = i.x; m[0][1] = j.x; m[0][2] = k.x;
	m[1][0] = i.y; m[1][1] = j.y; m[1][2] = k.y;
	m[2][0] = i.z; m[2][1] = j.z; m[2][2] = k.z;

	get_sample_pos(&ray->dir, x, y, sample);
	ray->dir.z = 2.0 / cam.fov;
	ray->dir.x *= RAY_MAG;
	ray->dir.y *= RAY_MAG;
	ray->dir.z *= RAY_MAG;

	dir = ray->dir;
	ray->dir.x = dir.x * m[0][0] + dir.y * m[0][1] + dir.z * m[0][2];
	ray->dir.y = dir.x * m[1][0] + dir.y * m[1][1] + dir.z * m[1][2];
	ray->dir.z = dir.x * m[2][0] + dir.y * m[2][1] + dir.z * m[2][2];

	ray->orig.x = cam.pos.x;
	ray->orig.y = cam.pos.y;
	ray->orig.z = cam.pos.z;
}


void get_sample_pos(struct vec3 *res, int x, int y, int sample)
{
	float px, py;
	struct vec3 jt;

	px = x;
	py = y;

	if(sample) {
		jitter(&jt, x, y, sample);
		px += jt.x;
		py += jt.y;
	}

	res->x = (px / xres * 2.0 - 1.0) * aspect;
	res->y = 1.0 - py / yres * 2.0;
	res->z = 0.0;
}

/* jitter function taken from Graphics Gems I. */
void jitter(struct vec3 *pt, int x, int y, int s) {
	pt->x = urand[(x + (y << 2) + irand[(x + s) & JITTER_MASK]) & JITTER_MASK].x;
	pt->y = urand[(y + (x << 2) + irand[(y + s) & JITTER_MASK]) & JITTER_MASK].y;
}

int ray_object(struct object *obj, struct ray *ray, struct spoint *sp)
{
	switch(obj->type) {
	case OBJ_SPHERE:
		if(ray_sphere(&obj->o.sph, ray, sp)) {
			if(sp) sp->obj = obj;
			return 1;
		}
		break;

	case OBJ_PLANE:
		if(ray_plane(&obj->o.pln, ray, sp)) {
			if(sp) sp->obj = obj;
			return 1;
		}
		break;
	}
	return 0;
}

/* Calculate ray-sphere intersection, and return {1, 0} to signify hit or no hit.
 * Also the surface point parameters like position, normal, etc are returned through
 * the sp pointer if it is not NULL.
 */
int ray_sphere(struct sphere *sph, struct ray *ray, struct spoint *sp)
{
	double a, b, c, d, sqrt_d, t1, t2;

	a = SQ(ray->dir.x) + SQ(ray->dir.y) + SQ(ray->dir.z);
	b = 2.0 * ray->dir.x * (ray->orig.x - sph->pos.x) +
				2.0 * ray->dir.y * (ray->orig.y - sph->pos.y) +
				2.0 * ray->dir.z * (ray->orig.z - sph->pos.z);
	c = SQ(sph->pos.x) + SQ(sph->pos.y) + SQ(sph->pos.z) +
				SQ(ray->orig.x) + SQ(ray->orig.y) + SQ(ray->orig.z) +
				2.0 * (-sph->pos.x * ray->orig.x - sph->pos.y * ray->orig.y - sph->pos.z * ray->orig.z) - SQ(sph->rad);

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

		sp->pos.x = ray->orig.x + ray->dir.x * sp->dist;
		sp->pos.y = ray->orig.y + ray->dir.y * sp->dist;
		sp->pos.z = ray->orig.z + ray->dir.z * sp->dist;

		sp->normal.x = (sp->pos.x - sph->pos.x) / sph->rad;
		sp->normal.y = (sp->pos.y - sph->pos.y) / sph->rad;
		sp->normal.z = (sp->pos.z - sph->pos.z) / sph->rad;

		reflect(&sp->vref, &ray->dir, &sp->normal);
		NORMALIZE(sp->vref);
	}
	return 1;
}

/* calculate ray-plane intersection... etc */
int ray_plane(struct plane *pln, struct ray *ray, struct spoint *sp)
{
	struct vec3 to_orig;
	double t;

	if(fabs(DOT(pln->norm, ray->dir)) < ERR_MARGIN) {
		return 0;
	}

	to_orig.x = ray->orig.x - pln->pos.x;
	to_orig.y = ray->orig.y - pln->pos.y;
	to_orig.z = ray->orig.z - pln->pos.z;
	t = -DOT(pln->norm, to_orig) / DOT(pln->norm, ray->dir);

	if(t < ERR_MARGIN || t > 1.0) {
		return 0;
	}

	if(sp) {
		sp->dist = t;

		sp->pos.x = ray->orig.x + ray->dir.x * t;
		sp->pos.y = ray->orig.y + ray->dir.y * t;
		sp->pos.z = ray->orig.z + ray->dir.z * t;

		sp->normal.x = pln->norm.x;
		sp->normal.y = pln->norm.y;
		sp->normal.z = pln->norm.z;

		reflect(&sp->vref, &ray->dir, &sp->normal);
		NORMALIZE(sp->vref);
	}
	return 1;
}

int ray_aabox(struct aabox *box, struct ray *ray)
{
	double tmin, tmax, tx1, tx2, ty1, ty2, tz1, tz2;

	tx1 = (box->vmin.x - ray->orig.x) * ray->invdir.x;
	tx2 = (box->vmax.x - ray->orig.x) * ray->invdir.x;

	tmin = MIN(tx1, tx2);
	tmax = MAX(tx1, tx2);

	ty1 = (box->vmin.y - ray->orig.y) * ray->invdir.y;
	ty2 = (box->vmax.y - ray->orig.y) * ray->invdir.y;

	tmin = MAX(tmin, MIN(ty1, ty2));
	tmax = MIN(tmax, MAX(ty1, ty2));

	tz1 = (box->vmin.z - ray->orig.z) * ray->invdir.z;
	tz2 = (box->vmax.z - ray->orig.z) * ray->invdir.z;

	tmin = MAX(tmin, MIN(tz1, tz2));
	tmax = MIN(tmax, MAX(tz1, tz2));

	return tmax >= tmin && tmax >= 0.0 && tmin <= 1.0;
}

int box_sphere(struct aabox *box, struct sphere *sph)
{
	struct vec3 p, b, bpos;
	double dmax;

	bpos.x = (box->vmin.x + box->vmax.x) * 0.5;
	bpos.y = (box->vmin.y + box->vmax.y) * 0.5;
	bpos.z = (box->vmin.z + box->vmax.z) * 0.5;

	b.x = (box->vmax.x - box->vmin.x) * 0.5;
	b.y = (box->vmax.y - box->vmin.y) * 0.5;
	b.z = (box->vmax.z - box->vmin.z) * 0.5;

	p.x = fabs(sph->pos.x - bpos.x) - b.x;
	p.y = fabs(sph->pos.y - bpos.y) - b.y;
	p.z = fabs(sph->pos.z - bpos.z) - b.z;

	if(p.x > 0.0 && p.y > 0.0 && p.z > 0.0) {
		double rsq = sph->rad * sph->rad;
		return p.x * p.x + p.y * p.y + p.z * p.z <= rsq;
	}

	dmax = MAX(p.x, p.y);
	if(p.z > dmax) dmax = p.z;

	return dmax <= sph->rad;
}

struct material *find_material(const char *name)
{
	int i;
	for(i=0; i<num_mtl; i++) {
		if(strcmp(materials[i].name, name) == 0) {
			return materials + i;
		}
	}
	return 0;
}

#define APPEND(arr, suffix) \
	do { \
		if(num_##suffix >= max_##suffix) { \
			int newsz = max_##suffix * 2; \
			void *tmp = realloc(arr, newsz * sizeof *arr); \
			if(!tmp) { \
				fprintf(stderr, "load_scene: failed to resize " #arr " to %d elements\n", newsz); \
				return -1; \
			} \
			arr = tmp; \
			max_##suffix = newsz; \
		} \
	} while(0)

char *expect_num(char *ptr, double *num)
{
	char *endp;
	*num = strtod(ptr, &endp);
	return endp == ptr ? 0 : endp;
}

char *expect_vec(char *ptr, struct vec3 *v)
{
	if(!(ptr = expect_num(ptr, &v->x))) return 0;
	if(!(ptr = expect_num(ptr, &v->y))) return 0;
	if(!(ptr = expect_num(ptr, &v->z))) return 0;
	return ptr;
}

char *expect_color(char *ptr, struct color *c)
{
	struct vec3 v;
	if(!(ptr = expect_vec(ptr, &v))) return 0;
	c->r = v.x;
	c->g = v.y;
	c->b = v.z;
	c->a = 1.0f;
	return ptr;
}

char *expect_str(char *ptr, char *buf, int bufsz)
{
	while(*ptr && isspace(*ptr)) ptr++;
	if(*ptr++ != '"') return 0;
	while(*ptr && *ptr != '"' && --bufsz > 0) {
		*buf++ = *ptr++;
	}
	*buf = 0;
	return *ptr ? ptr + 1 : ptr;
}

#define EXPECT_NUM(p, n) \
	do { \
		if(!(p = expect_num(p, n))) { \
			if(verbose) fprintf(stderr, "%d: expected number\n", lineno); \
			continue; \
		} \
	} while(0)

#define EXPECT_VEC(p, v) \
	do { \
		if(!(p = expect_vec(p, v))) { \
			if(verbose) fprintf(stderr, "%d: expected vector\n", lineno); \
			continue; \
		} \
	} while(0)

#define EXPECT_COLOR(p, c) \
	do { \
		if(!(p = expect_color(p, c))) { \
			if(verbose) fprintf(stderr, "%d: expected color\n", lineno); \
			continue; \
		} \
	} while(0)

#define EXPECT_STR(p, b, s) \
	do { \
		if(!(p = expect_str(p, b, s))) { \
			if(verbose) fprintf(stderr, "%d: expected string\n", lineno); \
			continue; \
		} \
	} while(0)


/* Load the scene from an extremely simple scene description file */
int load_scene(FILE *fp)
{
	int lineno = 0;
	char buf[256], name[64], *line, *ptr, *next;
	int max_obj, max_mtl, max_lt, max_inf;
	struct object *obj;
	struct light *lt;
	struct material *mtl;
	double num;

	max_obj = max_mtl = max_lt = max_inf = 8;
	objects = malloc(max_obj * sizeof *objects);
	materials = malloc(max_mtl * sizeof *materials);
	infobj = malloc(max_inf * sizeof *infobj);
	lights = malloc(max_lt * sizeof *lights);

	if(!objects || !materials || !lights || !infobj) {
		perror("load_scene: failed to allocate initial buffer");
		return -1;
	}

	while(fgets(buf, sizeof buf, fp)) {
		lineno++;
		if((ptr = strchr(buf, '#'))) *ptr = 0;
		ptr = buf;
		while(*ptr && isspace(*ptr)) ptr++;
		if(!*ptr || *ptr == '\n' || *ptr == '\r') continue;

		line = ptr++;
		switch(line[0]) {
		case 'M':
			APPEND(materials, mtl);
			mtl = materials + num_mtl;

			EXPECT_STR(ptr, name, sizeof name);
			EXPECT_COLOR(ptr, &mtl->kd);
			EXPECT_COLOR(ptr, &mtl->ks);
			EXPECT_NUM(ptr, &mtl->spow);
			EXPECT_NUM(ptr, &mtl->refl);
			EXPECT_NUM(ptr, &mtl->refr);
			EXPECT_NUM(ptr, &mtl->ior);

			if(!(mtl->name = malloc(strlen(name) + 1))) {
				perror("failed to allocate material name buffer");
				return -1;
			}
			strcpy(mtl->name, name);
			num_mtl++;
			break;

		case 's':
		case 'p':
			APPEND(objects, obj);
			obj = objects + num_obj;

			if(line[0] == 's') {
				obj->type = OBJ_SPHERE;
				EXPECT_VEC(ptr, &obj->o.sph.pos);
				EXPECT_NUM(ptr, &obj->o.sph.rad);
			} else {
				obj->type = OBJ_PLANE;
				EXPECT_VEC(ptr, &obj->o.pln.pos);
				EXPECT_VEC(ptr, &obj->o.pln.norm);
			}
			if((next = expect_str(ptr, name, sizeof name))) {
				if(!(mtl = find_material(name))) {
					fprintf(stderr, "%d: reference to unknown material: \"%s\"\n", lineno, name);
					continue;
				}
				obj->mat = *mtl;
			} else {
				obj->mat = def_mtl;
				if(!(next = expect_color(ptr, &obj->mat.kd))) goto mtldone;
				if(!(next = expect_num(next, &obj->mat.spow))) goto mtldone;
				if(!(next = expect_num(next, &obj->mat.refl))) goto mtldone;
				if(!(next = expect_num(next, &obj->mat.refr))) goto mtldone;
				expect_num(next, &obj->mat.ior);
			}
mtldone:	num_obj++;

			if(line[0] == 'p') {
				/* also add planes to the infinite objects list */
				APPEND(infobj, inf);
				infobj[num_inf++] = *obj;
			}
			break;

		case 'l':
			APPEND(lights, lt);
			lt = lights + num_lt;

			EXPECT_VEC(ptr, &lt->pos);
			if(!expect_color(ptr, &lt->col)) {
				lt->col.r = lt->col.g = lt->col.b = 1.0f;
			}
			num_lt++;
			break;

		case 'c':
			EXPECT_VEC(ptr, &cam.pos);
			EXPECT_NUM(ptr, &cam.fov);
			cam.fov = M_PI * cam.fov / 180.0;
			EXPECT_VEC(ptr, &cam.targ);
			break;

		case '?':
			while(*ptr && isspace(*ptr)) ptr++;
			if(memcmp(ptr, "gamma", 5) == 0) {
				if(!(override & OPT_GAMMA)) {
					ptr += 5;
					EXPECT_NUM(ptr, &num);
					inv_gamma = 1.0 / num;
				}

			} else if(memcmp(ptr, "samples", 7) == 0) {
				if(!(override & OPT_SAMPLES)) {
					ptr += 7;
					EXPECT_NUM(ptr, &num);
					rays_per_pixel = (int)num;
				}

			} else if(memcmp(ptr, "bgcolor", 7) == 0) {
				if(!(override & OPT_BGCOL)) {
					ptr += 7;
					EXPECT_COLOR(ptr, &bgcol);
				}

			} else if(memcmp(ptr, "ambient", 7) == 0) {
				ptr += 7;
				EXPECT_COLOR(ptr, &ambcol);

			} else {
				fprintf(stderr, "%d: unknown option: %s\n", lineno, ptr);
			}
			break;

		default:
			if(verbose) {
				fprintf(stderr, "%d: ignoring unknown command: %c\n", lineno, line[0]);
			}
		}
	}

	return 0;
}

int build_octree(struct octnode *node, int depth)
{
	int i;
	struct octnode *sub;
	double hx, hy, hz;

	if(node->num_obj < OCT_MAX_OBJ || depth >= OCT_MAX_DEPTH) {
		return 0;
	}

	hx = (node->aabb.vmin.x + node->aabb.vmax.x) * 0.5;
	hy = (node->aabb.vmin.y + node->aabb.vmax.y) * 0.5;
	hz = (node->aabb.vmin.z + node->aabb.vmax.z) * 0.5;

	for(i=0; i<8; i++) {
		if(!(sub = malloc(sizeof *sub)) || !(sub->obj = malloc(node->num_obj * sizeof *sub->obj))) {
			perror("failed to allocate octree node");
			free(sub);
			return -1;
		}
		sub->num_obj = 0;

		sub->aabb.vmin.x = LERP(node->aabb.vmin.x, hx, i & 1);
		sub->aabb.vmax.x = LERP(hx, node->aabb.vmax.x, i & 1);
		sub->aabb.vmin.y = LERP(node->aabb.vmin.y, hy, (i >> 1) & 1);
		sub->aabb.vmax.y = LERP(hy, node->aabb.vmax.y, (i >> 1) & 1);
		sub->aabb.vmin.z = LERP(node->aabb.vmin.z, hz, (i >> 2) & 1);
		sub->aabb.vmax.z = LERP(hz, node->aabb.vmax.z, (i >> 2) & 1);

		if(!populate_octnode(sub, node->obj, node->num_obj)) {
			/* no objects intersect this node's bbox, drop it */
			free(sub->obj);
			free(sub);
			sub = 0;
		} else {
			if(build_octree(sub, depth + 1) == -1) {
				return -1;
			}
		}
		node->sub[i] = sub;
	}

	free(node->obj);
	node->obj = 0;
	node->num_obj = 0;
	return 0;
}

int populate_octnode(struct octnode *node, struct object **objects, int num_obj)
{
	int i;
	struct object *obj;

	for(i=0; i<num_obj; i++) {
		obj = *objects++;

		if(obj->type == OBJ_SPHERE) {
			if(box_sphere(&node->aabb, &obj->o.sph)) {
				node->obj[node->num_obj++] = obj;
			}
		}
	}
	return node->num_obj;
}

int isect_octnode(struct octnode *node, struct ray *ray, struct spoint *spret)
{
	int i;
	struct spoint sp, sp0;

	if(!ray_aabox(&node->aabb, ray)) return 0;

	if(!spret) {
		/* find any intersection */

		if(node->obj) {
			for(i=0; i<node->num_obj; i++) {
				if(ray_object(node->obj[i], ray, 0)) {
					return 1;
				}
			}
			return 0;
		}

		for(i=0; i<8; i++) {
			if(node->sub[i] && isect_octnode(node->sub[i], ray, 0)) {
				return 1;
			}
		}
		return 0;
	}

	sp0.dist = DBL_MAX;
	sp0.obj = 0;

	/* find the nearest intersection */
	if(node->obj) {
		/* leaf node, test all objects */
		for(i=0; i<node->num_obj; i++) {
			if(ray_object(node->obj[i], ray, &sp) && sp.dist < sp0.dist) {
				sp0 = sp;
			}
		}
	} else {
		/* non-leaf node, recurse to subnodes */
		for(i=0; i<8; i++) {
			if(!node->sub[i]) continue;

			if(isect_octnode(node->sub[i], ray, &sp) && sp.dist < sp0.dist) {
				sp0 = sp;
			}
		}
	}

	if(sp0.obj) {
		*spret = sp0;
		return 1;
	}
	return 0;
}

void print_obj(void)
{
	int i;
	struct object *obj = objects;

	for(i=0; i<num_obj; i++) {
		if(obj->type == OBJ_SPHERE) {
			fprintf(stderr, "sphere: (%.2f %.2f %.2f) %.2f\n",
					obj->o.sph.pos.x, obj->o.sph.pos.y, obj->o.sph.pos.z, obj->o.sph.rad);
		} else {
			fprintf(stderr, " plane: (%.2f %.2f %.2f) & (%.2f %.2f %.2f)\n",
					obj->o.pln.pos.x, obj->o.pln.pos.y, obj->o.pln.pos.z,
					obj->o.pln.norm.x, obj->o.pln.norm.y, obj->o.pln.norm.z);
		}
		obj++;
	}
}

void print_progr(int x, int max)
{
	int i;
	int progr = 100 * x / max;
	int count = progr / 2;

	fputs("rendering [", stderr);
	for(i=0; i<count-1; i++) {
		fputc('=', stderr);
	}
	if(count) fputc('>', stderr);
	for(i=count; i<50; i++) {
		fputc(' ', stderr);
	}
	fprintf(stderr, "] %d%%\r", progr);
}

/* provide a millisecond-resolution timer for each system */
#if defined(unix) || defined(__unix__)
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#ifdef __bsd__
#include <sys/sysctl.h>
#endif

unsigned long get_msec(void)
{
	static struct timeval tv, tv0;

	gettimeofday(&tv, 0);
	if(tv0.tv_sec == 0) {
		tv0 = tv;
		return 0;
	}
	return (tv.tv_sec - tv0.tv_sec) * 1000 + (tv.tv_usec - tv0.tv_usec) / 1000;
}

int detect_nproc(void)
{
#if defined(__bsd__)
	/* BSD systems provide the num.processors through sysctl */
	int num, mib[] = {CTL_HW, HW_NCPU};
	size_t len = sizeof num;

	sysctl(mib, 2, &num, &len, 0, 0);
	return num;

#elif defined(__sgi)
	/* SGI IRIX flavour of the _SC_NPROC_ONLN sysconf */
	return sysconf(_SC_NPROC_ONLN);
#else
	/* Linux (and others?) have the _SC_NPROCESSORS_ONLN sysconf */
	return sysconf(_SC_NPROCESSORS_ONLN);
#endif	/* bsd/sgi/other */
}

#elif defined(__WIN32__) || defined(WIN32)
#include <windows.h>

unsigned long get_msec(void)
{
	return GetTickCount();
}

int detect_nproc(void)
{
	/* under windows we need to call GetSystemInfo */
	SYSTEM_INFO info;
	GetSystemInfo(&info);
	return info.dwNumberOfProcessors;
}

#else
#error "I don't know how to measure time on your platform"
#endif
