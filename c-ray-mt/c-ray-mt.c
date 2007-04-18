/* c-ray-mt - a simple multithreaded raytracing filter.
 * Copyright (C) 2006 John Tsiombikas <nuclear@siggraph.org>
 *
 * You are free to use, modify and redistribute this program under the
 * terms of the GNU General Public License v2 or (at your option) later.
 * see "http://www.gnu.org/licenses/gpl.txt" for details.
 * ---------------------------------------------------------------------
 * Usage:
 *   compile:  just type make
 *              (add any arch-specific optimizations for your compiler in CFLAGS first)
 *       run:  cat scene | ./c-ray-mt [-t num-threads] >foo.ppm
 *              (on broken systems such as windows try: c-ray-mt -i scene -o foo.ppm)
 *     enjoy:  display foo.ppm
 *              (with imagemagick, or use your favorite image viewer)
 * ---------------------------------------------------------------------
 * Scene file format:
 *   # sphere (many)
 *   s  x y z  rad   r g b   shininess   reflectivity
 *   # light (many)
 *   l  x y z
 *   # camera (one)
 *   c  x y z  fov   tx ty tz
 * ---------------------------------------------------------------------
 */
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <pthread.h>
#include <kdtree.h>

#define VER_MAJOR	2
#define VER_MINOR	0
#define VER_STR		"c-ray-mt v%d.%d\n"

#if !defined(unix) && !defined(__unix__)
#ifdef __MACH__
#define unix		1
#define __unix__	1
#endif	/* __MACH__ */
#endif	/* unix */

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

typedef struct vec3 {
	double x, y, z;
} vec3_t;

#define IOR_STACK_SZ	8
typedef struct ray {
	vec3_t orig, dir;
	double ior_stack[IOR_STACK_SZ];
	int ior_stack_top;
} ray_t;

struct material {
	vec3_t col;	/* color */
	double spow;		/* specular power */
	double refl;		/* reflection intensity */
	double refr;		/* refraction (transmission) intensity */

	int use_cook_tor;	/* use cook-torrance model */
	double roughness;	/* C&T */
	double specularity;	/* C&T */
	double ior;			/* C&T */
};

struct sphere {
	vec3_t pos;
	double rad;
	struct material mat;
	struct sphere *next;
};

struct spoint {
	ray_t ray;
	vec3_t pos, normal, view, vref;	/* position, normal, view and view reflection */
	double dist;		/* parametric distance of intersection along the ray */
};

struct camera {
	vec3_t pos, targ;
	double fov;
};

struct light {
	vec3_t pos;
	double radius;
	double intensity;
	int photons;
};

struct photon {
	vec3_t pos;
	vec3_t dir;
	float r, g, b;
};

struct thread_data {
	pthread_t tid;
	int sl_start, sl_count;

	uint32_t *pixels;
};

int photon_pass(void);
void render_scanline(int xsz, int ysz, int sl, uint32_t *fb, int samples);
struct sphere *find_nearest_hit(ray_t ray, struct spoint *sp);
vec3_t trace(ray_t ray, int depth);
vec3_t shade(struct sphere *obj, struct spoint *sp, int depth);
vec3_t reflect(vec3_t v, vec3_t n);
vec3_t refract(vec3_t v, vec3_t n, double from_ior, double to_ior);
vec3_t cross_product(vec3_t v1, vec3_t v2);
ray_t get_primary_ray(int x, int y, int sample);
vec3_t get_sample_pos(int x, int y, int sample);
vec3_t jitter(int x, int y, int s);
int ray_sphere(const struct sphere *sph, ray_t ray, struct spoint *sp);
void load_scene(FILE *fp);
unsigned long get_msec(void);
void sighandler(int sig);

void *thread_func(void *tdata);

#define LT_ENERGY		100
#define GATHER_DIST		0.25

#define MAX_LIGHTS		16				/* maximum number of lights */
#define RAY_MAG			1000.0			/* trace rays of this magnitude */
#define MAX_RAY_DEPTH	5				/* raytrace recursion limit */
#define FOV				0.78539816		/* field of view in rads (pi/4) */
#define HALF_FOV		(FOV * 0.5)
#define ERR_MARGIN		1e-6			/* an arbitrary error margin to avoid surface acne */
#define PI				3.1415926535897931

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
#define CLAMP(x, a, b)	MIN(MAX(x, a), b)
#define DOT(a, b)	((a).x * (b).x + (a).y * (b).y + (a).z * (b).z)
#define NORMALIZE(a)  do {\
	double len = sqrt(DOT(a, a));\
	(a).x /= len; (a).y /= len; (a).z /= len;\
} while(0);

/* global state */
int xres = 800;
int yres = 600;
int rays_per_pixel = 1;
double aspect = 1.333333;
struct sphere *obj_list;
struct light lights[MAX_LIGHTS];
int lnum = 0;
struct camera cam;

int thread_num = 1;
struct thread_data *threads;
int scanlines_done;

int start = 0;
pthread_mutex_t start_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t start_cond = PTHREAD_COND_INITIALIZER;

#define NRAN	1024
#define MASK	(NRAN - 1)
vec3_t urand[NRAN];
int irand[NRAN];

void *kd;	/* pointer to the kd-tree */
int photons_per_light = 0;

const char *usage = {
	"Usage: c-ray-mt [options]\n"
	"  Reads a scene file from stdin, writes the image to stdout, and stats to stderr.\n\n"
	"Options:\n"
	"  -t <num>   how many threads to use (default: 1)\n"
	"  -s WxH     where W is the width and H the height of the image\n"
	"  -r <rays>  shoot <rays> rays per pixel (antialiasing)\n"
	"  -i <file>  read from <file> instead of stdin\n"
	"  -o <file>  write to <file> instead of stdout\n"
	"  -p <num>   enable photon mapping and specify number of photons per light\n"
	"  -h         this help screen\n\n"
};



int main(int argc, char **argv) {
	int i;
	unsigned long rend_time, start_time;
	uint32_t *pixels;
	double sl, sl_per_thread;
	FILE *infile = stdin, *outfile = stdout;

	for(i=1; i<argc; i++) {
		if(argv[i][0] == '-' && argv[i][2] == 0) {
			char *sep;
			switch(argv[i][1]) {
			case 't':
				if(!isdigit(argv[++i][0])) {
					fprintf(stderr, "-t mus be followed by the number of worker threads to spawn\n");
					return EXIT_FAILURE;
				}
				thread_num = atoi(argv[i]);
				if(!thread_num) {
					fprintf(stderr, "invalid number of threads specified: %d\n", thread_num);
					return EXIT_FAILURE;
				}
				break;
					
			case 's':
				if(!isdigit(argv[++i][0]) || !(sep = strchr(argv[i], 'x')) || !isdigit(*(sep + 1))) {
					fputs("-s must be followed by something like \"640x480\"\n", stderr);
					return EXIT_FAILURE;
				}
				xres = atoi(argv[i]);
				yres = atoi(sep + 1);
				aspect = (double)xres / (double)yres;
				break;

			case 'i':
				if(!(infile = fopen(argv[++i], "rb"))) {
					fprintf(stderr, "failed to open input file %s: %s\n", argv[i], strerror(errno));
					return EXIT_FAILURE;
				}
				break;

			case 'o':
				if(!(outfile = fopen(argv[++i], "wb"))) {
					fprintf(stderr, "failed to open output file %s: %s\n", argv[i], strerror(errno));
					return EXIT_FAILURE;
				}
				break;

			case 'r':
				if(!isdigit(argv[++i][0])) {
					fputs("-r must be followed by a number (rays per pixel)\n", stderr);
					return EXIT_FAILURE;
				}
				rays_per_pixel = atoi(argv[i]);
				break;

			case 'p':
				if(!isdigit(argv[i + 1][0])) {
					if(argv[i + 1][0] == '-' && argv[i + 1][2] == 0) {
						photons_per_light = 1000;
						break;
					} else {
						fputs("-p must be followed by a number (photons per light)\n", stderr);
						return EXIT_FAILURE;
					}
				}
				photons_per_light = atoi(argv[++i]);
				break;

			case 'h':
				fputs(usage, stdout);
				return 0;
				
			default:
				fprintf(stderr, "unrecognized argument: %s\n", argv[i]);
				fputs(usage, stderr);
				return EXIT_FAILURE;
			}
		} else {
			fprintf(stderr, "unrecognized argument: %s\n", argv[i]);
			fputs(usage, stderr);
			return EXIT_FAILURE;
		}
	}

	if(!(pixels = malloc(xres * yres * sizeof *pixels))) {
		perror("pixel buffer allocation failed");
		return EXIT_FAILURE;
	}
	load_scene(infile);


	/* initialize the random number tables for the jitter */
	for(i=0; i<NRAN; i++) urand[i].x = (double)rand() / RAND_MAX - 0.5;
	for(i=0; i<NRAN; i++) urand[i].y = (double)rand() / RAND_MAX - 0.5;
	for(i=0; i<NRAN; i++) irand[i] = (int)(NRAN * ((double)rand() / RAND_MAX));

	if(thread_num > yres) {
		fprintf(stderr, "more threads than scanlines specified, reducing number of threads to %d\n", yres);
		thread_num = yres;
	}

	if(!(threads = malloc(thread_num * sizeof *threads))) {
		perror("failed to allocate thread table");
		return EXIT_FAILURE;
	}

	sl = 0.0;
	sl_per_thread = (double)yres / (double)thread_num;
	for(i=0; i<thread_num; i++) {
		threads[i].sl_start = (int)sl;
		sl += sl_per_thread;
		threads[i].sl_count = (int)sl - threads[i].sl_start;
		threads[i].pixels = pixels;
		
		if(pthread_create(&threads[i].tid, 0, thread_func, &threads[i]) != 0) {
			perror("failed to spawn thread");
			return EXIT_FAILURE;
		}
	}
	threads[thread_num - 1].sl_count = yres - threads[thread_num - 1].sl_start;

	fprintf(stderr, VER_STR, VER_MAJOR, VER_MINOR);
	signal(SIGINT, sighandler);

	if(photons_per_light > 0) {
		fprintf(stderr, "building the photon map (%d photons)... ", photons_per_light * lnum);
		if(photon_pass() == -1) {
			fprintf(stderr, "error while building the photon map\n");
			return EXIT_FAILURE;
		}
		fprintf(stderr, "done\n");
	}
	
	pthread_mutex_lock(&start_mutex);
	start_time = get_msec();
	start = 1;
	pthread_cond_broadcast(&start_cond);
	pthread_mutex_unlock(&start_mutex);

	for(i=0; i<thread_num; i++) {
		pthread_join(threads[i].tid, 0);
	}
	rend_time = get_msec() - start_time;
	
	/* output statistics to stderr */
	fprintf(stderr, "Rendering took: %lu seconds (%lu milliseconds)\n", rend_time / 1000, rend_time);

	/* output the image */
	fprintf(outfile, "P6\n%d %d\n255\n", xres, yres);
	for(i=0; i<xres * yres; i++) {
		fputc((pixels[i] >> RSHIFT) & 0xff, outfile);
		fputc((pixels[i] >> GSHIFT) & 0xff, outfile);
		fputc((pixels[i] >> BSHIFT) & 0xff, outfile);
	}
	fflush(outfile);

	if(infile != stdin) fclose(infile);
	if(outfile != stdout) fclose(outfile);

	if(kd) kd_free(kd);
	return 0;
}

void trace_photon(ray_t ray, float r, float g, float b, int iter)
{
	struct sphere *obj;
	struct spoint sp;

	if((obj = find_nearest_hit(ray, &sp))) {
		double rand_val = (double)rand() / (double)RAND_MAX;

		NORMALIZE(ray.dir);

		if(iter > MAX_RAY_DEPTH || rand_val > obj->mat.refl + obj->mat.refr) {
			struct photon *p = malloc(sizeof *p);
			p->pos = sp.pos;
			p->dir = ray.dir;
			p->r = r;
			p->g = g;
			p->b = b;
			
			kd_insert(kd, p->pos.x, p->pos.y, p->pos.z, p);
			return;
		}

		if(rand_val < obj->mat.refl) {
			ray.dir = reflect(ray.dir, sp.normal);
		} else {
			double from_ior, to_ior;

			if(DOT(ray.dir, sp.normal) > 0.0) {
				from_ior = obj->mat.ior;
				to_ior = 1.0;
				sp.normal.x = -sp.normal.x;
				sp.normal.y = -sp.normal.y;
				sp.normal.z = -sp.normal.z;
			} else {
				from_ior = 1.0;
				to_ior = obj->mat.ior;
			}

			ray.dir = refract(ray.dir, sp.normal, from_ior, to_ior);
		}
		
		ray.orig = sp.pos;
		ray.dir.x *= RAY_MAG;
		ray.dir.y *= RAY_MAG;
		ray.dir.z *= RAY_MAG;
		trace_photon(ray, r, g, b, iter + 1);
	}
}

/* called before rendering if photon mapping is enabled to build the photon map */
int photon_pass(void)
{
	int i, j;

	if(!(kd = kd_create())) {
		return -1;
	}
	kd_data_destructor(kd, free);

	for(i=0; i<lnum; i++) {
		double energy = (double)LT_ENERGY / (double)lights[i].photons;
		for(j=0; j<lights[i].photons; j++) {
			ray_t ray;

			ray.dir.x = RAY_MAG * (((double)rand() / (double)RAND_MAX) * 2.0 - 1.0);
			ray.dir.y = RAY_MAG * (((double)rand() / (double)RAND_MAX) * 2.0 - 1.0);
			ray.dir.z = RAY_MAG * (((double)rand() / (double)RAND_MAX) * 2.0 - 1.0);
			ray.orig = lights[i].pos;

			trace_photon(ray, energy, energy, energy, 0);
		}
	}

	return 0;
}

/* render a frame of xsz/ysz dimensions into the provided framebuffer */
void render_scanline(int xsz, int ysz, int sl, uint32_t *fb, int samples) {
	int i, s;
	double rcp_samples = 1.0 / (double)samples;

	for(i=0; i<xsz; i++) {
		double r, g, b;
		r = g = b = 0.0;
			
		for(s=0; s<samples; s++) {
			vec3_t col = trace(get_primary_ray(i, sl, s), 0);
			r += col.x;
			g += col.y;
			b += col.z;
		}

		r = r * rcp_samples;
		g = g * rcp_samples;
		b = b * rcp_samples;
			
		fb[sl * xsz + i] = ((uint32_t)(MIN(r, 1.0) * 255.0) & 0xff) << RSHIFT |
							((uint32_t)(MIN(g, 1.0) * 255.0) & 0xff) << GSHIFT |
							((uint32_t)(MIN(b, 1.0) * 255.0) & 0xff) << BSHIFT;
	}
}

struct sphere *find_nearest_hit(ray_t ray, struct spoint *sp)
{
	struct spoint nearest_sp;
	struct sphere *nearest_obj = 0;
	struct sphere *iter = obj_list->next;

	while(iter) {
		if(ray_sphere(iter, ray, sp)) {
			if(!nearest_obj || sp->dist < nearest_sp.dist) {
				nearest_obj = iter;
				nearest_sp = *sp;
			}
		}
		iter = iter->next;
	}

	*sp = nearest_sp;
	return nearest_obj;
}

/* trace a ray throught the scene recursively (the recursion happens through
 * shade() to calculate reflection rays if necessary).
 */
vec3_t trace(ray_t ray, int depth) {
	vec3_t col;
	struct spoint sp;
	struct sphere *obj;

	/* if we reached the recursion limit, bail out */
	if(depth >= MAX_RAY_DEPTH) {
		col.x = col.y = col.z = 0.0;
		return col;
	}
	
	if((obj = find_nearest_hit(ray, &sp))) {
		col = shade(obj, &sp, depth);
	} else {
		col.x = col.y = col.z = 0.0;
	}

	return col;
}

vec3_t calc_irradiance(vec3_t pos, vec3_t norm, float max_dist)
{
	void *kdres;
	vec3_t irrad = {0, 0, 0};
	int i, sz;
	float tmp;

	if(!(kdres = kd_nearest_range(kd, pos.x, pos.y, pos.z, max_dist))) {
		fprintf(stderr, "kd_nearest_range returned 0! failed to allocate memory?\n");
		exit(EXIT_FAILURE);
	}
	sz = kd_res_size(kd);

	if(sz < 8) {
		return irrad;
	}

	for(i=0; i<sz; i++) {
		struct photon *ph = kd_res_item_data(kdres);
		irrad.x += ph->r;
		irrad.y += ph->g;
		irrad.z += ph->b;
		kd_res_next(kdres);
	}

	tmp = (1.0 / PI) / (max_dist * max_dist);
	irrad.x *= tmp;
	irrad.y *= tmp;
	irrad.z *= tmp;

	return irrad;
}

/* Calculates direct illumination with the phong reflectance model.
 * Also handles reflections by calling trace again, if necessary.
 */
vec3_t shade(struct sphere *obj, struct spoint *sp, int depth) {
	int i, entering;
	vec3_t col = {0, 0, 0};

	if(DOT(sp->ray.dir, sp->normal) > 0.0) {
		sp->normal.x = -sp->normal.x;
		sp->normal.y = -sp->normal.y;
		sp->normal.z = -sp->normal.z;
		entering = 0;
	} else {
		entering = 1;
	}

	/* for all lights ... */
	for(i=0; i<lnum; i++) {
		double ispec, idiff;
		vec3_t ldir;
		ray_t shadow_ray;
		struct sphere *iter = obj_list->next;
		int in_shadow = 0;

		ldir.x = lights[i].pos.x - sp->pos.x;
		ldir.y = lights[i].pos.y - sp->pos.y;
		ldir.z = lights[i].pos.z - sp->pos.z;

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

#ifdef BUILD_COOK_TORRANCE
			if(obj->mat.use_cook_tor) {
				vec3_t half;
				double ndoth, ndotv, ndotl, vdoth, ndoth_sq, sin2_ang, tan2_ang;
				double geom_a, geom_b, geom;
				double ftmp, d_mf;
				double fres, c, g;
				double m = MAX(obj->mat.roughness, 0.0001);

				/* half vector */
				half.x = sp->view.x + ldir.x;
				half.y = sp->view.y + ldir.y;
				half.z = sp->view.z + ldir.z;
				NORMALIZE(half);

				/* calc various useful dot products */
				ndoth = MAX(DOT(sp->normal, half), 0.0);
				ndotv = MAX(DOT(sp->normal, sp->view), 0.0);
				ndotl = DOT(sp->normal, ldir);
				vdoth = MAX(DOT(sp->view, half), 0.0);
				ndoth_sq = SQ(ndoth);

				/* geometric term */
				geom_a = (2.0 * ndoth * ndotv) / vdoth;
				geom_b = (2.0 * ndoth * ndotl) / vdoth;
				geom = MIN(1.0, MIN(geom_a, geom_b));

				/* beckmann microfacet distribution term */
				sin2_ang = 1.0 - ndoth_sq;		/* sin^2(a) = 1.0 - cos^2(a) */
				tan2_ang = sin2_ang / ndoth_sq;	/* tan^2(a) = sin^2(a) / cos^2(a) */
				d_mf = exp(-tan2_ang / SQ(m)) / (SQ(m) * SQ(ndoth_sq));

				/* fresnel term */
				c = vdoth;
				g = sqrt(SQ(obj->mat.ior) + SQ(vdoth) - 1.0);
				ftmp = (c * (g + c) - 1.0) / (c * (g - c) + 1.0);
				fres = 0.5 * (SQ(g - c) / SQ(g + c)) * (1.0 + SQ(ftmp));

				/* specular and diffuse components */
				ispec = obj->mat.specularity * ((fres / PI) * (d_mf / ndotl) * (geom / ndotv));
				idiff = (1.0 - obj->mat.specularity) * MAX(DOT(sp->normal, ldir), 0.0);

				/* add calculated lighting to pixel color */
				col.x += idiff * obj->mat.col.x + ispec * obj->mat.col.x;
				col.y += idiff * obj->mat.col.y + ispec * obj->mat.col.y;
				col.z += idiff * obj->mat.col.z + ispec * obj->mat.col.z;

			} else {	/* use phong */
#endif
				idiff = MAX(DOT(sp->normal, ldir), 0.0);
				ispec = obj->mat.spow > 0.0 ? pow(MAX(DOT(sp->vref, ldir), 0.0), obj->mat.spow) : 0.0;

				col.x += idiff * obj->mat.col.x + ispec;
				col.y += idiff * obj->mat.col.y + ispec;
				col.z += idiff * obj->mat.col.z + ispec;
#ifdef BUILD_COOK_TORRANCE
			}
#endif
		}
	}

	/* use the photon map */
	if(kd) {
		vec3_t irrad;

		irrad = calc_irradiance(sp->pos, sp->normal, GATHER_DIST);
		col.x += irrad.x;
		col.y += irrad.y;
		col.z += irrad.z;
	}

	/* Also, if the object is reflective, spawn a reflection ray, and call trace()
	 * to calculate the light arriving from the mirror direction.
	 */
	if(obj->mat.refl > 0.0) {
		ray_t ray;
		vec3_t rcol;

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

	if(obj->mat.refr > 0.0) {
		ray_t ray;
		vec3_t rdir, rcol;

		double from_ior = entering ? 1.0 : obj->mat.ior;
		double to_ior = entering ? obj->mat.ior : 1.0;

		rdir = sp->ray.dir;
		NORMALIZE(rdir);

		ray.orig = sp->pos;
		ray.dir = refract(rdir, sp->normal, from_ior, to_ior);
		ray.dir.x *= RAY_MAG;
		ray.dir.y *= RAY_MAG;
		ray.dir.z *= RAY_MAG;

		rcol = trace(ray, depth + 1);
		col.x += rcol.x * obj->mat.refr;
		col.y += rcol.y * obj->mat.refr;
		col.z += rcol.z * obj->mat.refr;
	}

	return col;
}

/* calculate reflection vector */
vec3_t reflect(vec3_t v, vec3_t n) {
	vec3_t res;
	double dot = v.x * n.x + v.y * n.y + v.z * n.z;
	res.x = -(2.0 * dot * n.x - v.x);
	res.y = -(2.0 * dot * n.y - v.y);
	res.z = -(2.0 * dot * n.z - v.z);
	return res;
}

vec3_t refract(vec3_t v, vec3_t n, double from_ior, double to_ior)
{
	double cos_inc, ior, radical, beta;
	vec3_t neg_norm;

	neg_norm.x = -n.x;
	neg_norm.y = -n.y;
	neg_norm.z = -n.z;
	
	cos_inc = DOT(v, neg_norm);
	ior = from_ior / to_ior;
	
	radical = 1.0 + SQ(ior) * (SQ(cos_inc) - 1.0);
	if(radical < 0.0) {
		return reflect(v, n);
	}

	beta = ior * cos_inc - sqrt(radical);

	v.x = v.x * ior + n.x * beta;
	v.y = v.y * ior + n.y * beta;
	v.z = v.z * ior + n.z * beta;
	return v;
}

vec3_t cross_product(vec3_t v1, vec3_t v2) {
	vec3_t res;
	res.x = v1.y * v2.z - v1.z * v2.y;
	res.y = v1.z * v2.x - v1.x * v2.z;
	res.z = v1.x * v2.y - v1.y * v2.x;
	return res;
}

/* determine the primary ray corresponding to the specified pixel (x, y) */
ray_t get_primary_ray(int x, int y, int sample) {
	ray_t ray;
	double m[3][3];
	vec3_t i, j = {0, 1, 0}, k, dir, orig, foo;

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


vec3_t get_sample_pos(int x, int y, int sample) {
	vec3_t pt;
	static double sf = 0.0;

	if(sf == 0.0) {
		sf = 1.5 / (double)xres;
	}

	pt.x = ((double)x / (double)xres) - 0.5;
	pt.y = -(((double)y / (double)yres) - 0.65) / aspect;
	pt.z = 0;

	if(sample) {
		vec3_t jt = jitter(x, y, sample);
		pt.x += jt.x * sf;
		pt.y += jt.y * sf / aspect;
	}
	return pt;
}

/* jitter function taken from Graphics Gems I. */
vec3_t jitter(int x, int y, int s) {
	vec3_t pt;
	pt.x = urand[(x + (y << 2) + irand[(x + s) & MASK]) & MASK].x;
	pt.y = urand[(y + (x << 2) + irand[(y + s) & MASK]) & MASK].y;
	pt.z = 0;
	return pt;
}

/* Calculate ray-sphere intersection, and return {1, 0} to signify hit or no hit.
 * Also the surface point parameters like position, normal, etc are returned through
 * the sp pointer if it is not NULL.
 */
int ray_sphere(const struct sphere *sph, ray_t ray, struct spoint *sp) {
	double a, b, c, d, sqrt_d, t1, t2;
	
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

		sp->view.x = -ray.dir.x;
		sp->view.y = -ray.dir.y;
		sp->view.z = -ray.dir.z;
		sp->vref = reflect(ray.dir, sp->normal);
		NORMALIZE(sp->vref);
		NORMALIZE(sp->view);

		sp->ray = ray;
	}
	return 1;
}

/* Load the scene from an extremely simple scene description file */
#define DELIM	" \t\n"
void load_scene(FILE *fp) {
	char line[256], *ptr, type;

	obj_list = malloc(sizeof(struct sphere));
	obj_list->next = 0;
	
	while((ptr = fgets(line, 256, fp))) {
		int i;
		vec3_t pos, col;
		double rad, spow, refl, refr, ior, rough;
		
		while(*ptr == ' ' || *ptr == '\t') ptr++;
		if(*ptr == '#' || *ptr == '\n') continue;

		if(!(ptr = strtok(line, DELIM))) continue;
		type = *ptr;
		
		for(i=0; i<3; i++) {
			if(!(ptr = strtok(0, DELIM))) break;
			*((double*)&pos.x + i) = atof(ptr);
		}

		if(!(ptr = strtok(0, DELIM))) continue;
		rad = atof(ptr);

		if(type == 'l') {
			lights[lnum].pos = pos;
			lights[lnum].radius = rad;
			lights[lnum++].photons = photons_per_light;
			continue;
		}

		for(i=0; i<3; i++) {
			if(!(ptr = strtok(0, DELIM))) break;
			*((double*)&col.x + i) = atof(ptr);
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

		if(!(ptr = strtok(0, DELIM))) continue;
		refr = atof(ptr);

		if(!(ptr = strtok(0, DELIM))) {
			ior = 1.0;
		} else {
			ior = atof(ptr);
		}

		if(type == 'm') {

			if(!(ptr = strtok(0, DELIM))) continue;
			rough = atof(ptr);
		}

		if(type == 's' || type == 'm') {
			struct sphere *sph = malloc(sizeof *sph);
			sph->next = obj_list->next;
			obj_list->next = sph;

			sph->pos = pos;
			sph->rad = rad;
			sph->mat.col = col;
			sph->mat.spow = spow;
			sph->mat.refl = refl;
			sph->mat.refr = refr;

			sph->mat.specularity = spow;
			sph->mat.ior = ior;
			sph->mat.roughness = rough;

			sph->mat.use_cook_tor = (type == 'm') ? 1 : 0;
		} else {
			fprintf(stderr, "unknown type: %c\n", type);
		}
	}
}


/* provide a millisecond-resolution timer for each system */
#if defined(unix) || defined(__unix__)
#include <time.h>
#include <sys/time.h>
unsigned long get_msec(void) {
	static struct timeval timeval, first_timeval;
	
	gettimeofday(&timeval, 0);
	if(first_timeval.tv_sec == 0) {
		first_timeval = timeval;
		return 0;
	}
	return (timeval.tv_sec - first_timeval.tv_sec) * 1000 + (timeval.tv_usec - first_timeval.tv_usec) / 1000;
}
#elif defined(__WIN32__) || defined(WIN32)
#include <windows.h>
unsigned long get_msec(void) {
	return GetTickCount();
}
#else
#error "I don't know how to measure time on your platform"
#endif

void *thread_func(void *tdata) {
	int i;
	struct thread_data *td = (struct thread_data*)tdata;

	pthread_mutex_lock(&start_mutex);
	while(!start) {
		pthread_cond_wait(&start_cond, &start_mutex);
	}
	pthread_mutex_unlock(&start_mutex);

	for(i=0; i<td->sl_count; i++) {
		render_scanline(xres, yres, i + td->sl_start, td->pixels, rays_per_pixel);
		scanlines_done++;
	}

	return 0;
}


void sighandler(int sig)
{
	static unsigned int last_int;

	if(sig == SIGINT) {
		unsigned int msec;
		fprintf(stderr, "rendered: %d scanlines out of %d\n", scanlines_done, yres);
		if((msec = get_msec()) - last_int < 1000) {
			fprintf(stderr, "aborting\n");
			exit(0);
		}
		last_int = msec;
	}
}
