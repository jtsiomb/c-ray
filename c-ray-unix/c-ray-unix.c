#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "byteorder.h"

#define msleep(x)	usleep((x) * 1000);

struct vec3 {
	double x, y, z;
};

struct ray {
	struct vec3 orig, dir;
};

struct material {
	struct vec3 col;	/* color */
	double spow;		/* specular power */
	double refl;		/* reflection intensity */
};

struct sphere {
	struct vec3 pos;
	double rad;
	struct material mat;
	struct sphere *next;
};

struct spoint {
	struct vec3 pos, normal, vref;	/* position, normal and view reflection */
	double dist;		/* parametric distance of intersection along the ray */
};

struct camera {
	struct vec3 pos, targ;
	double fov;
};

void redraw(void);
void render(int xsz, int ysz, uint32_t *fb);
struct vec3 trace(struct ray ray, int depth);
struct vec3 shade(struct sphere *obj, struct spoint *sp, int depth);
struct vec3 reflect(struct vec3 v, struct vec3 n);
struct vec3 cross_product(struct vec3 v1, struct vec3 v2);
struct ray get_primary_ray(int x, int y);
int ray_sphere(const struct sphere *sph, struct ray ray, struct spoint *sp);
void load_scene(const char *fname);
unsigned long get_msec(void);

#define MAX_LIGHTS		16				/* maximum number of lights */
#define INTERVAL		33				/* minimum frame update interval */
#define RAY_MAG			1000.0			/* trace rays of this magnitude */
#define MAX_RAY_DEPTH	5				/* raytrace recursion limit */
#define FOV				0.78539816		/* field of view in rads (pi/4) */
#define HALF_FOV		(FOV * 0.5)
#define ERR_MARGIN		1e-6			/* an arbitrary error margin to avoid surface acne */

/* bit-shift ammount for packing each color into a 32bit uint */
#ifdef LITTLE_ENDIAN
#define RSHIFT	16
#define BSHIFT	0
#else
#define RSHIFT	0
#define BSHIFT	16
#endif	/* endianess */
#define GSHIFT	8	/* this is the same in both byte orders */

/* some helpful macros... */
#define SQ(x)		((x) * (x))
#define MAX(a, b)	((a) > (b) ? (a) : (b))
#define DOT(a, b)	((a).x * (b).x + (a).y * (b).y + (a).z * (b).z)
#define NORMALIZE(a)  do {\
	double len = sqrt(DOT(a, a));\
	(a).x /= len; (a).y /= len; (a).z /= len;\
} while(0);

/* global state */
int xres = 400;
int yres = 300;
double aspect = 1.333333;
struct sphere *obj_list;
struct vec3 lights[MAX_LIGHTS];
int lnum = 0;
struct camera cam;
unsigned long start_time;
unsigned long frames = 0;
uint32_t *fb;
Display *dpy;
Window win;
XImage *ximg;
GC gc;

int main(int argc, char **argv) {
	int i, fullscreen = 0, limit_fps = 1;
	int done = 0;
	
	/* X Window System variables */
	int scr;
	Window root;
	Visual *vis;
	XSetWindowAttributes xattr;
	Atom wm_del;
	XTextProperty wname;
	XClassHint chint;
	char *title = "X11 real-time ray tracer";
	int mapped = 0;
	

	for(i=1; i<argc; i++) {
		if(argv[i][0] == '-' && argv[i][2] == 0) {
			char *sep;
			switch(argv[i][1]) {
			case 's':
				if(!isdigit(argv[++i][0]) || !(sep = strchr(argv[i], 'x')) || !isdigit(*(sep + 1))) {
					fprintf(stderr, "-s must be followed by something like \"640x480\"\n");
					return EXIT_FAILURE;
				}
				xres = atoi(argv[i]);
				yres = atoi(sep + 1);
				aspect = (double)xres / (double)yres;
				break;

			case 'f':
				fullscreen = 1;
				break;

			case 'l':
				limit_fps = !limit_fps;
				break;
			}
		}
	}

	/* initialize SDL */
	if(!(dpy = XOpenDisplay(0))) {
		fputs("Could not connect to the X server\n", stderr);
		return EXIT_FAILURE;
	}
	scr = DefaultScreen(dpy);
	root = RootWindow(dpy, scr);
	vis = DefaultVisual(dpy, scr);

	xattr.background_pixel = xattr.border_pixel = BlackPixel(dpy, scr);
	xattr.colormap = DefaultColormap(dpy, scr);

	win = XCreateWindow(dpy, root, 0, 0, xres, yres, 0, CopyFromParent, InputOutput, vis, CWColormap | CWBackPixel | CWBorderPixel, &xattr);
	XSelectInput(dpy, win, StructureNotifyMask | ExposureMask | KeyPressMask);

	wm_del = XInternAtom(dpy, "WM_DELETE_WINDOW", True);
	XSetWMProtocols(dpy, win, &wm_del, 1);

	XStringListToTextProperty(&title, 1, &wname);
	XSetWMName(dpy, win, &wname);
	XFree(wname.value);

	chint.res_name = "raytracer";
	chint.res_class = "rtwin";
	XSetClassHint(dpy, win, &chint);

	XMapWindow(dpy, win);
	XFlush(dpy);

	fb = malloc(xres * yres * sizeof(uint32_t));
	ximg = XCreateImage(dpy, vis, 24, ZPixmap, 0, (void*)fb, xres, yres, 8, 0);

	gc = XCreateGC(dpy, win, 0, 0);
	
	/* load the scene */
	load_scene("scene");
	start_time = get_msec();

	while(!done) {
		if(XPending(dpy)) {
			XEvent event;
			XNextEvent(dpy, &event);

			switch(event.type) {
			case MapNotify:
				mapped = 1;
				break;

			case UnmapNotify:
				mapped = 0;
				break;

			case Expose:
				if(event.xexpose.count == 0 && mapped) {
					redraw();
				}
				break;

			case KeyPress:
				{
					KeySym sym = XLookupKeysym((XKeyEvent*)&event.xkey, 0);
					if((sym & 0xff) == 27) {
						done = 1;
					}
				}
				break;

			case ClientMessage:
				{
					char *atom_name = XGetAtomName(dpy, event.xclient.message_type);
					if(!strcmp(atom_name, "WM_PROTOCOLS")) {
						done = 1;
					}
					XFree(atom_name);
				}
				break;

			default:
				break;
			}
		} else {
			redraw();
		}
	}

	XFreeGC(dpy, gc);
	XDestroyImage(ximg);
	XDestroyWindow(dpy, win);
	XCloseDisplay(dpy);

	{
		double running = (double)(get_msec() - start_time) / 1000.0;
		double fps = (double)frames / (double)running;	
		printf("running time: %f sec, frames rendered: %lu, fps: %f\n", running, frames, fps);
	}
	
	return 0;
}

void redraw(void) {
	render(xres, yres, fb);
	XPutImage(dpy, win, gc, ximg, 0, 0, 0, 0, ximg->width, ximg->height);
	XFlush(dpy);
	frames++;
}

/* render a frame of xsz/ysz dimensions into the provided framebuffer */
void render(int xsz, int ysz, uint32_t *fb) {
	int i, j;
	unsigned long time = get_msec() - start_time;
	double t = (double)time / 1000.0;

	/* perform some crude animation */
	double y = obj_list->next->pos.y;
	obj_list->next->pos.y = obj_list->next->pos.y + fabs(sin(t));
	lights[1].z = sin(t) * 400.0;
	lights[1].x = cos(t * 0.5) * 400.0;

	cam.pos.x = sin(t * 0.3) * 17;
	cam.pos.z = -cos(t * 0.3) * 17;
	cam.pos.y = sin(t * 0.6) * 2 + 4;
	
	/* for each pixel, trace a ray through the scene, then pack the color
	 * and put it into the framebuffer. (assumes contiguous scanlines).
	 */
	for(j=0; j<ysz; j++) {
		for(i=0; i<xsz; i++) {
			struct vec3 col = trace(get_primary_ray(i, j), 0);
			if(col.x > 1.0) col.x = 1.0;
			if(col.y > 1.0) col.y = 1.0;
			if(col.z > 1.0) col.z = 1.0;
			
			*fb++ = (uint8_t)(col.x * 255.0) << RSHIFT |
					(uint8_t)(col.y * 255.0) << GSHIFT |
					(uint8_t)(col.z * 255.0) << BSHIFT;
		}
	}

	obj_list->next->pos.y = y;
}

/* trace a ray throught the scene recursively (the recursion happens through
 * shade() to calculate reflection rays if necessary).
 */
struct vec3 trace(struct ray ray, int depth) {
	struct vec3 col;
	struct spoint sp, nearest_sp;
	struct sphere *nearest_obj = 0;
	struct sphere *iter = obj_list->next;

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
		double ispec, idiff;
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

struct vec3 reflect(struct vec3 v, struct vec3 n) {
	struct vec3 res;
	double dot = v.x * n.x + v.y * n.y + v.z * n.z;
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
	ray.dir.x = ((double)x / (double)xres) - 0.5;
	ray.dir.y = -(((double)y / (double)yres) - 0.65) / aspect;
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

/* Calculate ray-sphere intersection, and return {1, 0} to signify hit or no hit.
 * Also the surface point parameters like position, normal, etc are returned through
 * the sp pointer if it is not NULL.
 */
int ray_sphere(const struct sphere *sph, struct ray ray, struct spoint *sp) {
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

		sp->vref = reflect(ray.dir, sp->normal);
		NORMALIZE(sp->vref);
	}
	return 1;
}

/* Load the scene from an extremely simple scene description file */
#define DELIM	" \t\n"
void load_scene(const char *fname) {
	char line[256], *ptr, type;
	FILE *fp = fopen(fname, "r");

	obj_list = malloc(sizeof(struct sphere));
	obj_list->next = 0;
	
	while((ptr = fgets(line, 256, fp))) {
		int i;
		struct vec3 pos, col;
		double rad, spow, refl;
		
		while(*ptr == ' ' || *ptr == '\t') ptr++;
		if(*ptr == '#' || *ptr == '\n') continue;

		if(!(ptr = strtok(line, DELIM))) continue;
		type = *ptr;
		
		for(i=0; i<3; i++) {
			if(!(ptr = strtok(0, DELIM))) break;
			*((double*)&pos.x + i) = atof(ptr);
		}

		if(type == 'l') {
			lights[lnum++] = pos;
			continue;
		}

		if(!(ptr = strtok(0, DELIM))) continue;
		rad = atof(ptr);

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

		if(type == 's') {
			struct sphere *sph = malloc(sizeof *sph);
			sph->next = obj_list->next;
			obj_list->next = sph;

			sph->pos = pos;
			sph->rad = rad;
			sph->mat.col = col;
			sph->mat.spow = spow;
			sph->mat.refl = refl;
		} else {
			fprintf(stderr, "unknown type: %c\n", type);
		}
	}

	fclose(fp);
}

unsigned long get_msec(void) {
	static struct timeval timeval, first_timeval;
	
	gettimeofday(&timeval, 0);

	if(first_timeval.tv_sec == 0) {
		first_timeval = timeval;
		return 0;
	}
	return (timeval.tv_sec - first_timeval.tv_sec) * 1000 + (timeval.tv_usec - first_timeval.tv_usec) / 1000;
}
