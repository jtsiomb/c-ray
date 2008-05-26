#include <stdio.h>
#include <stdlib.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include "sdr.h"
#include "scene.h"

void redraw(void);
void reshape(int x, int y);
void keyb(unsigned char key, int x, int y);
void set_cam_matrix(void);
int next_pow2(int x);
unsigned int create_geom_tex(void);

#define RAY_MAG			1000.0			/* trace rays of this magnitude */
#define MAX_RAY_DEPTH	5				/* raytrace recursion limit */
#define FOV				0.78539816		/* field of view in rads (pi/4) */
#define HALF_FOV		(FOV * 0.5)
#define ERR_MARGIN		1e-6			/* an arbitrary error margin to avoid surface acne */

int xsz, ysz;
unsigned int prog;
unsigned int geom_tex;

float aspect = 1.333333;

int main(int argc, char **argv)
{
	load_scene(stdin);

	glutInitWindowSize(800, 600);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutCreateWindow("GPU Ray tracing");

	xsz = glutGet(GLUT_WINDOW_WIDTH);
	ysz = glutGet(GLUT_WINDOW_HEIGHT);

	glutDisplayFunc(redraw);
	glutKeyboardFunc(keyb);
	glutReshapeFunc(reshape);

	glewInit();

	if(!(geom_tex = create_geom_tex())) {
		fprintf(stderr, "failed to create geometry texture\n");
		return 1;
	}

	if(!(prog = create_program_load("ray.vs.glsl", "ray.ps.glsl"))) {
		fprintf(stderr, "failed to load GPU program\n");
		return 1;
	}
	bind_program(prog);

	set_uniform_float(prog, "fov", 0.78539816);
	set_uniform_float(prog, "aspect", 1.33333333);
	set_uniform_int(prog, "obj_num", obj_num);
	set_uniform_int(prog, "geom_tex_height", next_pow2(obj_num));

	glutMainLoop();
	return 0;
}

void redraw(void)
{
	set_cam_matrix();

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, geom_tex);

	glBegin(GL_QUADS);
	glTexCoord2f(0, 0);
	glVertex2f(-1, -1);
	glTexCoord2f(1, 0);
	glVertex2f(1, -1);
	glTexCoord2f(1, 1);
	glVertex2f(1, 1);
	glTexCoord2f(0, 1);
	glVertex2f(-1, 1);
	glEnd();

	glDisable(GL_TEXTURE_2D);

	glutSwapBuffers();
}

void reshape(int x, int y)
{
	glViewport(0, 0, x, y);
	xsz = x;
	ysz = y;
}

void keyb(unsigned char key, int x, int y)
{
	if(key == 27) {
		exit(0);
	}
}

void set_cam_matrix(void)
{
	vec3_t i, j, k;
	mat4_t m;

	j = v3_cons(0, 1, 0);

	k.x = cam.targ.x - cam.pos.x;
	k.y = cam.targ.y - cam.pos.y;
	k.z = cam.targ.z - cam.pos.z;
	k = v3_normalize(k);

	i = v3_cross(j, k);
	j = v3_cross(k, i);

	m4_identity(m);
	m[0][0] = i.x; m[0][1] = j.x; m[0][2] = k.x;
	m[1][0] = i.y; m[1][1] = j.y; m[1][2] = k.y;
	m[2][0] = i.z; m[2][1] = j.z; m[2][2] = k.z;

	if(set_uniform_mat4(prog, "cam_mat", m) == -1 ||
			set_uniform_vec3(prog, "cam_pos", cam.pos) == -1) {
		fprintf(stderr, "failed to set camera properties\n");
	}
}

int next_pow2(int x)
{
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x + 1;
}

#define DTEX_WIDTH	4
unsigned int create_geom_tex(void)
{
	unsigned int tex;
	struct sphere *obj;
	float *data, *dptr;
	int ysz;

	/* create texture image for the scene data */
	ysz = next_pow2(obj_num);

	if(!(data = malloc(DTEX_WIDTH * ysz * 4 * sizeof *data))) {
		perror("failed to allocate memory");
		return 0;
	}
	memset(data, 0, DTEX_WIDTH * ysz * 4 * sizeof *data);
	dptr = data;

	obj = obj_list->next;
	while(obj) {
		/* first pixel: sphere position/radius */
		*dptr++ = obj->pos.x;
		*dptr++ = obj->pos.y;
		*dptr++ = obj->pos.z;
		*dptr++ = obj->rad;
		/* second pixel: material color */
		*dptr++ = obj->mat.col.x;
		*dptr++ = obj->mat.col.y;
		*dptr++ = obj->mat.col.z;
		*dptr++ = 1.0;
		/* third pixel: other material attributes */
		*dptr++ = obj->mat.spow;
		*dptr++ = obj->mat.refl;

		dptr += 6;	/* skip rest of the row */

		obj = obj->next;
	}

	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_2D, tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, DTEX_WIDTH, ysz, 0, GL_RGBA, GL_FLOAT, data);

	free(data);

	return tex;
}
