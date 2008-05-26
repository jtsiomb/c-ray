#ifndef _SDR_H_
#define _SDR_H_

#include <vmath.h>

#ifdef __cplusplus
extern "C" {
#endif	/* __cplusplus */

/* ---- shaders ---- */
unsigned int create_vertex_shader(const char *src);
unsigned int create_pixel_shader(const char *src);
unsigned int create_shader(const char *src, unsigned int sdr_type);
void free_shader(unsigned int sdr);

unsigned int load_vertex_shader(const char *fname);
unsigned int load_pixel_shader(const char *fname);
unsigned int load_shader(const char *src, unsigned int sdr_type);

unsigned int get_vertex_shader(const char *fname);
unsigned int get_pixel_shader(const char *fname);
unsigned int get_shader(const char *fname, unsigned int sdr_type);

int add_shader(const char *fname, unsigned int sdr);
int remove_shader(const char *fname);

/* ---- gpu programs ---- */
unsigned int create_program(void);
unsigned int create_program_link(unsigned int vs, unsigned int ps);
unsigned int create_program_load(const char *vfile, const char *pfile);
void free_program(unsigned int sdr);

void attach_shader(unsigned int prog, unsigned int sdr);
int link_program(unsigned int prog);
int bind_program(unsigned int prog);

int set_uniform_int(unsigned int prog, const char *name, int val);
int set_uniform_float(unsigned int prog, const char *name, float val);
int set_uniform_vec2(unsigned int prog, const char *name, vec2_t val);
int set_uniform_vec3(unsigned int prog, const char *name, vec3_t val);
int set_uniform_vec4(unsigned int prog, const char *name, vec4_t val);
int set_uniform_mat4(unsigned int prog, const char *name, mat4_t val);

int get_attrib_loc(unsigned int prog, const char *name);
void set_attrib_vec3(int attr_loc, vec3_t val);
void set_attrib_float3(int attr_loc, float x, float y, float z);

#ifdef __cplusplus
}
#endif	/* __cplusplus */

#endif	/* _SDR_H_ */
