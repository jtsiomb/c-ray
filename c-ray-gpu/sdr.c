/*
Copyright 2007 John Tsiombikas <nuclear@siggraph.org>

This file is part of the pixelshow 2007 invitation demo.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with the program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <GL/glew.h>

#if defined(unix) || defined(__unix__)
#include <unistd.h>
#include <sys/stat.h>
#endif	/* unix */

#include "sdr.h"

unsigned int create_vertex_shader(const char *src)
{
	return create_shader(src, GL_VERTEX_SHADER);
}

unsigned int create_pixel_shader(const char *src)
{
	return create_shader(src, GL_FRAGMENT_SHADER);
}

unsigned int create_shader(const char *src, unsigned int sdr_type)
{
	unsigned int sdr;
	int success, info_len;
	char *info_str = 0;
	GLenum err;

	sdr = glCreateShader(sdr_type);
	assert(glGetError() == GL_NO_ERROR);
	glShaderSource(sdr, 1, &src, 0);
	err = glGetError();
	assert(err == GL_NO_ERROR);
	glCompileShader(sdr);
	assert(glGetError() == GL_NO_ERROR);

	glGetShaderiv(sdr, GL_COMPILE_STATUS, &success);
	assert(glGetError() == GL_NO_ERROR);
	glGetShaderiv(sdr, GL_INFO_LOG_LENGTH, &info_len);
	assert(glGetError() == GL_NO_ERROR);

	if(info_len) {
		if((info_str = malloc(info_len + 1))) {
			glGetShaderInfoLog(sdr, info_len, 0, info_str);
			assert(glGetError() == GL_NO_ERROR);
		}
	}

	if(success) {
		fprintf(stderr, info_str ? "done: %s\n" : "done\n", info_str);
	} else {
		fprintf(stderr, info_str ? "failed: %s\n" : "failed\n", info_str);
		glDeleteShader(sdr);
		sdr = 0;
	}

	free(info_str);
	return sdr;
}

void free_shader(unsigned int sdr)
{
	glDeleteShader(sdr);
}

unsigned int load_vertex_shader(const char *fname)
{
	return load_shader(fname, GL_VERTEX_SHADER);
}

unsigned int load_pixel_shader(const char *fname)
{
	return load_shader(fname, GL_FRAGMENT_SHADER);
}

unsigned int load_shader(const char *fname, unsigned int sdr_type)
{
#if defined(unix) || defined(__unix__)
	struct stat st;
#endif
	unsigned int sdr;
	size_t filesize;
	FILE *fp;
	char *src;

	if(!(fp = fopen(fname, "r"))) {
		fprintf(stderr, "failed to open shader %s: %s\n", fname, strerror(errno));
		return 0;
	}

#if defined(unix) || defined(__unix__)
	fstat(fileno(fp), &st);
	filesize = st.st_size;
#else
	fseek(fp, 0, SEEK_END);
	filesize = ftell(fp);
	fseek(fp, 0, SEEK_SET);
#endif	/* unix */

	if(!(src = malloc(filesize + 1))) {
		fclose(fp);
		return 0;
	}
	fread(src, 1, filesize, fp);
	src[filesize] = 0;
	fclose(fp);

	fprintf(stderr, "compiling %s shader: %s... ", (sdr_type == GL_VERTEX_SHADER ? "vertex" : "pixel"), fname);
	sdr = create_shader(src, sdr_type);

	free(src);
	return sdr;
}


unsigned int get_vertex_shader(const char *fname)
{
	return get_shader(fname, GL_VERTEX_SHADER);
}

unsigned int get_pixel_shader(const char *fname)
{
	return get_shader(fname, GL_FRAGMENT_SHADER);
}

unsigned int get_shader(const char *fname, unsigned int sdr_type)
{
	unsigned int sdr;
#if 0
	if((res = get_resource(sdrman, fname))) {
		/* TODO: validate that the shader matches sdr_type? */
		return (uintptr_t)res;
	}
#endif

	if(!(sdr = load_shader(fname, sdr_type))) {
		return 0;
	}
	add_shader(fname, sdr);
	return sdr;
}

int add_shader(const char *fname, unsigned int sdr)
{
	return 0;/*add_resource(sdrman, fname, (void*)(uintptr_t)sdr);*/
}

int remove_shader(const char *fname)
{
	return 0;/*remove_resource(sdrman, fname);*/
}


/* ---- gpu programs ---- */

unsigned int create_program(void)
{
	unsigned int prog = glCreateProgram();
	assert(glGetError() == GL_NO_ERROR);
	return prog;
}

unsigned int create_program_link(unsigned int vs, unsigned int ps)
{
	unsigned int prog;

	if(!(prog = create_program())) {
		return 0;
	}
	
	attach_shader(prog, vs);
	assert(glGetError() == GL_NO_ERROR);
	attach_shader(prog, ps);
	assert(glGetError() == GL_NO_ERROR);
	
	if(link_program(prog) == -1) {
		free_program(prog);
		return 0;
	}
	return prog;
}

unsigned int create_program_load(const char *vfile, const char *pfile)
{
	unsigned int vs, ps;

	if(!(vs = get_vertex_shader(vfile)) || !(ps = get_pixel_shader(pfile))) {
		return 0;
	}
	return create_program_link(vs, ps);
}

void free_program(unsigned int sdr)
{
	glDeleteProgram(sdr);
}

void attach_shader(unsigned int prog, unsigned int sdr)
{
	glAttachShader(prog, sdr);
	assert(glGetError() == GL_NO_ERROR);
}

int link_program(unsigned int prog)
{
	int linked, info_len, retval = 0;
	char *info_str = 0;

	glLinkProgram(prog);
	assert(glGetError() == GL_NO_ERROR);
	glGetProgramiv(prog, GL_LINK_STATUS, &linked);
	assert(glGetError() == GL_NO_ERROR);
	glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &info_len);
	assert(glGetError() == GL_NO_ERROR);

	if(info_len) {
		if((info_str = malloc(info_len + 1))) {
			glGetProgramInfoLog(prog, info_len, 0, info_str);
			assert(glGetError() == GL_NO_ERROR);
		}
	}

	if(linked) {
		fprintf(stderr, info_str ? "linking done: %s\n" : "linking done\n", info_str);
	} else {
		fprintf(stderr, info_str ? "linking failed: %s\n" : "linking failed\n", info_str);
		retval = -1;
	}

	free(info_str);
	return retval;
}

int bind_program(unsigned int prog)
{
	GLenum err;

	glUseProgram(prog);
	if(prog && (err = glGetError()) != GL_NO_ERROR) {
		/* maybe the program is not linked, try to link first */
		if(err == GL_INVALID_OPERATION) {
			if(link_program(prog) == -1) {
				return -1;
			}
			glUseProgram(prog);
			return glGetError() == GL_NO_ERROR ? 0 : -1;
		}
		return -1;
	}
	return 0;
}

/* ugly but I'm not going to write the same bloody code over and over */
#define BEGIN_UNIFORM_CODE \
	int loc, curr_prog; \
	glGetIntegerv(GL_CURRENT_PROGRAM, &curr_prog); \
	if(curr_prog != prog && bind_program(prog) == -1) { \
		return -1; \
	} \
	if((loc = glGetUniformLocation(prog, name)) != -1)

#define END_UNIFORM_CODE \
	if(curr_prog != prog) { \
		bind_program(curr_prog); \
	} \
	return loc == -1 ? -1 : 0

int set_uniform_int(unsigned int prog, const char *name, int val)
{
	BEGIN_UNIFORM_CODE {
		glUniform1i(loc, val);
	}
	END_UNIFORM_CODE;
}

int set_uniform_float(unsigned int prog, const char *name, float val)
{
	BEGIN_UNIFORM_CODE {
		glUniform1f(loc, val);
	}
	END_UNIFORM_CODE;
}

int set_uniform_vec2(unsigned int prog, const char *name, vec2_t val)
{
	BEGIN_UNIFORM_CODE {
		glUniform2f(loc, val.x, val.y);
	}
	END_UNIFORM_CODE;
}

int set_uniform_vec3(unsigned int prog, const char *name, vec3_t val)
{
	BEGIN_UNIFORM_CODE {
		glUniform3f(loc, val.x, val.y, val.z);
	}
	END_UNIFORM_CODE;
}

int set_uniform_vec4(unsigned int prog, const char *name, vec4_t val)
{
	BEGIN_UNIFORM_CODE {
		glUniform4f(loc, val.x, val.y, val.z, val.w);
	}
	END_UNIFORM_CODE;
}

int set_uniform_mat4(unsigned int prog, const char *name, mat4_t val)
{
	BEGIN_UNIFORM_CODE {
		glUniformMatrix4fv(loc, 1, 1, (float*)val);
	}
	END_UNIFORM_CODE;
}


int get_attrib_loc(unsigned int prog, const char *name)
{
	int loc, curr_prog;
	
	glGetIntegerv(GL_CURRENT_PROGRAM, &curr_prog);
	if(curr_prog != prog && bind_program(prog) == -1) {
		return -1;
	}
	
	loc = glGetAttribLocation(prog, (char*)name);
	
	if(curr_prog != prog) {
		bind_program(curr_prog);
	}
	return loc;
}

void set_attrib_vec3(int attr_loc, vec3_t val)
{
	glVertexAttrib3f(attr_loc, val.x, val.y, val.z);
}

void set_attrib_float3(int attr_loc, float x, float y, float z)
{
	glVertexAttrib3f(attr_loc, x, y, z);
}
