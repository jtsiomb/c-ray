#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "scene.h"

struct sphere *obj_list;
int obj_num;
vec3_t lights[MAX_LIGHTS];
int lnum;
struct camera cam;

/* Load the scene from an extremely simple scene description file */
#define DELIM	" \t\n"
void load_scene(FILE *fp)
{
	char line[256], *ptr, type;

	if(!(obj_list = malloc(sizeof(struct sphere)))) {
		perror("failed to allocate memory");
		exit(1);
	}
	obj_list->next = 0;
	
	while((ptr = fgets(line, 256, fp))) {
		int i;
		vec3_t pos, col;
		float rad, spow, refl;
		
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

		if(type == 's') {
			struct sphere *sph;
			if(!(sph = malloc(sizeof *sph))) {
				perror("failed to allocate memory");
				exit(1);
			}
			sph->next = obj_list->next;
			obj_list->next = sph;

			sph->pos = pos;
			sph->rad = rad;
			sph->mat.col = col;
			sph->mat.spow = spow;
			sph->mat.refl = refl;

			obj_num++;
		} else {
			fprintf(stderr, "unknown type: %c\n", type);
		}
	}
}


