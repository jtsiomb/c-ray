#include <math.h>
#include <float.h>
#include <assert.h>
#include "color.h"
#include "gfx.h"

typedef unsigned long uint32_t;
typedef unsigned short uint16_t;
typedef unsigned char uint8_t;


static int find_color(int r, int g, int b);
static int nearest_color(int r, int g, int b);

static struct palcol {
	unsigned char r, g, b;
} palette[256];

static int free_col_idx = 0;


void init_palette(void)
{
	unsigned int i;

	logfoo("setting up fixed palette\n");

	for(i=0; i<256; i++) {
		int r = i & 0xe0;
		int g = (i & 0x1c) << 3;
		int b = (i & 3) << 6;

		palette[i].r = r;
		palette[i].g = g;
		palette[i].b = b;

		set_palette(i, r, g, b);
	}

	free_col_idx = 256;
}


int alloc_color(int r, int g, int b)
{
	int idx = find_color(r, g, b);
	if(idx >= 0) {
		return idx;
	}

	if(free_col_idx >= 256) {
		return nearest_color(r, g, b);
	}
	
	idx = free_col_idx++;

	palette[idx].r = r;
	palette[idx].g = g;
	palette[idx].b = b;

	logfoo("alloc_color: %3d -> %3d %3d %3d\n", idx, r, g, b);

	set_palette(idx, r, g, b);
	return idx;
}

void lookup_color(int idx, int *r, int *g, int *b)
{
	*r = palette[idx].r;
	*g = palette[idx].g;
	*b = palette[idx].b;
}

static int find_color(int r, int g, int b)
{
	int i;

	for(i=0; i<256; i++) {
		if(palette[i].r == r && palette[i].g == g && palette[i].b == b) {
			return i;
		}
	}
	return -1;
}

static int nearest_color(int r, int g, int b)
{
	float fr, fg, fb, dsq0 = FLT_MAX;
	int i, nearest;

	fr = (float)r;
	fg = (float)g;
	fb = (float)b;

	for(i=0; i<256; i++) {
		float dr = fabs(fr - (float)palette[i].r);
		float dg = fabs(fg - (float)palette[i].g);
		float db = fabs(fb - (float)palette[i].b);

		float dist_sq = dr * dr + dg * dg + db * db;
		if(dist_sq < dsq0) {
			dsq0 = dist_sq;
			nearest = i;
		}
	}
	return nearest;
}
