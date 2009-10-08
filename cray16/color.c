#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <errno.h>
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


int init_palette(const char *palfile)
{
	FILE *fp;
	int i, num_col, line_num;
	char line[64];

	if(!palfile) {
		logfoo("setting up fixed palette\n");
		num_col = 256;

		for(i=0; i<num_col; i++) {
			int r = i & 0xe0;
			int g = (i & 0x1c) << 3;
			int b = (i & 3) << 6;

			palette[i].r = r;
			palette[i].g = g;
			palette[i].b = b;

			set_palette(i, r, g, b);
		}
		
		free_col_idx = 256;
		return 0;
	}

	if(!(fp = fopen(palfile, "r"))) {
		logfoo("failed to open palette: %s: %s\n", palfile, strerror(errno));
		return -1;
	}

	num_col = 0;
	line_num = 1;
	while(fgets(line, sizeof line, fp)) {
		int red, green, blue;

		if(sscanf(line, "%d %d %d", &red, &green, &blue) != 3) {
			logfoo("palette file error at line: %d: %s\n", line_num, line);
			fclose(fp);
			return -1;
		}

		if(num_col >= 256) {
			num_col++;
			continue;
		}

		if(red > 255) red = 255;
		if(green > 255) green = 255;
		if(blue > 255) blue = 255;

		palette[num_col].r = red;
		palette[num_col].g = green;
		palette[num_col].b = blue;

		logfoo("pal[%d] = %3d %3d %3d\n", red, green, blue);

		set_palette(num_col++, red, green, blue);
	}

	fclose(fp);

	if(num_col > 256) {
		logfoo("palette file %s contained %d colors, using the first 256\n", palfile, num_col);
		num_col = 256;
	} else {
		logfoo("loaded %d colors from palette file %s\n", num_col, palfile);
	}

	free_col_idx = 256;
	return 0;
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
