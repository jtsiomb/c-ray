#ifndef COLOR_H_
#define COLOR_H_

int init_palette(const char *palfile);
int alloc_color(int r, int g, int b);
void lookup_color(int idx, int *r, int *g, int *b);

#endif	/* COLOR_H_ */
