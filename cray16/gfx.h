#ifndef GFX_H_
#define GFX_H_

extern unsigned char far *vidmem;

int set_video_mode(int xsz, int ysz, int bpp);

void save_vga(void);
void restore_vga(void);

void put_pixel_modex(int x, int y, int col);
void set_palette(int idx, int r, int g, int b);

#define set_plane(p)	set_plane_mask(1 << (p))
void set_plane_mask(unsigned char pmask);

#endif	/* GFX_H_ */
