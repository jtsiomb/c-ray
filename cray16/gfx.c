#include <stdio.h>
#include <dos.h>
#include "gfx.h"

#define ATTR_ADDR		0x3c0
#define ATTR_DATA_WR	0x3c0
#define ATTR_DATA_RD	0x3c1

#define GFX_ADDR		0x3ce
#define GFX_DATA		0x3cf

#define SEQ_ADDR		0x3c4
#define SEQ_DATA		0x3c5

#define CRTC_ADDR		0x3d4
#define CRTC_DATA		0x3d5

#define DAC_STATE		0x3c7
#define DAC_ADDR_RD		0x3c7
#define DAC_ADDR_WR		0x3c8
#define DAC_DATA		0x3c9

#define MISC_OUTP_WR	0x3c2
#define MISC_OUTP_RD	0x3cc

#define INP_STATUS		0x3da

enum {
	SEQ_RESET,
	SEQ_CLOCK_MODE,
	SEQ_MAP_MASK,
	SEQ_CHAR_MAP,
	SEQ_MEM_MODE,

	NUM_SEQ_REGS
};

enum {
	CRTC_HTOTAL,
	CRTC_HDISP_END,
	CRTC_HBLANK_START,
	CRTC_HBLANK_END,
	CRTC_HRETRACE_START,
	CRTC_HRETRACE_END,
	CRTC_VTOTAL,
	CRTC_OVERFLOW,
	CRTC_PRESET_ROW_SCAN,
	CRTC_MAX_SCANLINE,
	CRTC_CURSOR_START,
	CRTC_CURSOR_END,
	CRTC_START_ADDR_HIGH,
	CRTC_START_ADDR_LOW,
	CRTC_CURSOR_LOC_HIGH,
	CRTC_CURSOR_LOC_LOW,
	CRTC_VRETRACE_START,
	CRTC_VRETRACE_END,
	CRTC_VDISP_END,
	CRTC_OFFSET,
	CRTC_UNDERLINE_LOC,
	CRTC_VBLANK_START,
	CRTC_VBLANK_END,
	CRTC_MODE_CTL,
	CRTC_LINE_CMP,

	NUM_CRTC_REGS
};

enum {
	GFX_SR,
	GFX_ENABLE_SR,
	GFX_COLOR_CMP,
	GFX_DATA_ROT,
	GFX_READ_MAP,
	GFX_MODE,
	GFX_MISC,
	GFX_COLOR_DONT_CARE,
	GFX_BIT_MASK,

	NUM_GFX_REGS
};

enum {
	ATTR_PALETTE_BASE,	/* 0 - 15 */
	ATTR_MODE_CTL = 16,
	ATTR_OVERSCAN_COLOR,
	ATTR_COLOR_PLANE,
	ATTR_HPAN,
	ATTR_COLOR_SEL,

	NUM_ATTR_REGS
};

static struct vga_regs {
	unsigned char attr[NUM_ATTR_REGS];
	unsigned char gfx[NUM_GFX_REGS];
	unsigned char seq[NUM_SEQ_REGS];
	unsigned char crtc[NUM_CRTC_REGS];
	unsigned char misc_outp;
	unsigned char dac[256 * 3];
} vga;

unsigned char far *vidmem = MK_FP(0xa000, 0);

static void set_modex(void);

int set_video_mode(int xsz, int ysz, int bpp)
{
	if(xsz != 320 || bpp != 8) {
		fprintf(stderr, "currently only 320-wide 8bpp modes supported\n");
		return -1;
	}

	switch(ysz) {
	case 200:
		save_vga();
		asm {
			mov ax, 0x13
			int 0x10
		}
		break;

	case 240:
		save_vga();
		set_modex();
		break;

	default:
		fprintf(stderr, "currently only 200 and 240 heights supported\n");
		return -1;
	}
	return 0;
}

void save_vga(void)
{
	int i, j;
	unsigned char val;

	/* save attribute registers */
	for(i=0; i<NUM_ATTR_REGS; i++) {
		asm {
			mov dx, INP_STATUS
			in al, dx

			mov dx, ATTR_ADDR
			mov ax, i
			out dx, ax
			mov dx, ATTR_DATA_RD
			in al, dx
			mov val, al
		}
		vga.attr[i] = val;
	}

	/* save graphics registers */
	for(i=0; i<NUM_GFX_REGS; i++) {
		asm {
			mov dx, GFX_ADDR
			mov ax, i
			out dx, al
			inc dx
			in al, dx
			mov val, al
		}
		vga.gfx[i] = val;
	}

	/* save sequencer registers */
	for(i=0; i<NUM_SEQ_REGS; i++) {
		asm {
			mov dx, SEQ_ADDR
			mov ax, i
			out dx, al
			inc dx
			in al, dx
			mov val, al
		}
		vga.seq[i] = val;
	}

	/* save crtc registers */
	for(i=0; i<NUM_CRTC_REGS; i++) {
		asm {
			mov dx, CRTC_ADDR
			mov ax, i
			out dx, al
			inc dx
			in al, dx
			mov val, al
		}
		vga.crtc[i] = val;
	}

	/* save the miscellaneous output register */
	asm {
		mov dx, MISC_OUTP_RD
		in al, dx
		mov val, al
	}
	vga.misc_outp = val;

	/* save palette */
	for(i=0; i<256; i++) {
		asm {
			mov dx, DAC_ADDR_RD
			mov ax, i
			out dx, al
		}
		for(j=0; j<3; j++) {
			asm {
				mov dx, DAC_DATA
				in al, dx
				mov val, al
			}
			vga.dac[i * 3 + j] = val;
		}
	}
}

void restore_vga(void)
{
	int i, j;
	unsigned char val;

	/* restore attribute registers */
	for(i=0; i<NUM_ATTR_REGS; i++) {
		val = vga.attr[i];
		asm {
			mov dx, INP_STATUS
			in al, dx

			mov dx, ATTR_ADDR
			mov ax, i
			out dx, al
			mov al, val
			out dx, al
		}
	}

	/* restore graphics registers */
	for(i=0; i<NUM_GFX_REGS; i++) {
		val = vga.gfx[i];
		asm {
			mov dx, GFX_ADDR
			mov ax, i
			out dx, al
			inc dx
			mov al, val
			out dx, al
		}
	}

	/* restore sequencer registers */
	for(i=0; i<NUM_SEQ_REGS; i++) {
		val = vga.seq[i];
		asm {
			mov dx, SEQ_ADDR
			mov ax, i
			out dx, al
			inc dx
			mov al, val
			out dx, al
		}
	}

	/* restore crtc registers */
	for(i=0; i<NUM_CRTC_REGS; i++) {
		val = vga.crtc[i];
		asm {
			mov dx, CRTC_ADDR
			mov ax, i
			out dx, al
			inc dx
			mov al, val
			out dx, al
		}
	}

	/* restore the miscellaneous output register */
	val = vga.misc_outp;
	asm {
		mov dx, MISC_OUTP_WR
		mov al, val
		out dx, al
	}

	/* restore palette */
	for(i=0; i<256; i++) {
		asm {
			mov dx, DAC_ADDR_WR
			mov ax, i
			out dx, al
		}
		for(j=0; j<3; j++) {
			val = vga.dac[i * 3 + j];
			asm {
				mov dx, DAC_DATA
				mov al, val
				out dx, al
			}
		}
	}
}


void put_pixel_modex(int x, int y, int col)
{
	unsigned char far *ptr = vidmem + y * 80 + (x >> 2);

	/* select plane */
	unsigned char plane = x & 3;
	unsigned short bits = (1 << plane) | (SEQ_MAP_MASK << 8);
	asm {
		mov dx, SEQ_ADDR
		mov ax, bits
		out dx, ax
	}

	*ptr = (unsigned char)col;
}

void set_palette(int idx, int r, int g, int b)
{
	asm {
		mov dx, DAC_ADDR_WR
		mov ax, idx
		out dx, al
		inc dx
		mov ax, r
		shr ax, 2
		out dx, al
		mov ax, g
		shr ax, 2
		out dx, al
		mov ax, b
		shr ax, 2
		out dx, al
	}
}

void set_plane_mask(unsigned char pmask)
{
	unsigned short bits = (pmask << 8) | SEQ_MAP_MASK;
	asm {
		mov dx, SEQ_ADDR
		mov ax, bits
		out dx, ax
	}
}

static void set_modex(void)
{
	static unsigned short crt_params[] = {
		0x0d06,		/* vertical total */
		0x3e07,		/* overflow */
		0x4109,		/* cell height (2 to double scan) */
		0xea10,		/* vsync start */
		0xac11,		/* vsync end and protect cr0-cr7 */
		0xdf12,		/* vertical displayed */
		0x0014,		/* turn off dword mode */
		0xe715,		/* vblank start */
		0x0616,		/* vblank end */
		0xe317		/* turn on byte mode */
	};

	unsigned int crt_params_offs = (unsigned int)crt_params;
	unsigned int crt_params_len = sizeof crt_params;

	asm {
		/* start by entering mode 13h */
		mov ax, 13h
		int 10h

		/* disable chaining */
		mov dx, SEQ_ADDR
		mov ax, 0604h
		out dx, ax

		/* syncronous reset while setting misc */
		mov ax, 0100h
		out dx, ax

		/* select 25mhz dot clock & 60hz refresh */
		mov dx, MISC_OUTP_WR
		mov al, 0e3h
		out dx, al

		/* restart sequencer */
		mov dx, SEQ_ADDR
		mov ax, 0300h
		out dx, ax

		mov dx, CRTC_ADDR	/* reprogram CRTC */
		mov al, 11h			/* vsync end reg contains reg write protect bit */
		out dx, al
		inc dx				/* CRTC data register */
		in al, dx			/* get vsync end register */
		and al, 7fh			/* remove write protect */
		out dx, al
		dec dx				/* CRTC index */
		cld
		mov si, crt_params_offs
		mov cx, crt_params_len
	}
paramloop:
	asm {
		lodsw				/* get next CRT index/data pair */
		out dx, ax			/* set the next CRT index/data pair */
		loop paramloop

		mov dx, SEQ_ADDR
		mov ax, 0f02h
		out dx, ax			/* enable writes to all 4 planes */
		mov ax, 0a000h
		mov es, ax
		xor di, di
		xor ax, ax
		mov cx, 8000h
		rep stosw
	}
}
