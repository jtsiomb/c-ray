obj := c-ray-unix.o
bin := c-ray-unix

opt := -O3 -ffast-math -mfpmath=sse
CFLAGS := $(opt) -std=gnu89 -pedantic -Wall -I/usr/X11R6/include

$(bin): $(obj)
	$(CC) -o $@ $(obj) -L/usr/X11R6/lib -lX11 -lm

.PHONY: clean
clean:
	$(RM) $(obj) $(bin)
