obj := c-ray.o
bin := c-ray

opt := -O3 -ffast-math
CC := gcc
CFLAGS := $(opt) -std=c99 -pedantic -Wall `sdl-config --cflags`

$(bin): $(obj)
	$(CC) -o $@ $(obj) -lm `sdl-config --libs`

.PHONY: clean
clean:
	$(RM) $(obj) $(bin)
