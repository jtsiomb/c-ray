obj := c-ray.o
bin := c-ray

opt := -O3 -ffast-math -mfpmath=sse
CFLAGS := $(opt) -std=gnu89 -pedantic -Wall `sdl-config --cflags`

$(bin): $(obj)
	$(CC) -o $@ $(obj) `sdl-config --libs`

.PHONY: clean
clean:
	$(RM) $(obj) $(bin)
