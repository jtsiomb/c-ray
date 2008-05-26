obj = c-ray-gpu.o sdr.o scene.o
bin = c-ray-gpu

CC = gcc
CFLAGS = -pedantic -Wall -g `pkg-config --cflags vmath`
LDFLAGS = -lGL -lglut -lm `pkg-config --libs vmath` -lGLEW

$(bin): $(obj)
	$(CC) -o $@ $(obj) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) $(bin)
