obj = c-ray-rl.o
bin = c-ray-rl

CC = gcc
CFLAGS = -O3 -ffast-math

$(bin): $(obj)
	$(CC) -o $@ $(obj) -lm

.PHONY: clean
clean:
	rm -f $(obj) $(bin)
