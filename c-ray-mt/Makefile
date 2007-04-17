obj = c-ray-mt.o
bin = c-ray-mt

CC = gcc
CFLAGS = -O3 -ffast-math -DBUILD_COOK_TORRANCE

$(bin): $(obj)
	$(CC) -o $@ $(obj) -lm -lpthread

.PHONY: clean
clean:
	rm -f $(obj) $(bin)

.PHONY: install
install:
	cp $(bin) /usr/local/bin/$(bin)

.PHONY: uninstall
uninstall:
	rm -f /usr/local/bin/$(bin)
