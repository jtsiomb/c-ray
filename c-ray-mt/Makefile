obj = c-ray-mt.o
bin = c-ray-mt

CC = gcc
CFLAGS = -pedantic -Wall -O3 -g -ffast-math -DBUILD_COOK_TORRANCE
LDFLAGS = -L/usr/local/lib -lm -lpthread -lkdtree

$(bin): $(obj)
	$(CC) -o $@ $(obj) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) $(bin)

.PHONY: install
install:
	cp $(bin) /usr/local/bin/$(bin)

.PHONY: uninstall
uninstall:
	rm -f /usr/local/bin/$(bin)
