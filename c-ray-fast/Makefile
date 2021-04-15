PREFIX = /usr/local

obj = c-ray-fast.o
bin = c-ray-fast

CC = gcc
CFLAGS = -pedantic -Wall -g -O3 -ffast-math
LDFLAGS = -lpthread -lm

$(bin): $(obj)
	$(CC) -o $@ $(obj) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) $(bin)

.PHONY: install
install:
	cp $(bin) $(DESTDIR)$(PREFIX)/$(bin)

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/$(bin)
