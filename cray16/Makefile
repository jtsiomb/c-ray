obj = cray16.obj color.obj gfx.obj log.obj
bin = cray16.exe

CC = tcc

# 286 instructions, large memory model
CFLAGS = -1 -f287 -G

$(bin): $(obj)
	$(CC) $(obj)

.SUFFIXES: .c .obj

.c.obj:
	$(CC) $(CFLAGS) -c $<

clean:
	del *.obj
	del $(bin)
