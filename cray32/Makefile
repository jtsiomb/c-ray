obj = cray32.obj gfx.obj vbe.obj dpmi.obj
bin = cray32.exe

CC = wcc386
CFLAGS = -fp5 -ot
LD = wlink

$(bin): $(obj)
	$(LD) name $@ file { $(obj) } $(LDFLAGS)

.c.obj: .autodepend
	$(CC) $(CFLAGS) $<

clean: .symbolic
	del *.obj
	del $(bin)
