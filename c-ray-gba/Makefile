obj := crt0.o c-ray-gba.o
bin := c-ray-gba.gba

CC := arm-agb-elf-gcc
AS := arm-agb-elf-as
OBJCOPY := arm-agb-elf-objcopy

opt := -O3 -fomit-frame-pointer -ffast-math -marm -mcpu=arm7tdmi
libdir := /home/nuclear/code/gba/libgba

CFLAGS := $(opt) -std=c99 -pedantic -Wall -I$(libdir)/src

$(bin): $(obj) $(libdir)/libgba.a
	$(CC) -o /tmp/out.elf $(obj) $(libdir)/libgba.a -Tlnkscript -nostartfiles -lm
	$(OBJCOPY) -O binary /tmp/out.elf /tmp/out.flat
	gba-header -i /tmp/out.flat -o $@ -fix -title c-ray -gcode 1337 -mcode 12
	chmod a+r $@


-include $(obj:.o=.d)

%.d: %.c
	@rm -f $@; $(CC) -MM $(CFLAGS) $< > $@

.PHONY: clean
clean:
	rm -f $(obj) $(obj:.o=.d)
