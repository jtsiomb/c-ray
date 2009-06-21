#include <stdio.h>
#include <stdarg.h>
#include "log.h"

void logfoo(const char *fmt, ...)
{
	static FILE *fp;
	va_list ap;

	if(!fp) {
		if(!(fp = fopen("logfile.log", "wb"))) {
			fprintf(stderr, "failed to open log file, logging to stderr\n");
			fp = stderr;
		}
	}

	va_start(ap, fmt);
	vfprintf(fp, fmt, ap);
	va_end(ap);

	fflush(fp);
}
