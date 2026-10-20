CC ?= gcc
CFLAGS ?= -O2 -Wall -Wextra
LDLIBS = -luhd -lfftw3f -lm

all: nr_scan

nr_scan: nr_scan.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f nr_scan *.csv *.png *.gp *.wav

.PHONY: all clean
