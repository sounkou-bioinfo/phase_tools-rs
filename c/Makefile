CC ?= gcc
CFLAGS ?= -O3 -Wall -Wextra -std=c11
PREFIX ?= $(HOME)/.local
BINDIR ?= $(PREFIX)/bin
INSTALL ?= install

# Prefer pkg-config when htslib development files are installed. Override for
# static/bundled builds, e.g.:
#   make HTS_CPPFLAGS='-I/path/to/htslib' HTS_LIBS='/path/to/libhts.a -lz -lm -lpthread'
HTS_CPPFLAGS_PKG := $(shell pkg-config --cflags htslib 2>/dev/null)
HTS_LIBS_PKG     := $(shell pkg-config --libs htslib 2>/dev/null)

ifeq ($(strip $(HTS_CPPFLAGS_PKG)),)
  ifneq ($(wildcard /usr/lib/R/site-library/Rhtslib/include/htslib/vcf.h),)
    HTS_CPPFLAGS ?= -I/usr/lib/R/site-library/Rhtslib/include
  else ifneq ($(wildcard /usr/lib/R/site-library/vcfppR/include/htslib/vcf.h),)
    HTS_CPPFLAGS ?= -I/usr/lib/R/site-library/vcfppR/include
  else
    HTS_CPPFLAGS ?=
  endif
else
  HTS_CPPFLAGS ?= $(HTS_CPPFLAGS_PKG)
endif

ifeq ($(strip $(HTS_LIBS_PKG)),)
  ifneq ($(wildcard /usr/lib/x86_64-linux-gnu/libhts.so.3),)
    HTS_LIBS ?= -Wl,-l:libhts.so.3
  else ifneq ($(wildcard /lib/x86_64-linux-gnu/libhts.so.3),)
    HTS_LIBS ?= -Wl,-l:libhts.so.3
  else
    HTS_LIBS ?= -lhts
  endif
else
  HTS_LIBS ?= $(HTS_LIBS_PKG)
endif

.PHONY: all clean test install uninstall

all: phase_mnv

phase_mnv: phase_mnv.c
	$(CC) $(CFLAGS) $(HTS_CPPFLAGS) $< $(HTS_LIBS) -o $@

clean:
	rm -f phase_mnv

install: phase_mnv
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 0755 phase_mnv $(BINDIR)/phase_mnv

uninstall:
	rm -f $(BINDIR)/phase_mnv

test: phase_mnv
	../tests/test_phase_mnv.sh ./phase_mnv
