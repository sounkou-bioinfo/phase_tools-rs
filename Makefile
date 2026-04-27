CARGO ?= cargo
PREFIX ?= $(HOME)/.local
BINDIR ?= $(PREFIX)/bin
INSTALL ?= install
TARGET ?=
TARGET_ARG := $(if $(TARGET),--target $(TARGET),)
RELEASE_BIN ?= target/release/phase_mnv_rs
STATIC_BIN ?= target/$(shell $(CARGO) -vV | sed -n 's/^host: //p')/release/phase_mnv_rs

.PHONY: release static-release install install-static clean test c c-test c-static byte-test

release:
	$(CARGO) build --release $(TARGET_ARG)

static-release:
	./scripts/build_rust_static.sh

install: release
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 0755 $(RELEASE_BIN) $(BINDIR)/phase_mnv_rs

install-static: static-release
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 0755 $(STATIC_BIN) $(BINDIR)/phase_mnv_rs

clean:
	$(CARGO) clean
	$(MAKE) -C c clean

test: release
	./tests/test_phase_mnv.sh $(RELEASE_BIN)

c:
	$(MAKE) -C c

c-test:
	$(MAKE) -C c test

c-static:
	./scripts/build_c_static.sh

byte-test: release c-test
	./test_byte_identical.sh
