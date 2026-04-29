CARGO ?= cargo
PREFIX ?= $(HOME)/.local
BINDIR ?= $(PREFIX)/bin
INSTALL ?= install
TARGET ?=
TARGET_ARG := $(if $(TARGET),--target $(TARGET),)
RELEASE_BIN ?= target/release/phase_mnv_rs
PHASE_COMPARE_BIN ?= target/release/phase_compare
STATIC_BIN ?= target/$(shell $(CARGO) -vV | sed -n 's/^host: //p')/release/phase_mnv_rs
STATIC_PHASE_COMPARE_BIN ?= target/$(shell $(CARGO) -vV | sed -n 's/^host: //p')/release/phase_compare

.PHONY: release static-release install install-static clean test negative-test c c-test c-negative-test c-static byte-test compare-whatshap-phase readme readme-external-example check-readme

release:
	$(CARGO) build --release $(TARGET_ARG)

static-release:
	./scripts/build_rust_static.sh

install: release
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 0755 $(RELEASE_BIN) $(BINDIR)/phase_mnv_rs
	$(INSTALL) -m 0755 $(PHASE_COMPARE_BIN) $(BINDIR)/phase_compare

install-static: static-release
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 0755 $(STATIC_BIN) $(BINDIR)/phase_mnv_rs
	$(INSTALL) -m 0755 $(STATIC_PHASE_COMPARE_BIN) $(BINDIR)/phase_compare

clean:
	$(CARGO) clean
	$(MAKE) -C c clean

test: release
	./tests/test_phase_mnv.sh $(RELEASE_BIN)
	./tests/test_output_formats.sh $(RELEASE_BIN)
	./tests/test_bcftools_norm.sh $(RELEASE_BIN)
	./tests/test_all_sites.sh $(RELEASE_BIN)
	./tests/test_bam_phase.sh $(RELEASE_BIN)
	./tests/test_phase_compare.sh $(PHASE_COMPARE_BIN)
	./tests/test_negative.sh $(RELEASE_BIN)

negative-test: release
	./tests/test_negative.sh $(RELEASE_BIN)

c:
	$(MAKE) -C c

c-test:
	$(MAKE) -C c test

c-negative-test: c
	./tests/test_negative.sh c/phase_mnv

c-static:
	./scripts/build_c_static.sh

byte-test: release c-test
	./test_byte_identical.sh

compare-whatshap-phase: release
	./scripts/compare_whatshap_phase.sh

readme: release c
	Rscript -e 'invisible(suppressWarnings(knitr::knit("README.Rmd", "README.md", quiet = TRUE)))'
	perl -0pi -e 's/\A(# phase_mnv_rs)\n{3,}/$$1\n\n/' README.md

readme-external-example:
	PHASE_MNV_RUN_EXTERNAL=1 $(MAKE) readme

check-readme:
	env -u PHASE_MNV_RUN_EXTERNAL -u PHASE_MNV_EXAMPLE_VCF -u PHASE_MNV_EXAMPLE_REF -u PHASE_MNV_EXAMPLE_SAMPLE $(MAKE) readme
	git diff --exit-code README.md
