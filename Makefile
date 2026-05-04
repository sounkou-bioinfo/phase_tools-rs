CARGO ?= cargo
PREFIX ?= $(HOME)/.local
BINDIR ?= $(PREFIX)/bin
INSTALL ?= install
TARGET ?=
TARGET_ARG := $(if $(TARGET),--target $(TARGET),)
RELEASE_BIN ?= target/release/phase_mnv_rs
PHASE_COMPARE_BIN ?= target/release/phase_compare
UNPHASE_BIN ?= target/release/unphase_vcf
FERMI_LITE_BIN ?= target/release/fermi_lite_assemble
BAM_ERROR_MODEL_BIN ?= target/release/bam_error_model
PHASE_ADJUDICATE_BIN ?= target/release/phase_adjudicate
BAM_CONTAMINATION_BIN ?= target/release/bam_contamination
STATIC_BIN ?= target/$(shell $(CARGO) -vV | sed -n 's/^host: //p')/release/phase_mnv_rs
STATIC_PHASE_COMPARE_BIN ?= target/$(shell $(CARGO) -vV | sed -n 's/^host: //p')/release/phase_compare
STATIC_UNPHASE_BIN ?= target/$(shell $(CARGO) -vV | sed -n 's/^host: //p')/release/unphase_vcf
STATIC_FERMI_LITE_BIN ?= target/$(shell $(CARGO) -vV | sed -n 's/^host: //p')/release/fermi_lite_assemble
STATIC_BAM_ERROR_MODEL_BIN ?= target/$(shell $(CARGO) -vV | sed -n 's/^host: //p')/release/bam_error_model
STATIC_PHASE_ADJUDICATE_BIN ?= target/$(shell $(CARGO) -vV | sed -n 's/^host: //p')/release/phase_adjudicate
STATIC_BAM_CONTAMINATION_BIN ?= target/$(shell $(CARGO) -vV | sed -n 's/^host: //p')/release/bam_contamination

.PHONY: release static-release install install-static clean test negative-test compare-whatshap-phase readme readme-external-example check-readme

release:
	$(CARGO) build --release $(TARGET_ARG)

static-release:
	./scripts/build_rust_static.sh

install: release
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 0755 $(RELEASE_BIN) $(BINDIR)/phase_mnv_rs
	$(INSTALL) -m 0755 $(PHASE_COMPARE_BIN) $(BINDIR)/phase_compare
	$(INSTALL) -m 0755 $(UNPHASE_BIN) $(BINDIR)/unphase_vcf
	$(INSTALL) -m 0755 $(FERMI_LITE_BIN) $(BINDIR)/fermi_lite_assemble
	$(INSTALL) -m 0755 $(BAM_ERROR_MODEL_BIN) $(BINDIR)/bam_error_model
	$(INSTALL) -m 0755 $(PHASE_ADJUDICATE_BIN) $(BINDIR)/phase_adjudicate
	$(INSTALL) -m 0755 $(BAM_CONTAMINATION_BIN) $(BINDIR)/bam_contamination

install-static: static-release
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 0755 $(STATIC_BIN) $(BINDIR)/phase_mnv_rs
	$(INSTALL) -m 0755 $(STATIC_PHASE_COMPARE_BIN) $(BINDIR)/phase_compare
	$(INSTALL) -m 0755 $(STATIC_UNPHASE_BIN) $(BINDIR)/unphase_vcf
	$(INSTALL) -m 0755 $(STATIC_FERMI_LITE_BIN) $(BINDIR)/fermi_lite_assemble
	$(INSTALL) -m 0755 $(STATIC_BAM_ERROR_MODEL_BIN) $(BINDIR)/bam_error_model
	$(INSTALL) -m 0755 $(STATIC_PHASE_ADJUDICATE_BIN) $(BINDIR)/phase_adjudicate
	$(INSTALL) -m 0755 $(STATIC_BAM_CONTAMINATION_BIN) $(BINDIR)/bam_contamination

clean:
	$(CARGO) clean

test: release
	./tests/test_unphase_vcf.sh $(UNPHASE_BIN)
	./tests/test_phase_mnv.sh $(RELEASE_BIN)
	./tests/test_output_formats.sh $(RELEASE_BIN)
	./tests/test_bcftools_norm.sh $(RELEASE_BIN)
	./tests/test_all_sites.sh $(RELEASE_BIN)
	./tests/test_bam_phase.sh $(RELEASE_BIN)
	./tests/test_phase_compare.sh $(PHASE_COMPARE_BIN)
	./tests/test_fermi_lite.sh $(FERMI_LITE_BIN)
	./tests/test_bam_error_model.sh $(BAM_ERROR_MODEL_BIN)
	./tests/test_phase_adjudicate.sh $(PHASE_ADJUDICATE_BIN)
	./tests/test_bam_contamination.sh $(BAM_CONTAMINATION_BIN)
	./tests/test_negative.sh $(RELEASE_BIN)

negative-test: release
	./tests/test_negative.sh $(RELEASE_BIN)

compare-whatshap-phase: release
	./scripts/compare_whatshap_phase.sh

readme: release
	Rscript -e 'invisible(suppressWarnings(knitr::knit("README.Rmd", "README.md", quiet = TRUE)))'
	Rscript -e 'invisible(suppressWarnings(knitr::knit("docs/cli.Rmd", "docs/cli.md", quiet = TRUE)))'
	perl -0pi -e 's/\A(# phase_tools-rs)\n{3,}/$$1\n\n/' README.md
	perl -0pi -e 's/\A(# phase_tools-rs CLI help)\n{3,}/$$1\n\n/' docs/cli.md

readme-external-example:
	PHASE_MNV_RUN_EXTERNAL=1 $(MAKE) readme

check-readme:
	env -u PHASE_MNV_RUN_EXTERNAL -u PHASE_MNV_EXAMPLE_VCF -u PHASE_MNV_EXAMPLE_REF -u PHASE_MNV_EXAMPLE_SAMPLE $(MAKE) readme
	git diff --exit-code README.md docs/cli.md
