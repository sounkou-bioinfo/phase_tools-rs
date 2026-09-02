CARGO ?= cargo
LAKE ?= lake
PREFIX ?= $(HOME)/.local
BINDIR ?= $(PREFIX)/bin
INSTALL ?= install
BIN ?= target/release/phase-tools

.PHONY: all fmt test lint proof release smoke install clean

all: test proof

fmt:
	$(CARGO) fmt --check

test:
	$(CARGO) test --all-targets

lint:
	$(CARGO) clippy --all-targets -- -D warnings

proof:
	$(LAKE) build

release:
	$(CARGO) build --release

smoke: release
	$(BIN) targets --assay wes --validated-enrichment
	rm -f /tmp/phase-tools-hba.cert
	$(BIN) hba \
		--assay wgs \
		--evidence examples/hba/evidence.synthetic.tsv \
		--hypotheses examples/hba/hypotheses.synthetic.tsv \
		--min-margin 10 \
		--certificate /tmp/phase-tools-hba.cert
	$(BIN) verify --certificate /tmp/phase-tools-hba.cert

install: release
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 0755 $(BIN) $(BINDIR)/phase-tools

clean:
	$(CARGO) clean
	$(LAKE) clean
