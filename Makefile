#@(#)Makefile  2026-06-06  A.J.Travis and A.Douglas

#
# pique: parallel identification of QTL's using EMMAX
#
# A software pipeline for performing GWAS (Genome Wide Association Studies)
#

# installation directory
DIR = /usr/local/pique

BIN = bin/GWAS_manhattanplots bin/pique-input bin/pique-run
ETC = etc/profile.d/EIGENSOFT.sh etc/profile.d/pique.sh

# packages installed via apt (R handled separately)
PACKAGES = emmax plink2 eigensoft libparallel-forkmanager-perl libreadonly-perl r-base-core r-cran-forecast

.PHONY: help all check-deps install-deps install clean clobber

help:
	@echo 'Type "make check-deps"        to check dependencies without installing'
	@echo 'Type "sudo make install-deps" to install dependencies only'
	@echo 'Type "sudo make install"      to install dependencies and pique'

all: install-deps
	install -d $(DIR)/bin $(DIR)/doc
	install -C -o root -g root bin/* $(DIR)/bin
	install -C -o root -g root etc/profile.d/* /etc/profile.d

# check dependencies without installing
check-deps:
	@echo "==> Checking dependencies..."
	@missing=""; \
	for pkg in $(PACKAGES); do \
		dpkg -s $$pkg >/dev/null 2>&1 || missing="$$missing $$pkg"; \
	done; \
	if [ -n "$$missing" ]; then \
		echo "Missing packages:$$missing"; \
	else \
		echo "All apt dependencies are installed."; \
	fi

# install dependencies
install-deps:
	@echo "==> Checking and installing apt dependencies..."
	@missing=""; \
	for pkg in $(PACKAGES); do \
		dpkg -s $$pkg >/dev/null 2>&1 || missing="$$missing $$pkg"; \
	done; \
	if [ -n "$$missing" ]; then \
		echo "Installing:$$missing"; \
		apt-get install -y $$missing; \
	else \
		echo "All apt dependencies already installed."; \
	fi

# install pique
install: install-deps
	@echo "==> Installing PIQUE to $(DIR)"
	install -d $(DIR)/bin $(DIR)/doc
	install -C -o root -g root bin/* $(DIR)/bin
	install -C -o root -g root etc/profile.d/* /etc/profile.d
	@echo "==> Done. Run: source /etc/profile.d/pique.sh"

%.pdf: %.odt
	lowriter --headless --convert-to pdf $< --outdir $$(dirname $<)

clean:
	rm -f *.log

clobber: clean
