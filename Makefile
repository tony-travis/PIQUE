#@(#)Makefile  2021-05-22  A.J.Travis and A.Douglas

#
# pique: parallel identification of QTL's using EMMAX
#
# A software pipeline for performing GWAS (Genome Wide Association Studies)
#

# installation directory
DIR = /usr/local/pique

EMMAX = /usr/bin/emmax
PLINK = /usr/bin/p-link
EIGENSTRAT = /usr/lib/eigensoft/smartpca
R = /usr/bin/R
FORECAST = /usr/lib/R/site-library/forecast
PARALLEL = /usr/share/perl5/Parallel
READONLY = /usr/share/perl5/Readonly.pm
BIN = bin/GWAS_manhattanplots bin/pique-input bin/pique-run
ETC = etc/profile.d/EIGENSOFT.sh etc/profile.d/pique.sh

TARGETS = $(EMMAX) $(PLINK) $(EIGENSTRAT) $(R) $(FORECAST) $(READONLY) $(PARALLEL) $(BIN) $(ETC)

help:
	@echo 'Type "sudo make install" to install "pique"'

all: $(TARGETS)

# install pique
install: $(TARGETS)
	install -d $(DIR)/bin $(DIR)/doc
	install -C -o root -g root bin/* $(DIR)/bin
	install -C -o root -g root etc/profile.d/* /etc/profile.d

# perl scripts
#pique-input: bin/pique-input
#pique-run: bin/pique-run

# install "EMMAX"
$(EMMAX):
	apt install emmax

emmax-beta-07Mar2010.tar.gz:
	wget http://genetics.cs.ucla.edu/emmax/$@

# install PLINK
$(PLINK):
	apt install plink

plink_linux_x86_64.zip:
	wget https://www.cog-genomics.org/static/bin/plink161202/$@

# install "EIGENSOFT"
#$(EIGENSTRAT): EIG5.0.2.tar.gz
#	tar xvf $<
#	install -d /usr/local/EIGENSOFT/bin
#	install -C -o root -g root EIG5.0.2/bin/* /usr/local/EIGENSOFT/bin
#
#EIG5.0.2.tar.gz:
#	wget http://cdn1.sph.harvard.edu/wp-content/uploads/sites/181/2014/05/$@
$(EIGENSTRAT):
	apt install eigensoft

# install R
$(R):
	echo "Please install R"...

# install CRAN forecast library
$(FORECAST):
	apt install r-cran-forecast

# install Parallel::ForkManager
$(PARALLEL):
	apt install libparallel-forkmanager-perl

# install Perl Readonly module
$(READONLY):
	apt install libreadonly-perl

%.pdf: %.odt
	lowriter --headless --convert-to pdf $< --outdir $$(dirname $<)

clean:
	rm -f *.log

clobber: clean
	rm -rf EIG5.0.2*  emmax-beta-07Mar2010*
