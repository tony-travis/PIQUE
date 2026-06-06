#@(#)Makefile  2026-06-06  A.J.Travis and A.Douglas, last modified by A.J.Travis

#
# PIQUE - Parallel Identification Of QTL's using EMMAX
#

PLINK = plink2
PIQUE-INPUT = pique-input
PIQUE-RUN = pique-run
PHENO_GROUP = rice_phenotype_group.txt
PHENO = rice_phenotype.txt
PLOT_FILES = plot_files.txt
IPREFIX = sativas413
OPREFIX = sativas_GWAS
HPREFIX = sativas_GWAS_het
SNPHET_PREFIX  = sativas_GWAS_snphet
SNPZHET_PREFIX = sativas_GWAS_snpzhet
TPREFIX = sativas_GWAS_trans
GPREFIX = sativas_group_GWAS
GROUPS = ADMIX AUS IND TEJ TRJ
THREADS = 8
OPT = -d -v
TEST := $(dir $(realpath $(firstword $(MAKEFILE_LIST))))
PATH := $(realpath $(TEST)/../bin):$(PATH)

help:
	@echo "Targets:"
	@echo "\tinput"
	@echo "\trun"
	@echo "\tinput_het"
	@echo "\trun_het"
	@echo "\tinput_snphet"
	@echo "\trun_snphet"
	@echo "\tinput_snpzhet"
	@echo "\trun_snpzhet"
	@echo "\ttrans"
	@echo "\tgroup"
	@echo "\ttest_vcf"
	@echo "\tplot"
	@echo "\tall"

all:  input run input_het run_het input_snphet run_snphet input_snpzhet run_snpzhet trans group test_vcf plot

input: $(IPREFIX).ped
	rm -rf $(OPREFIX)
	$(PIQUE-INPUT) $(OPT) -i $(IPREFIX) -o $(OPREFIX) -p $(PHENO) -k -c -e 5
	touch $@

input_het: $(IPREFIX).ped
	rm -rf $(HPREFIX)
	$(PIQUE-INPUT) $(OPT) -i $(IPREFIX) -o $(HPREFIX) -p $(PHENO) -k -c -e 5 -z 0.999
	touch $@

run: input
	$(PIQUE-RUN) $(OPT) -i $(OPREFIX) -k IBS -c $(OPREFIX).covar -n $(THREADS)
	touch $@

run_het: input_het
	$(PIQUE-RUN) $(OPT) -i $(HPREFIX) -k IBS -c $(HPREFIX).covar -n $(THREADS)
	touch $@

input_snphet: $(IPREFIX).ped
	rm -rf $(SNPHET_PREFIX)
	$(PIQUE-INPUT) $(OPT) -i $(IPREFIX) -o $(SNPHET_PREFIX) -p $(PHENO) -k -c -e 5 -y 0.005
	touch $@

run_snphet: input_snphet
	$(PIQUE-RUN) $(OPT) -i $(SNPHET_PREFIX) -k IBS -c $(SNPHET_PREFIX).covar -n $(THREADS)
	touch $@

input_snpzhet: $(IPREFIX).ped
	rm -rf $(SNPZHET_PREFIX)
	$(PIQUE-INPUT) $(OPT) -i $(IPREFIX) -o $(SNPZHET_PREFIX) -p $(PHENO) -k -c -e 5 -y 0.005 -z 0.999
	touch $@

run_snpzhet: input_snpzhet
	$(PIQUE-RUN) $(OPT) -i $(SNPZHET_PREFIX) -k IBS -c $(SNPZHET_PREFIX).covar -n $(THREADS)
	touch $@

trans: $(IPREFIX).ped
	rm -rf $(TPREFIX)
	$(PIQUE-INPUT) $(OPT) -i $(IPREFIX) -o $(TPREFIX) -p $(PHENO) -k -c -e 5
	$(PIQUE-RUN) $(OPT) -i $(TPREFIX) -k IBS -c $(TPREFIX).covar -n $(THREADS) -t BC

group: $(IPREFIX).ped
	rm -rf $(GPREFIX)
	$(PIQUE-INPUT) $(OPT) -i $(IPREFIX) -o $(GPREFIX) -p $(PHENO_GROUP) -k -c -e 5
	$(PIQUE-RUN) $(OPT) -i $(GPREFIX) -k IBS -c $(GPREFIX).covar -n $(THREADS)
	for group in $(GROUPS); do \
	    prefix=$(GPREFIX)\_$$group; \
	    $(PIQUE-RUN) $(OPT) -i $$prefix -k IBS -n $(THREADS); \
	done

$(IPREFIX).ped $(IPREFIX).map:
	tar xf RiceDiversity_44K_Genotypes_PLINK.tar.xz
	mv -i ./RiceDiversity_44K_Genotypes_PLINK/sativas* .
	rmdir RiceDiversity_44K_Genotypes_PLINK

test_vcf: test_vcf.vcf
	rm -rf $@_GWAS
	$(PIQUE-INPUT) $(OPT) -i $@ -f vcf -o $@_GWAS -p $(PHENO) -k -c -e 5
	$(PIQUE-RUN) -i $@_GWAS -k IBS -c $@_GWAS.covar -n $(THREADS)

test_vcf.vcf: $(IPREFIX).ped $(IPREFIX).map
	$(PLINK) --pedmap $(IPREFIX) --recode vcf --out $(basename $@)

plot: $(PLOT_FILES) run
	GWAS_manhattanplots -b 0.1 -q -y 12 -i $<

clean:
	rm -rf RiceDiversity_44K_Genotypes_PLINK sativas413*
	rm -rf __MACOSX
	rm -rf $(OPREFIX)* $(HPREFIX)* $(SNPHET_PREFIX)* $(SNPZHET_PREFIX)* $(TPREFIX)* $(GPREFIX)*
	rm -rf test_vcf*
	rm -f input run input_het run_het input_snphet run_snphet input_snpzhet run_snpzhet
	rm -f *.bh *.png
	rm -f .pversion *.log

clobber: clean
	rm -f $(IPREFIX)_recode12.*
