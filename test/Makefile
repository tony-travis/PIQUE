#@(#)Makefile  2022-01-27  A.J.Travis and A.Douglas

#
# PIQUE - Parallel Identification Of QTL's using EMMAX
#

WGET = wget --no-check-certificate 
PLINK = plink1.9
PIQUE-INPUT = pique-input
PIQUE-RUN = pique-run
PHENO_GROUP = rice_phenotype_group.txt
PHENO = rice_phenotype.txt
PLOT_FILES = plot_files.txt
IPREFIX = sativas413
OPREFIX = sativas_GWAS
TPREFIX = sativas_GWAS_trans
GPREFIX = sativas_group_GWAS
GROUPS = ADMIX AUS IND TEJ TRJ
THREADS = 16
OPT = -d -v
TEST := $(dir $(realpath $(firstword $(MAKEFILE_LIST))))
PATH := $(realpath $(TEST)/../bin):$(PATH)

help:
	@echo "Targets:"
	@echo "\tinput"
	@echo "\trun"
	@echo "\ttrans"
	@echo "\tgroup"
	@echo "\ttest_vcf"
	@echo "\tld"
	@echo "\tblocks"
	@echo "\tplot"
	@echo "\tall"

all:  input run trans group test_vcf ld blocks plot

input: $(IPREFIX).ped
	$(PIQUE-INPUT) $(OPT) -i $(IPREFIX) -o $(OPREFIX) -p $(PHENO) -k -c -e 5
	touch input

run: input
	$(PIQUE-RUN) $(OPT) -i $(OPREFIX) -k IBS -c $(OPREFIX).covar -n $(THREADS)
	touch run

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
	$(WGET) http://ricediversity.org/data/sets/44kgwas/RiceDiversity.44K.MSU6.Genotypes_PLINK.zip
	unzip RiceDiversity.44K.MSU6.Genotypes_PLINK.zip
	mv -i ./RiceDiversity_44K_Genotypes_PLINK/sativas* .
	rmdir RiceDiversity_44K_Genotypes_PLINK
	rm RiceDiversity.44K.MSU6.Genotypes_PLINK.zip

test_vcf: test_vcf.vcf
	rm -rf $@_GWAS
	$(PIQUE-INPUT) $(OPT) -i $@ -f vcf -o $@_GWAS -p $(PHENO) -k -c -e 5
	$(PIQUE-RUN) -i $@_GWAS -k IBS -c $@_GWAS.covar -n $(THREADS)

test_vcf.vcf: $(IPREFIX).ped $(IPREFIX).map
	$(PLINK) --file $(IPREFIX) --recode vcf --out $(basename $@)

ld: input
	$(PLINK) --file $(IPREFIX)_recode12 --r --out $(OPREFIX)_r
	$(PLINK) --file $(IPREFIX)_recode12 --r2 --out $(OPREFIX)_r2

blocks:
	$(PLINK) --file $(IPREFIX)_recode12 --blocks no-pheno-req --out $(OPREFIX)_blocks

plot: $(PLOT_FILES) run
	GWAS_manhattanplots -b 0.1 -q -y 12 -i $<

clean:
	rm -rf RiceDiversity_44K_Genotypes_PLINK sativas413*
	rm -rf __MACOSX
	rm -rf $(OPREFIX)* $(TPREFIX)* $(GPREFIX)*
	rm -rf test_vcf*
	rm -f input run
	rm -f *.bh *.png
	rm -f .pversion *.log

clobber: clean
	rm -f RiceDiversity.44K.MSU6.Genotypes_PLINK.zip
	rm -f $(IPREFIX)_recode12.*
