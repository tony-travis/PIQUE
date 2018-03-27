#@(#)Makefile  2017-10-05  A.J.Travis and A.Douglas

#
# PIQUE - Parallel Identification Of QTL's using EMMAX
#

WGET = wget --no-check-certificate 
PLINK = p-link
PIQUE-INPUT = pique-input
PIQUE-RUN = pique-run
PHENO_GROUP = rice_phenotype_group.txt
PHENO = rice_phenotype.txt
IPREFIX = sativas413
OPREFIX = sativas_GWAS
TPREFIX = sativas_GWAS_trans
GPREFIX = sativas_group_GWAS
GROUPS = ADMIX AUS IND TEJ TRJ
THREADS = 16
OPT = -d -v
TEST := $(dir $(realpath $(firstword $(MAKEFILE_LIST))))
PATH := $(realpath $(TEST)/../bin):$(PATH)

all:  input run

debug:
	which GWAS_manhattanplots

input: $(IPREFIX).ped
	$(PIQUE-INPUT) $(OPT) -i $(IPREFIX) -o $(OPREFIX) -p $(PHENO) -k -c -e 5

run: $(OPREFIX)
	$(PIQUE-RUN) $(OPT) -i $(OPREFIX) -k IBS -c $(OPREFIX).covar -n $(THREADS)

trans: $(IPREFIX).ped
	$(PIQUE-INPUT) $(OPT) -i $(IPREFIX) -o $(TPREFIX) -p $(PHENO) -k -c -e 5
	$(PIQUE-RUN) $(OPT) -i $(TPREFIX) -k IBS -c $(TPREFIX).covar -n $(THREADS) -t BC

group: $(IPREFIX).ped
	$(PIQUE-INPUT) $(OPT) -i $(IPREFIX) -o $(GPREFIX) -p $(PHENO_GROUP) -k -c -e 5
	$(PIQUE-RUN) $(OPT) -i $(GPREFIX) -k IBS -c $(GPREFIX).covar -n $(THREADS)
	for group in $(GROUPS); do \
	    prefix=$(GPREFIX)\_$$group; \
	    $(PIQUE-RUN) $(OPT) -i $$prefix -k IBS -n $(THREADS); \
	done

$(IPREFIX).ped:
	$(WGET) http://ricediversity.org/data/sets/44kgwas/RiceDiversity.44K.MSU6.Genotypes_PLINK.zip
	unzip RiceDiversity.44K.MSU6.Genotypes_PLINK.zip
	mv -i ./RiceDiversity_44K_Genotypes_PLINK/sativas* .
	rmdir RiceDiversity_44K_Genotypes_PLINK
	rm RiceDiversity.44K.MSU6.Genotypes_PLINK.zip

test: test.vcf
	$(PIQUE-INPUT) $(OPT) -i $@ -f vcf -o $@_vcf -p $(PHENO) -k -c -e 5
	$(PIQUE-RUN) -i $@_vcf -k IBS -c $@_vcf.covar -t $(THREADS)

test.vcf: $(IPREFIX).ped $(IPREFIX).map
	$(PLINK) --file $(IPREFIX) --recode-vcf --out $(basename $@)

ld:
	$(PLINK) --file $(IPREFIX)_recode12 --r --out $(OPREFIX)_r
	$(PLINK) --file $(IPREFIX)_recode12 --r2 --out $(OPREFIX)_r2

blocks:
	$(PLINK) --file $(IPREFIX)_recode12 --blocks no-pheno-req --out $(OPREFIX)_blocks

diff:
	-diff -x '.??*' -x Makefile -x '*.log' -rq ../ref .

plot: files.txt
	GWAS_manhattanplots -b 0.1 -q -y 12 -i files.txt

clean:
	rm -f RiceDiversity.44K.MSU6.Genotypes_PLINK.zip
	rm -rf RiceDiversity_44K_Genotypes_PLINK sativas413*
	rm -rf __MACOSX
	rm -rf $(OPREFIX)* $(GPREFIX)*
	rm -f .pversion *.log

clobber: clean
	rm -f RiceDiversity.44K.MSU6.Genotypes_PLINK.zip
	rm -f $(IPREFIX)_recode12.*
