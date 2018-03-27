#@(#)Makefile 2017-01-13  A.J.Travis and A.Douglas

#
# PIQUE documentation
#

TARGETS = pique_manual.pdf pique_tutorial.pdf

all: $(TARGETS)

%.pdf: %.odt
	lowriter --headless --convert-to pdf $<

clean:
	rm -f *.log

clobber: clean
	rm -f $(TARGETS)
