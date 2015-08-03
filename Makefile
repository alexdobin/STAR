# You want latexmk to *always* run, because make does not have all the info.
.PHONY: 

INPUT_DIR = ./extras/doc-latex
OUTPUT_DIR = ./doc

all : STAR manual

STAR :
	cd source && $(MAKE) $@
manual : $(OUTPUT_DIR)/STARmanual.pdf


%.tex : %.raw
	./raw2tex $< > $@
%.tex : %.dat
	./dat2tex $< > $@

$(OUTPUT_DIR)/STARmanual.pdf : $(INPUT_DIR)/STARmanual.tex $(INPUT_DIR)/parametersDefault.tex
	cd $(INPUT_DIR) && latexmk -pdf -output-directory=../../$(OUTPUT_DIR) -jobname=STARmanual -pdflatex='pdflatex -halt-on-error %O %S -synctex=1 -interaction=nonstopmode --src-specials' -quiet -f -use-make $(basename $(notdir $<))

clean : $(INPUT_DIR)/STARmanual.tex
	cd $(INPUT_DIR) && latexmk -C -output-directory=../../$(OUTPUT_DIR) -jobname=STARmanual $(basename $(notdir $<))
	cd source && $(MAKE) $@
