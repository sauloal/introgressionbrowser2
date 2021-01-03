USE_PYPY=1

ifdef USE_PYPY
$(info USING PYPY)
PIP=pyp
PYTHON=pypy
else
$(info USING PYTHON)
PIP=pip
PYTHON=python
endif

$(info PIP    $(PIP))
$(info PYTHON $(PYTHON))
$(info )
$(info )


.PHONY: help

help:
	@echo "help"
	@echo "display"
	@echo "tree"
	@echo ""
	@echo "runs"
	@echo "  data/150_VCFs_2.50.tar.gz_ib (150)"
	@echo "  data/360_merged_2.50.vcf.gz_ib (360)"
	@echo "  data/annotated_tomato_150.vcf.gz_ib (anno)"
	@echo ""
	@echo "examples"
	@echo "  data/150_VCFs_2.50.tar.gz"
	@echo "  data/360_merged_2.50.vcf.gz"
	@echo ""
	@echo "requirements"



.PHONY: display
display:
	$(PYTHON) ./reader_ui.py data/



.PHONY: tree
tree:
	tree -h --du --charset ASCII data/



.PHONY: runs
runs: data/150_VCFs_2.50.tar.gz_ib data/360_merged_2.50.vcf.gz_ib data/annotated_tomato_150.vcf.gz_ib

data/150_VCFs_2.50.tar.gz_ib: data/150_VCFs_2.50.tar.gz
	$(PYTHON) ./reader.py $<

data/360_merged_2.50.vcf.gz_ib: data/360_merged_2.50.vcf.gz
	$(PYTHON) ./reader.py $<

data/annotated_tomato_150.vcf.gz_ib: data/annotated_tomato_150.vcf.gz data/annotated_tomato_150.vcf.gz.tsv
	$(PYTHON) ./reader.py --rename-tsv data/annotated_tomato_150.vcf.gz.tsv data/annotated_tomato_150.vcf.gz

.PHONY: 150 360 anno
150: data/150_VCFs_2.50.tar.gz_ib
360: data/360_merged_2.50.vcf.gz_ib
anno: data/annotated_tomato_150.vcf.gz_ib



.PHONY: examples
examples: data/150_VCFs_2.50.tar.gz data/360_merged_2.50.vcf.gz

data/150_VCFs_2.50.tar.gz:
	wget --no-clobber https://s3.eu-central-1.amazonaws.com/saulo.ibrowser/150_VCFs_2.50.tar.gz   -O data/150_VCFs_2.50.tar.gz.tmp   && mv data/150_VCFs_2.50.tar.gz.tmp   data/150_VCFs_2.50.tar.gz

data/360_merged_2.50.vcf.gz:
	wget --no-clobber https://github.com/sauloal/introgressionbrowser2/releases/download/tests/360_merged_2.50.vcf.gz -O data/360_merged_2.50.vcf.gz.tmp && mv data/360_merged_2.50.vcf.gz.tmp data/360_merged_2.50.vcf.gz



.PHONY: requirements
requirements:
	$(PIP) install -r requirements.txt
