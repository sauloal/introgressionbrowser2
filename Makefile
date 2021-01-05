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
	@echo "  run_360"
	@echo "  run_370 (360 debug)"
	@echo "  run_150_anno"
	@echo "  run_160_anno (150 anno debug)"
	@echo "  run_360_anno"
	@echo "  run_370_anno (360 anno debug)"
	@echo ""
	@echo "clean"
	@echo "  360_clean"
	@echo "  370_clean"
	@echo "  150_anno_clean"
	@echo "  160_anno_clean"
	@echo "  360_anno_clean"
	@echo "  370_anno_clean"
	@echo ""
	@echo "downloads"
	@echo "  150_download"
	@echo "  360_download"
	@echo "  370_download (360 symlink)"
	@echo "  150_anno_download"
	@echo "  160_anno_download (150 anno symlink)"
	@echo "  360_anno_download"
	@echo "  370_anno_download (360 anno symlink)"
	@echo ""
	@echo "requirements"



.PHONY: display
display:
	$(PYTHON) ./reader_ui.py data/



.PHONY: tree
tree:
	tree -h --du --charset ASCII data/



.PHONY: runs
runs: run_360 run_370 run_150_anno run_160_anno run_360_anno run_370_anno

.PHONY: run_360 run_370 run_150_anno run_160_anno run_360_anno run_370_anno
run_360: data/360_merged_2.50.vcf.gz_ib
run_370: data/370_merged_2.50.vcf.gz_ib
run_150_anno: data/annotated_tomato_150.vcf.gz_ib
run_160_anno: data/annotated_tomato_160.vcf.gz_ib
run_360_anno: data/annotated_tomato_360.vcf.gz_ib
run_370_anno: data/annotated_tomato_370.vcf.gz_ib


data/360_merged_2.50.vcf.gz_ib: data/360_merged_2.50.vcf.gz
	$(PYTHON) ./reader.py $^

data/370_merged_2.50.vcf.gz_ib: data/370_merged_2.50.vcf.gz
	$(PYTHON) ./reader.py --debug $^

data/annotated_tomato_150.vcf.gz_ib: data/annotated_tomato_150.vcf.gz.tsv data/annotated_tomato_150.vcf.gz
	$(PYTHON) ./reader.py --rename-tsv $^

data/annotated_tomato_160.vcf.gz_ib: data/annotated_tomato_160.vcf.gz.tsv data/annotated_tomato_160.vcf.gz
	$(PYTHON) ./reader.py --debug --rename-tsv $^

data/annotated_tomato_360.vcf.gz_ib: data/annotated_tomato_360.vcf.gz
	$(PYTHON) ./reader.py $^

data/annotated_tomato_370.vcf.gz_ib: data/annotated_tomato_370.vcf.gz
	$(PYTHON) ./reader.py --debug $^






.PHONY: clean
clean: 360_clean 370_clean 150_anno_clean 160_anno_clean 360_anno_clean 370_anno_clean

.PHONY: 360_clean 370_clean 150_anno_clean 160_anno_clean 360_anno_clean 370_anno_clean
360_clean:
	rm -rfv data/360_merged_2.50.vcf.gz_ib      || true

370_clean:
	rm -rfv data/370_merged_2.50.vcf.gz_ib      || true

150_anno_clean:
	rm -rfv data/annotated_tomato_150.vcf.gz_ib || true

160_anno_clean:
	rm -rfv data/annotated_tomato_160.vcf.gz_ib || true

360_anno_clean:
	rm -rfv data/annotated_tomato_360.vcf.gz_ib || true

370_anno_clean:
	rm -rfv data/annotated_tomato_370.vcf.gz_ib || true



.PHONY: downloads
downloads: 150_download 360_download 370_download 150_anno_download 160_anno_download 360_anno_download 370_anno_download

.PHONY: 150_download 360_download 370_download 150_anno_download 160_anno_download 360_anno_download 370_anno_download

150_download: data/150_VCFs_2.50.tar.gz
360_download: data/360_merged_2.50.vcf.gz
370_download: data/370_merged_2.50.vcf.gz
150_anno_download: data/annotated_tomato_150.vcf.gz
160_anno_download: data/annotated_tomato_160.vcf.gz
360_anno_download: data/annotated_tomato_360.vcf.gz
370_anno_download: data/annotated_tomato_370.vcf.gz

data/150_VCFs_2.50.tar.gz:
	wget --no-clobber https://s3.eu-central-1.amazonaws.com/saulo.ibrowser/$(notdir $@)   -O $@.tmp || (rm -v $@.tmp && exit 1)
	mv $@.tmp $@

data/360_merged_2.50.vcf.gz:
	wget --no-clobber https://github.com/sauloal/introgressionbrowser2/releases/download/tests/$(notdir $@) -O $@.tmp || (rm -v $@.tmp && exit 1)
	mv $@.tmp $@

data/370_merged_2.50.vcf.gz: data/360_merged_2.50.vcf.gz
	ln -s $< $@

data/annotated_tomato_150.vcf.gz:
	wget --no-clobber https://s3.eu-central-1.amazonaws.com/saulo.ibrowser/$(notdir $@)   -O $@.tmp || (rm -v $@.tmp && exit 1)
	mv $@.tmp $@

data/annotated_tomato_160.vcf.gz: data/annotated_tomato_150.vcf.gz
	ln -s $< $@

data/annotated_tomato_360.vcf.gz:
	wget --no-clobber https://s3.eu-central-1.amazonaws.com/saulo.ibrowser/$(notdir $@)   -O $@.tmp || (rm -v $@.tmp && exit 1)
	mv $@.tmp $@

data/annotated_tomato_370.vcf.gz: data/annotated_tomato_360.vcf.gz
	ln -s $< $@

data/annotated_tomato_160.vcf.gz.tsv: data/annotated_tomato_150.vcf.gz.tsv
	ln -s $< $@



.PHONY: requirements
requirements:
	$(PIP) install -r requirements.txt
