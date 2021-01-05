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
	@echo "  360_run"
	@echo "  360_debug_run"
	@echo "  150_anno_run"
	@echo "  150_anno_debug_run"
	@echo "  360_anno_run"
	@echo "  360_anno_debug_run"
	@echo ""
	@echo "clean"
	@echo "  360_clean"
	@echo "  360_debug_clean"
	@echo "  150_anno_clean"
	@echo "  150_anno_debug_clean"
	@echo "  360_anno_clean"
	@echo "  360_anno_debug_clean"
	@echo ""
	@echo "downloads"
	@echo "  150_download"
	@echo "  360_download"
	@echo "  360_debug_download (360 symlink)"
	@echo "  150_anno_download"
	@echo "  150_anno_debug_download (150 anno symlink)"
	@echo "  360_anno_download"
	@echo "  360_anno_debug_download (360 anno symlink)"
	@echo ""
	@echo "requirements"



.PHONY: display
display:
	$(PYTHON) ./reader_ui.py data/



.PHONY: tree
tree:
	tree -h --du --charset ASCII data/



.PHONY: runs
runs: 360_run 360_debug_run 150_anno_run 150_anno_debug_run 360_anno_run 360_anno_debug_run

.PHONY: 360_run 360_debug_run 150_anno_run 150_anno_debug_run 360_anno_run 360_anno_debug_run
360_run: data/360_merged_2.50.vcf.gz_ib
360_debug_run: data/370_merged_2.50.vcf.gz_ib
150_anno_run: data/annotated_tomato_150.vcf.gz_ib
150_anno_debug_run: data/annotated_tomato_160.vcf.gz_ib
360_anno_run: data/annotated_tomato_360.vcf.gz_ib
360_anno_debug_run: data/annotated_tomato_370.vcf.gz_ib


data/360_merged_2.50.vcf.gz_ib: data/360_merged_2.50.vcf.gz
	$(PYTHON) ./reader.py $^ || (rm -rfv $@ && exit 1)

data/370_merged_2.50.vcf.gz_ib: data/370_merged_2.50.vcf.gz
	$(PYTHON) ./reader.py --debug $^ || (rm -rfv $@ && exit 1)

data/annotated_tomato_150.vcf.gz_ib: data/annotated_tomato_150.vcf.gz.tsv data/annotated_tomato_150.vcf.gz
	$(PYTHON) ./reader.py --rename-tsv $^ || (rm -rfv $@ && exit 1)

data/annotated_tomato_160.vcf.gz_ib: data/annotated_tomato_160.vcf.gz.tsv data/annotated_tomato_160.vcf.gz
	$(PYTHON) ./reader.py --debug --rename-tsv $^ || (rm -rfv $@ && exit 1)

data/annotated_tomato_360.vcf.gz_ib: data/annotated_tomato_360.vcf.gz
	$(PYTHON) ./reader.py $^ || (rm -rfv $@ && exit 1)

data/annotated_tomato_370.vcf.gz_ib: data/annotated_tomato_370.vcf.gz
	$(PYTHON) ./reader.py --debug $^ || (rm -rfv $@ && exit 1)






.PHONY: clean
clean: 360_clean 360_debug_clean 150_anno_clean 150_anno_debug_clean 360_anno_clean 360_anno_debug_clean

.PHONY: 360_clean 360_debug_clean 150_anno_clean 150_anno_debug_clean 360_anno_clean 360_anno_debug_clean
360_clean:
	rm -rfv data/360_merged_2.50.vcf.gz_ib      || true

360_debug_clean:
	rm -rfv data/370_merged_2.50.vcf.gz_ib      || true

150_anno_clean:
	rm -rfv data/annotated_tomato_150.vcf.gz_ib || true

150_anno_debug_clean:
	rm -rfv data/annotated_tomato_160.vcf.gz_ib || true

360_anno_clean:
	rm -rfv data/annotated_tomato_360.vcf.gz_ib || true

360_anno_debug_clean:
	rm -rfv data/annotated_tomato_370.vcf.gz_ib || true



.PHONY: downloads
downloads: 150_download 360_download 360_debug_download 150_anno_download 150_anno_debug_download 360_anno_download 360_anno_debug_download

.PHONY: 150_download 360_download 360_debug_download 150_anno_download 150_anno_debug_download 360_anno_download 360_anno_debug_download

150_download: data/150_VCFs_2.50.tar.gz
360_download: data/360_merged_2.50.vcf.gz
360_debug_download: data/370_merged_2.50.vcf.gz
150_anno_download: data/annotated_tomato_150.vcf.gz
150_anno_debug_download: data/annotated_tomato_160.vcf.gz
360_anno_download: data/annotated_tomato_360.vcf.gz
360_anno_debug_download: data/annotated_tomato_370.vcf.gz

data/150_VCFs_2.50.tar.gz:
	wget --no-clobber https://s3.eu-central-1.amazonaws.com/saulo.ibrowser/$(notdir $@)   -O $@.tmp || (rm -v $@.tmp && exit 1)
	mv $@.tmp $@

data/360_merged_2.50.vcf.gz:
	wget --no-clobber https://github.com/sauloal/introgressionbrowser2/releases/download/tests/$(notdir $@) -O $@.tmp || (rm -v $@.tmp && exit 1)
	mv $@.tmp $@

data/370_merged_2.50.vcf.gz: data/360_merged_2.50.vcf.gz
	cd $(dir $<) && ln -s $(notdir $<) $(notdir $@)

data/annotated_tomato_150.vcf.gz:
	wget --no-clobber https://s3.eu-central-1.amazonaws.com/saulo.ibrowser/$(notdir $@)   -O $@.tmp || (rm -v $@.tmp && exit 1)
	mv $@.tmp $@

data/annotated_tomato_160.vcf.gz: data/annotated_tomato_150.vcf.gz
	cd $(dir $<) && ln -s $(notdir $<) $(notdir $@)

data/annotated_tomato_360.vcf.gz:
	wget --no-clobber https://github.com/sauloal/introgressionbrowser2/releases/download/tests/$(notdir $@)   -O $@.tmp || (rm -v $@.tmp && exit 1)
	mv $@.tmp $@

data/annotated_tomato_370.vcf.gz: data/annotated_tomato_360.vcf.gz
	cd $(dir $<) && ln -s $(notdir $<) $(notdir $@)

data/annotated_tomato_160.vcf.gz.tsv: data/annotated_tomato_150.vcf.gz.tsv
	cd $(dir $<) && ln -s $(notdir $<) $(notdir $@)



.PHONY: requirements
requirements:
	$(PIP) install -r requirements.txt
