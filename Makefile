.PHONY: examples

examples: data/150_VCFs_2.50.tar.gz data/360_merged_2.50.vcf.gz

data/150_VCFs_2.50.tar.gz:
	wget --no-clobber https://s3.eu-central-1.amazonaws.com/saulo.ibrowser/150_VCFs_2.50.tar.gz   -O data/150_VCFs_2.50.tar.gz.tmp   && mv data/150_VCFs_2.50.tar.gz.tmp   data/150_VCFs_2.50.tar.gz

data/360_merged_2.50.vcf.gz:
	wget --no-clobber https://github.com/sauloal/introgressionbrowser2/releases/download/tests/360_merged_2.50.vcf.gz -O data/360_merged_2.50.vcf.gz.tmp && mv data/360_merged_2.50.vcf.gz.tmp data/360_merged_2.50.vcf.gz

