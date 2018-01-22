BASE=/workdir/sc2457/
PL:= /workdir/sc2457/alleleDB/alleledb_pipeline
SNPS:=$(BASE)/SNP/1000genome_vol1.ftp.release.20130502/snp.call.all
CNVS:=$(BASE)/alleleseq.gersteinlab.org/NA12878_diploid_dec16.2012.alleleseq.input/rd_4327183snps_na12878_hg19.txt
BNDS:=hits.bed
MAPS:=$(BASE)/alleleseq.gersteinlab.org/NA12878_diploid/NA12878_diploid_2015_feb5_3versions/1kgp3-svs-pass_NA12878_hg19_150109_w_transcriptome/%s_NA12878.map
FDR_SIMS:=5
FDR_CUTOFF:=0.1
Min_count:=1

PREFIX:=NULL

sourcefiles:=$(PREFIX)     
countfiles:=$(PREFIX).cnt  

## bowtie files that contain aligned reads
MATBOWTIE:=$(PREFIX).mat.bowtie
PATBOWTIE:=$(PREFIX).pat.bowtie

## file that contains read IDs to be removed from bowtie files above
READS2FILTER:=originalmatpatreads.toremove.ids



all: $(PREFIX)_interestingHets.txt


#$(countfiles): $(countfiles).gz
#	zcat $(countfiles).gz > $(countfiles)


#check:
#	@echo $(sourcefiles)


$(PREFIX)_counts.txt: $(countfiles)
	python $(PL)/CombineSnpCounts.py $(Min_count) $(SNPS) $(BNDS) $(CNVS) $(PREFIX)_counts.txt $(PREFIX)_counts.log $(countfiles)

# calculate false discovery rates
$(PREFIX)_FDR.txt: $(PREFIX)_counts.txt
	python $(PL)/FalsePos.py $(PREFIX)_counts.txt $(FDR_SIMS) $(FDR_CUTOFF) > $(PREFIX)_FDR.txt

$(PREFIX)_interestingHets.txt: $(PREFIX)_counts.txt $(PREFIX)_FDR.txt
	awk -f $(PL)/finalFilter.awk thresh=$(shell awk 'END {print $$6}' $(PREFIX)_FDR.txt) < $(PREFIX)_counts.txt > $(PREFIX)_interestingHets.txt

clean:
	@rm -f $(PREFIX)_FDR.txt $(PREFIX)_interestingHets.txt $(PREFIX)_counts.txt

cleanall: clean
	@rm -f *.cnt

.DELETE_ON_ERROR:
