BASE=/workdir/sc2457/
PL:= $(BASE)/alleleDB/alleledb_pipeline_mouse
SNPS:=$(BASE)/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.alleleSeqInput.snp
CNVS:=$(BASE)/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.alleleSeqInput.cnv
BNDS:=hits.bed
MAPS:=$(BASE)/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/%s_P.CAST.EiJ_M.C57BL.6J.map
FDR_SIMS:=5
FDR_CUTOFF:=0.1

PREFIX:=NULL

sourcefiles:=$(PREFIX)     
countfiles_plus:=$(PREFIX).plus.cnt  
countfiles_minus:=$(PREFIX).minus.cnt  

## bowtie files that contain aligned reads
MATBOWTIE:=$(PREFIX).mat.bowtie
PATBOWTIE:=$(PREFIX).pat.bowtie

MATBOWTIE_PLUS:=$(PREFIX).mat.bowtie_plus
PATBOWTIE_PLUS:=$(PREFIX).pat.bowtie_plus
MATBOWTIE_MINUS:=$(PREFIX).mat.bowtie_minus
PATBOWTIE_MINUS:=$(PREFIX).pat.bowtie_minus

## file that contains read IDs to be removed from bowtie files above
READS2FILTER:=originalmatpatreads.toremove.ids

target:interestingHets_plus.txt

#Generate Strand Specific bowtie output

$(PATBOWTIE_PLUS):$(PATBOWTIE_MINUS)
$(PATBOWTIE_MINUS):$(PATBOWTIE)
	python $(PL)/filter_reads_out.py $(PATBOWTIE) - $(READS2FILTER) | python ${PL}/seperate_strand_of_bowtie_output_alleleDB.py $(PATBOWTIE)

$(MATBOWTIE_PLUS):$(MATBOWTIE_MINUS)
$(MATBOWTIE_MINUS):$(MATBOWTIE)
	python $(PL)/filter_reads_out.py $(MATBOWTIE) - $(READS2FILTER) | python ${PL}/seperate_strand_of_bowtie_output_alleleDB.py $(MATBOWTIE) 


$(countfiles_plus): $(PATBOWTIE_PLUS) $(MATBOWTIE_PLUS)
	bash -c "python $(PL)/MergeBowtie.py \
           $(PATBOWTIE_PLUS) $(MATBOWTIE_PLUS) \
           $(MAPS) | python $(PL)/SnpCounts.py $(SNPS) - $(MAPS) $@"

$(countfiles_minus): $(PATBOWTIE_MINUS) $(MATBOWTIE_MINUS)
	bash -c "python $(PL)/MergeBowtie.py \
           $(PATBOWTIE_MINUS) $(MATBOWTIE_MINUS) \
           $(MAPS) | python $(PL)/SnpCounts.py $(SNPS) - $(MAPS) $@"

check:
	@echo $(sourcefiles)


counts_plus.txt: $(countfiles_plus) $(countfiles_minus)
	python $(PL)/CombineSnpCounts.py 5 $(SNPS) $(BNDS) $(CNVS) counts_plus.txt counts_plus.log $(countfiles_plus)
	python $(PL)/CombineSnpCounts.py 5 $(SNPS) $(BNDS) $(CNVS) counts_minus.txt counts_minus.log $(countfiles_minus)

# calculate false discovery rates
FDR_plus.txt: counts_plus.txt
	python $(PL)/FalsePos.py counts_plus.txt $(FDR_SIMS) $(FDR_CUTOFF) > FDR_plus.txt
	python $(PL)/FalsePos.py counts_minus.txt $(FDR_SIMS) $(FDR_CUTOFF) > FDR_minus.txt

interestingHets_plus.txt: counts_plus.txt FDR_plus.txt
	awk -f $(PL)/finalFilter.awk thresh=$(shell awk 'END {print $$6}' FDR_plus.txt) < counts_plus.txt > interestingHets_plus.txt
	awk -f $(PL)/finalFilter.awk thresh=$(shell awk 'END {print $$6}' FDR_minus.txt) < counts_minus.txt > interestingHets_minus.txt

clean:
	@rm -f FDR*.txt interestingHets*.txt counts*.txt

cleanall: clean
	@rm -f *.cnt

.DELETE_ON_ERROR:

