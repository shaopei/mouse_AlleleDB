	## align mat and pat
	echo ${FASTQ_PATH}/${FASTQ}.fastq.gz
	zcat ${FASTQ_PATH}/${FASTQ}.fastq.gz | gzip -c | ${PL}/alleledb_filter_input.sh ${PL} - | bowtie --best --strata -p 20 -v 2 -m 1 -f ${PGENOME_PATH}/AltRefMother/AltRefMother - > ${FASTQ}.mat.bowtie 2> ${FASTQ}.mat.log & 
	zcat ${FASTQ_PATH}/${FASTQ}.fastq.gz | gzip -c | ${PL}/alleledb_filter_input.sh ${PL} - | bowtie --best --strata -p 20 -v 2 -m 1 -f ${PGENOME_PATH}/AltRefFather/AltRefFather - > ${FASTQ}.pat.bowtie 2> ${FASTQ}.pat.log &
	

zcat DIPLO_CGATA_R1.fastq.gz | python /workdir/sc2457/alleleDB/alleledb_pipeline/fastq2result.py - - |  python /workdir/sc2457/alleleDB/alleledb_pipeline/filter_query.py - - | python ${PL}/alleledb_ConvertTags.py


	/workdir/sc2457/alleleDB/alleledb_pipeline/alleledb_filter_input.sh /workdir/sc2457/alleleDB/alleledb_pipeline/ LEP_ZYG_ATGCA_R1.fastq.gz > LEP_ZYG_ATGCA_R1_filtered.fasta
	/workdir/sc2457/alleleDB/alleledb_pipeline/alleledb_filter_input.sh /workdir/sc2457/alleleDB/alleledb_pipeline/ LEP_ZYG_ATGCA_R2.fastq.gz > LEP_ZYG_ATGCA_R2_filtered.fasta
	
	/workdir/sc2457/alleleDB/alleledb_pipeline/alleledb_filter_input.sh /workdir/sc2457/alleleDB/alleledb_pipeline/ DIPLO_CGATA_R2.fastq.gz > DIPLO_CGATA_R2_filtered.fasta

bowtie --best --strata -p 20 -v 2 -m 1 -f /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/working/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/AltRefMother/AltRefMother DIPLO_CGATA_R2_filtered.fasta > DIPLO_CGATA_R2.mat.bowtie


	bowtie --best --strata -p 20 -v 2 -m 1 -f /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/working/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/AltRefMother/AltRefMother -1 DIPLO_CGATA_R1_filtered.fasta -2 DIPLO_CGATA_R2_filtered.fasta > DIPLO_CGATA_PE.mat.bowtie
	[sc2457@cbsudanko mouse_AlleleSpecific]$ bowtie --best --strata -p 20 -v 2 -m 1 -f /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/working/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/AltRefMother/AltRefMother -1 DIPLO_CGATA_R1_filtered.fasta -2 DIPLO_CGATA_R2_filtered.fasta > DIPLO_CGATA_PE.mat.bowtie
# reads processed: 3506970
# reads with at least one reported alignment: 3779 (0.11%)
# reads that failed to align: 3502174 (99.86%)
# reads with alignments suppressed due to -m: 1017 (0.03%)
Reported 3779 paired-end alignments to 1 output stream(s)

bowtie --best --strata -p 20 -v 3 -m 1 -f /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/working/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/AltRefMother/AltRefMother -1 DIPLO_CGATA_R1_filtered.fasta -2 DIPLO_CGATA_R2_filtered.fasta > DIPLO_CGATA_PE_v3.mat.bowtie
# reads processed: 3506970
# reads with at least one reported alignment: 6539 (0.19%)
# reads that failed to align: 3498487 (99.76%)
# reads with alignments suppressed due to -m: 1944 (0.06%)
Reported 6539 paired-end alignments to 1 output stream(s)

bowtie -p 20 -v 2 -m 1 -f /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/working/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/AltRefMother/AltRefMother -1 DIPLO_CGATA_R1_filtered.fasta -2 DIPLO_CGATA_R2_filtered.fasta > DIPLO_CGATA_PE_v2Nobest.mat.bowtie
# reads processed: 3506970
# reads with at least one reported alignment: 9854 (0.28%)
# reads that failed to align: 3496077 (99.69%)
# reads with alignments suppressed due to -m: 1039 (0.03%)
Reported 9854 paired-end alignments to 1 output stream(s)

######################
## 1-alignment #######
######################
Wed Jul 26 16:49:36 EDT 2017
/workdir/sc2457/mouse_AlleleSpecific/DIPLO_CGATA_R1.fastq.gz
14070452
/workdir/sc2457/mouse_AlleleSpecific/DIPLO_CGATA_R1.fastq.gz 14070452
2020309 DIPLO_CGATA_R1.mat.bowtie
2018035 DIPLO_CGATA_R1.pat.bowtie
Wed Jul 26 16:52:55 EDT 2017
# reads processed: 3513004
# reads with at least one reported alignment: 2020309 (57.51%)
# reads that failed to align: 966356 (27.51%)
# reads with alignments suppressed due to -m: 526339 (14.98%)
Reported 2020309 alignments to 1 output stream(s)
# reads processed: 3513004
# reads with at least one reported alignment: 2018035 (57.44%)
# reads that failed to align: 967478 (27.54%)
# reads with alignments suppressed due to -m: 527491 (15.02%)
Reported 2018035 alignments to 1 output stream(s)