# modified vcf4.2 to vcf4.0 file for vcf2diploid
# keep two mouse strains of interest (column)
# keep PASS filtered (row) 
# zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz | awk 'BEGIN{OFS="\t"} (NR <104){print $0}; (NR >103){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$25,$11}' \
#| awk 'BEGIN{OFS="\t"} (NR <105){print $0}; (NR >104 && $7=="PASS"){print $0}' |sed 's/1\/1/1|1/g' | sed 's/0\/0/0|0/g' > CAST_129S1.mgp.v5.snps.dbSNP142.vcf

# filter out genotypes such as 0|., .|.,1|.,... only keep genotypes 1|1, 0|0, 0|1, 1|0, 2|1, 1|2,...
# filter out snps where both pareants are identical to ref genome. ie 0|0, 0|0 (child also 0|0)
# only keep snps where both parents have homolog alleles
# add artificial F1 Genotype

#####CAST and B6
# PersonalGenome_inCurrentDir_make.bsh
cd /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels
zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz |python ../make_vcf_for_vcf2diploid.py P.CAST_M.B6_F1hybrid.snps.forAlleleSeq.vcf 104 25 &
zcat mgp.v5.merged.indels.dbSNP142.normed.vcf.gz | python ../make_vcf_for_vcf2diploid.py P.CAST_M.B6_F1hybrid.indels.forAlleleSeq.vcf 70 25 &
samtools view -s 0.1 -b CAST_EiJ.bam > CAST_EiJ.subsample.bam &
wait
wc -l P.CAST_M.B6_F1hybrid.snps.forAlleleSeq.vcf  #17494363 - 104 = 17494259 
tail -n 17494259 P.CAST_M.B6_F1hybrid.snps.forAlleleSeq.vcf >> P.CAST_M.B6_F1hybrid.indels.forAlleleSeq.vcf
mv P.CAST_M.B6_F1hybrid.indels.forAlleleSeq.vcf P.CAST_M.B6_F1hybrid.indelsNsnps.forAlleleSeq.vcf

cd ..
# indels and snps
make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels \
VCF_sampleID=P.CAST.EiJ_M.C57BL.6J MAT_SAMPLE_NAME=C57BL.6J PAT_SAMPLE_NAME=CAST.EiJ \
OUTPUT_SAMPLE_NAME=P.CAST_M.B6_indelsNsnps_CAST.bam \
FILE_NAME_VCF=P.CAST_M.B6_F1hybrid.indelsNsnps.forAlleleSeq.vcf \
FILE_NAME_BAM=CAST_EiJ.bam 

# only snps
make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels \
VCF_sampleID=P.CAST.EiJ_M.C57BL.6J MAT_SAMPLE_NAME=C57BL.6J PAT_SAMPLE_NAME=CAST.EiJ \
OUTPUT_SAMPLE_NAME=P.CAST_M.B6_snps_CAST.subsample.bam \
FILE_NAME_VCF=P.CAST_M.B6_F1hybrid.snps.forAlleleSeq.vcf \
FILE_NAME_BAM=CAST_EiJ.subsample.bam

###### CAST and 129
cd /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels
zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz |python ../make_vcf_for_vcf2diploid.py P.CAST_M.129S1_F1hybrid.snps.forAlleleSeq.vcf 104 25 11 &
zcat mgp.v5.merged.indels.dbSNP142.normed.vcf.gz | python ../make_vcf_for_vcf2diploid.py P.CAST_M.129S1_F1hybrid.indels.forAlleleSeq.vcf 70 25 11 &
samtools view -s 0.1 -b 129S1_SvImJ.bam > 129S1_SvImJ.subsample.bam &
wait
wc -l P.CAST_M.129S1_F1hybrid.snps.forAlleleSeq.vcf #19514554
tail -n 19514450 P.CAST_M.129S1_F1hybrid.snps.forAlleleSeq.vcf >> P.CAST_M.129S1_F1hybrid.indels.forAlleleSeq.vcf
mv P.CAST_M.129S1_F1hybrid.indels.forAlleleSeq.vcf P.CAST_M.129S1_F1hybrid.indelsNsnps.forAlleleSeq.vcf

cd ..
# CAST and 129
# use CAST_EiJ.subsample.bam
make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels \
VCF_sampleID=P.CAST.EiJ_M.129S1.SvImJ MAT_SAMPLE_NAME=129S1.SvImJ PAT_SAMPLE_NAME=CAST.EiJ \
OUTPUT_SAMPLE_NAME=P.CAST_M.129S1_indelsNsnps_CAST.subsample.bam \
FILE_NAME_VCF=P.CAST_M.129S1_F1hybrid.indelsNsnps.forAlleleSeq.vcf \
FILE_NAME_BAM=CAST_EiJ.subsample.bam > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.129S1_indelsNsnps_CAST.subsample.bam.log 2>&1

# use 129S1_SvImJ.subsample.bam
make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels \
VCF_sampleID=P.CAST.EiJ_M.129S1.SvImJ MAT_SAMPLE_NAME=129S1.SvImJ PAT_SAMPLE_NAME=CAST.EiJ \
OUTPUT_SAMPLE_NAME=P.CAST_M.129S1_indelsNsnps_129S1.subsample.bam \
FILE_NAME_VCF=P.CAST_M.129S1_F1hybrid.indelsNsnps.forAlleleSeq.vcf \
FILE_NAME_BAM=129S1_SvImJ.subsample.bam > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.129S1_indelsNsnps_129S1.subsample.bam.log 2>&1
# here


