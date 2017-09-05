

# modified vcf4.2 to vcf4.0 file for vcf2diploid
# keep two mouse strains of interest (column)
# keep PASS filtered (row) 
#zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz | awk 'BEGIN{OFS="\t"} (NR <104){print $0}; (NR >103){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$25,$11}' \
#| awk 'BEGIN{OFS="\t"} (NR <105){print $0}; (NR >104 && $7=="PASS"){print $0}' |sed 's/1\/1/1|1/g' | sed 's/0\/0/0|0/g' > CAST_129S1.mgp.v5.snps.dbSNP142.vcf

# filter out genotypes such as 0|., .|.,1|.,... only keep genotypes 1|1, 0|0, 0|1, 1|0
# filter out snps where both pareants are identical to ref genome. ie 0|0, 0|0 (child also 0|0)
# only keep snps where both parents have homolog alleles
# add artificial F1 Genotype

#####CAST and B6
zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz |python ../make_vcf_for_vcf2diploid.py P.CAST_M.B6_F1hybrid.snps.forAlleleSeq.vcf 104 25 &
zcat mgp.v5.merged.indels.dbSNP142.normed.vcf.gz | python ../make_vcf_for_vcf2diploid.py P.CAST_M.B6_F1hybrid.indels.forAlleleSeq.vcf 70 25 &
wc -l P.CAST_M.B6_F1hybrid.snps.forAlleleSeq.vcf  #17494363
tail -n 17494259 P.CAST_M.B6_F1hybrid.snps.forAlleleSeq.vcf >> P.CAST_M.B6_F1hybrid.indels.forAlleleSeq.vcf
mv P.CAST_M.B6_F1hybrid.indels.forAlleleSeq.vcf P.CAST_M.B6_F1hybrid.indelsNsnps.forAlleleSeq.vcf

# indels and snps
make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels \
VCF_sampleID=P.CAST.EiJ_M.C57BL.6J MAT_SAMPLE_NAME=C57BL.6J PAT_SAMPLE_NAME=CAST.EiJ \
OUTPUT_SAMPLE_NAME=P.CAST_M.B6_indelsNsnps_CAST.bam \
FILE_NAME_VCF=P.CAST_M.B6_F1hybrid.indelsNsnps.forAlleleSeq.vcf \
FILE_NAME_BAM=CAST_EiJ.bam 
# only snps
make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels \
VCF_sampleID=P.CAST.EiJ_M.C57BL.6J MAT_SAMPLE_NAME=C57BL.6J PAT_SAMPLE_NAME=CAST.EiJ \
OUTPUT_SAMPLE_NAME=P.CAST_M.B6_snps_CAST.bam \
FILE_NAME_VCF=P.CAST_M.B6_F1hybrid.snps.forAlleleSeq.vcf \
FILE_NAME_BAM=CAST_EiJ.bam 


###### CAST and 129
zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz |python ../make_vcf_for_vcf2diploid.py P.CAST_M.129S1_F1hybrid.snps.forAlleleSeq.vcf 104 25 11
zcat mgp.v5.merged.indels.dbSNP142.normed.vcf.gz | python ../make_vcf_for_vcf2diploid.py P.CAST_M.129S1_F1hybrid.indels.forAlleleSeq.vcf 70 25 11
# here
wc -l P.CAST_M.129S1_F1hybrid.snps.forAlleleSeq.vcf
tail -n 19514450 P.CAST_M.129S1_F1hybrid.snps.forAlleleSeq.vcf >> P.CAST_M.129S1_F1hybrid.indels.forAlleleSeq.vcf
mv P.CAST_M.129S1_F1hybrid.indels.forAlleleSeq.vcf P.CAST_M.129S1_F1hybrid.indelsNsnps.forAlleleSeq.vcf

# CAST and 129
# use CAST_EiJ.bam
make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels \
VCF_sampleID=P.CAST.EiJ_M.129S1.SvImJ MAT_SAMPLE_NAME=129S1.SvImJ PAT_SAMPLE_NAME=CAST.EiJ \
OUTPUT_SAMPLE_NAME=P.CAST_M.129S1_indelsNsnps_CAST.bam \
FILE_NAME_VCF=P.CAST_M.129S1_F1hybrid.indelsNsnps.forAlleleSeq.vcf \
FILE_NAME_BAM=CAST_EiJ.bam

# use 129S1_SvImJ.bam
make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels \
VCF_sampleID=P.CAST.EiJ_M.129S1.SvImJ MAT_SAMPLE_NAME=129S1.SvImJ PAT_SAMPLE_NAME=CAST.EiJ \
OUTPUT_SAMPLE_NAME=P.CAST_M.129S1_indelsNsnps_129S1.bam \
FILE_NAME_VCF=P.CAST_M.129S1_F1hybrid.indelsNsnps.forAlleleSeq.vcf \
FILE_NAME_BAM=129S1_SvImJ.bam





######### underneath is  stock ######### 




# CAST and 129
# use 129S1_SvImJ.bam
make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels \
VCF_sampleID=P.CAST.EiJ_M.129S1.SvImJ \
OUTPUT_SAMPLE_NAME=P.CAST_M.129S1_indelsNsnps_129S1.bam \
FILE_NAME_VCF=P.CAST_M.129S1_F1hybrid.indelsNsnps.forAlleleSeq.vcf \
FILE_NAME_BAM=129S1_SvImJ.bam



#####CAST and B6

make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels \
VCF_sampleID=P.CAST.EiJ_M.C57BL.6J \
OUTPUT_SAMPLE_NAME=P.CAST_M.B6_indelsNsnps_CAST.bam \
FILE_NAME_VCF=P.CAST_M.B6_F1hybrid.indelsNsnps.forAlleleSeq.vcf \
FILE_NAME_BAM=CAST_EiJ.bam






java -Xmx100000000000 -jar /workdir/sc2457/tools/vcf2diploid_v0.2.6a/vcf2diploid.jar -id P.CAST.EiJ_M.C57BL.6J -pass -chr /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/GRCm38_68.fa -vcf REL-1505-SNPs_Indels/P.CAST_M.B6_F1hybrid.snpsNindels.forAlleleSeq.vcf -outDir REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels >& REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/vcf2diploid.log

bowtie-build --offrate 2 REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/1_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/2_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/3_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/4_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/5_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/6_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/7_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/8_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/9_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/10_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/11_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/12_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/13_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/14_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/15_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/16_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/17_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/18_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/19_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/X_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/Y_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/MT_P.CAST.EiJ_M.C57BL.6J_maternal.fa REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/P.CAST_M.B6_snpsNindels_maternal > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/bowtie_build.maternal.log
bowtie-build --offrate 2 REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/1_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/2_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/3_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/4_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/5_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/6_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/7_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/8_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/9_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/10_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/11_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/12_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/13_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/14_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/15_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/16_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/17_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/18_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/19_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/X_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/Y_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/MT_P.CAST.EiJ_M.C57BL.6J_paternal.fa REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/P.CAST_M.B6_snpsNindels_paternal > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/bowtie_build.paternal.log
/workdir/sc2457/tools/vcf2diploid_v0.2.6a/vcf2snp -c P.CAST.EiJ_M.C57BL.6J -p 1 -r 1 -s 0  REL-1505-SNPs_Indels/P.CAST_M.B6_F1hybrid.snpsNindels.forAlleleSeq.vcf \
> REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/P.CAST_M.B6_snpsNindels.snp

cat REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/P.CAST_M.B6_snpsNindels.snp | awk '{print $1"\t"($2-1000)"\t"($2+1000)"\t"$1"_"($2-1000)"_"($2+1000)}' \
| grep -v "-" > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/P.CAST_M.B6_snpsNindels.snp.bed
intersectBed -sorted -wo -bed -a /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1502-BAM/CAST_EiJ.bam -b REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/P.CAST_M.B6_snpsNindels.snp.bed  | awk '{print $16}' | sort | uniq -c > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/P.CAST_M.B6_snpsNindels.snp.bed.counts



cat REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/P.CAST_M.B6_snpsNindels.snp.bed.counts | awk '{ total += $1; count++ } END { print total/count }' > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/P.CAST_M.B6_snpsNindels.snp.bed.meancount
awk: cmd. line:1: fatal: division by zero attempted
make: *** [REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_snpsNindels/P.CAST_M.B6_snpsNindels.alleleSeqInput.snp] Error 2


#make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels VCF_sampleID=P.CAST.EiJ_M.129S1.SvImJ OUTPUT_SAMPLE_NAME=P.CAST_M.129S1_indels FILE_NAME_VCF=P.CAST_M.129S1_F1hybrid.indels.forAlleleSeq.vcf.gz




#####


# inside
# make -f makePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels OUTPUT_SAMPLE_NAME=129S1_SvImJ_CAST_EiJ FILE_NAME_VCF=129S1_CAST_F1hybrid.forAlleleSeq.vcf

mkdir CAST_EiJ.mgp.v5/PersonalGenome_CAST
mkdir CAST_EiJ.mgp.v5/PersonalGenome_CAST/AltRefMother
mkdir CAST_EiJ.mgp.v5/PersonalGenome_CAST/AltRefFather


java -Xmx100000000000 -jar /workdir/sc2457/tools/vcf2diploid_v0.2.6a/vcf2diploid.jar -id CAST -pass -chr /workdir/sc2457/mouse_AlleleSpecific/mouse_genome/GRCm38_68.fa -vcf CAST_EiJ.mgp.v5/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz -outDir CAST_EiJ.mgp.v5/PersonalGenome_CAST >& CAST_EiJ.mgp.v5/PersonalGenome_CAST/vcf2diploid.log




make -f makePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels OUTPUT_SAMPLE_NAME=129S1_SvImJ_CAST_EiJ FILE_NAME_VCF=129S1_CAST_F1hybrid.forAlleleSeq.vcf
make -f makePersonalGenome.mk DATA_DIR=/path/to/your/VCFdir OUTPUT_SAMPLE_NAME=NA12878 FILE_NAME_BAM=filename.in.DATA_DIR.bam FILE_NAME_VCF=filename.in.DATA_DIR.vcf
mkdir REL-1505-SNPs_Indels/PersonalGenome_129S1_SvImJ_CAST_EiJ
mkdir REL-1505-SNPs_Indels/PersonalGenome_129S1_SvImJ_CAST_EiJ/AltRefMother
mkdir REL-1505-SNPs_Indels/PersonalGenome_129S1_SvImJ_CAST_EiJ/AltRefFather
java -Xmx100000000000 -jar /workdir/sc2457/tools/vcf2diploid_v0.2.6a/vcf2diploid.jar -id 129S1_SvImJ_CAST_EiJ -pass -chr /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/GRCm38_68.fa -vcf REL-1505-SNPs_Indels/129S1_CAST_F1hybrid.forAlleleSeq.vcf -outDir REL-1505-SNPs_Indels/PersonalGenome_129S1_SvImJ_CAST_EiJ >& REL-1505-SNPs_Indels/PersonalGenome_129S1_SvImJ_CAST_EiJ/vcf2diploid.log


#make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels OUTPUT_SAMPLE_NAME=B6_CAST FILE_NAME_VCF=B6_CAST_F1hybrid.forAlleleSeq.vcf
mkdir REL-1505-SNPs_Indels/PersonalGenome_B6_CAST
mkdir REL-1505-SNPs_Indels/PersonalGenome_B6_CAST/AltRefMother
mkdir REL-1505-SNPs_Indels/PersonalGenome_B6_CAST/AltRefFather
java -Xmx100000000000 -jar /workdir/sc2457/tools/vcf2diploid_v0.2.6a/vcf2diploid.jar -id B6_CAST -pass -chr /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/GRCm38_68.fa -vcf REL-1505-SNPs_Indels/B6_CAST_F1hybrid.forAlleleSeq.vcf -outDir REL-1505-SNPs_Indels/PersonalGenome_B6_CAST >& REL-1505-SNPs_Indels/PersonalGenome_B6_CAST/vcf2diploid.log


#make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels OUTPUT_SAMPLE_NAME=P.CAST_M.B6 FILE_NAME_VCF=P.CAST_M.B6_F1hybrid.forAlleleSeq.vcf

