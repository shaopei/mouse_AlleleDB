make -f makeMousePersonalGenome.mk DATA_DIR=REL-1505-SNPs_Indels VCF_sampleID=P.CAST.EiJ_M.C57BL.6J OUTPUT_SAMPLE_NAME=P.CAST_M.B6_indelsNsnps_CAST.bam \
FILE_NAME_VCF=P.CAST_M.B6_F1hybrid.indelsNsnps.forAlleleSeq.vcf FILE_NAME_BAM=CAST_EiJ.bam

java -Xmx100000000000 -jar /workdir/sc2457/tools/vcf2diploid_v0.2.6a/vcf2diploid.jar -id P.CAST.EiJ_M.C57BL.6J \
-pass -chr /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/GRCm38_68.fa \
-vcf REL-1505-SNPs_Indels/P.CAST_M.B6_F1hybrid.indelsNsnps.forAlleleSeq.vcf \
-outDir REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam >& REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/vcf2diploid.log


bowtie-build --offrate 2 REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/1_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/2_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/3_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/4_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/5_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/6_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/7_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/8_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/9_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/10_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/11_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/12_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/13_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/14_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/15_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/16_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/17_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/18_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/19_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/X_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/Y_P.CAST.EiJ_M.C57BL.6J_maternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/MT_P.CAST.EiJ_M.C57BL.6J_maternal.fa REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam_maternal > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/bowtie_build.maternal.log


bowtie-build --offrate 2 REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/1_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/2_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/3_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/4_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/5_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/6_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/7_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/8_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/9_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/10_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/11_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/12_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/13_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/14_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/15_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/16_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/17_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/18_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/19_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/X_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/Y_P.CAST.EiJ_M.C57BL.6J_paternal.fa,REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/MT_P.CAST.EiJ_M.C57BL.6J_paternal.fa REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam_paternal > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/bowtie_build.paternal.log

# add dad (-d) and mom(-m)
/workdir/sc2457/tools/vcf2diploid_v0.2.6a/vcf2snp -c P.CAST.EiJ_M.C57BL.6J -p 1 -r 1 -s 0 -d CAST.EiJ -m C57BL.6J REL-1505-SNPs_Indels/P.CAST_M.B6_F1hybrid.indelsNsnps.forAlleleSeq.vcf > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp

source ~/.bashrc

cat REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp | awk '{print $1"\t"($2-1000)"\t"($2+1000)"\t"$1"_"($2-1000)"_"($2+1000)}' | grep -v "-" > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.bed
intersectBed -sorted -wo -bed -a REL-1505-SNPs_Indels/CAST_EiJ.bam -b REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.bed \
| awk '{print $16}' | sort | uniq -c > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.bed.counts

ERROR: chromomsome sort ordering for file REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.bed is inconsistent with other files. Record was:
10      3100362 3102362 10_3100362_3102362
sort: write failed: /tmp/sorthAEAxN: No space left on device
cat REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.bed.counts | awk '{ total += $1; count++ } END { print total/count }' > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.bed.meancount
awk: cmd. line:1: fatal: division by zero attempted
make: *** [REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/P.CAST_M.B6_indelsNsnps_CAST.bam.alleleSeqInput.snp] Error 2


make -f makeMousePersonalGenome_2ndpart.mk DATA_DIR=REL-1505-SNPs_Indels VCF_sampleID=P.CAST.EiJ_M.C57BL.6J OUTPUT_SAMPLE_NAME=P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo FILE_NAME_VCF=P.CAST_M.B6_F1hybrid.indelsNsnps.forAlleleSeq.vcf FILE_NAME_BAM=CAST_EiJ.bam
cat REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.snp.bed.counts | awk '{ total += $1; count++ } END { print total/count }' > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.snp.bed.meancount
cat REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.snp.bed.counts | sort -n -k 1 | awk '{ lines[NR]=$0; } END { print lines[int(NR/2)+1] }' | awk '{print $1}' > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.snp.bed.mediancount
cat REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.snp.bed.counts | awk '{getline avg<"REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.snp.bed.mediancount"; print $2"\t"$1/avg }' > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/tmp.cnv
cat REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.snp | awk '{print $1"_"($2-1000)"_"($2+1000)"\t"$0}' | sort -k 1 > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/tmp.snp
## Combine the snp data and the cnv data by region ID
awk 'NR==FNR {h[$1] = $0; next} {print h[$1]"\t"$0}' REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/tmp.snp REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/tmp.cnv > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.snpANDcnv


## Create the snp and cnv input files for alleleseq
cat REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.snpANDcnv | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.alleleSeqInput.snp
echo -e "chrm\tsnppos\trd" > REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.alleleSeqInput.cnv
cat REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.snpANDcnv | awk '{print $2"\t"$3"\t"$10}' >> REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo/P.CAST_M.B6_indelsNsnps_CAST.bam_NoHomo.alleleSeqInput.cnv




