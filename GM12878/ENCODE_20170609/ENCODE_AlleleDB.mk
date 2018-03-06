# generate a lot of temp file
# make sure have enough disk space before start
# prepare RunAlleleDN.sh


Initial_steps:
#download_fastq:
	# download chipseq fastq.gz and bed file from ENCODE using the files.txt
	xargs -n 1 curl -O -L < files.txt
	#xargs -n 1 curl -O -L < bed_files.txt
#download_bed:
#examine_md5sum:
	# examine md5sum to see if file is correctly downloaded
	md5sum *.fastq.gz > ENCODE_chipseq_md5sum
	python md5sum_examiner.py metadata.tsv ENCODE_chipseq_md5sum
#Run_AlleleDB:
	# run Allele DB pipeline
	# examine RunAlleleDN.sh to modify the path to diploid genome, snp.bed, ...
	sh RunAlleleDN.sh

#Organize_counts.txt:
	#use ENCODE_alleleDB.sh

	

	


