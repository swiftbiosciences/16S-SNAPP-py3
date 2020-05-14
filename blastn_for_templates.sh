#!/bin/bash
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

#Run blastn three times with different DBs and queries to make increasingly
#more complete sets of matches between reads and reference

##Format a blastDB of asv for second blast in order to capture a more complete
#set of matches
echo -e "\nMaking ASV blast DB...\n    Starts: $(date)">>$log
start=$(date +%s.%N)
$MAKEBLASTDB \
    -dbtype nucl \
    -in asv_uniq.fasta  \
    -out asv

echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Making ASV blastDB Runtime: $runtime sec" >> $log

##Blast the pre-formatted reference DB to find the hits above the cutoff
echo -e "\nRunning first blast...\n    Starts: $(date)">>$log
start=$(date +%s.%N)
$BLASTN \
    -num_threads 4 \
    -db ${RDP_FULL_BLAST} \
    -query asv_uniq.fasta \
    -out blastn_1.txt \
    -outfmt "6 std qlen slen" \
    -perc_identity 99 \
    -qcov_hsp_perc 0.99

echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    First Runtime: $runtime sec" >> $log

#obtain top hits of each query sequence from all samples
echo -e "\nFiltering the first blast...\n    Starts: $(date)">>$log
start=$(date +%s.%N)

cat blastn_1.txt \
    | awk -F '\t' '{if($3 >= 99 && ($4/$13) > 0.99) print $2}' \
    | sort -u \
    > reflist.txt
echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Parsing first blast Runtime: $runtime sec" >> $log

#Extract hit sequences from the non-redundant (100%) reference set
template_count=$(< reflist.txt wc -l)
start=$(date +%s.%N)
echo -e "\nExtracting ${template_count} preliminary candidate template sequences\n    Starts: $(date)">>$log
java -jar ${RDPHOME}/ReadSeq.jar \
    select-seqs reflist.txt \
    refset.fasta fasta Y ${RDP_FULL_SEQ}

echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Extraction of preliminary template sequences Runtime: $runtime sec" >> $log

#run blast using preliminary candidate templetes against asv DB in order to obtain a
#more complete set of template-asv matches
echo -e "\nRunning second blast...\n    Starts: $(date)">>$log
start=$(date +%s.%N)
$BLASTN \
    -num_threads 4 \
    -db asv \
    -max_target_seqs 2000 \
    -query refset.fasta \
    -out blastn_2.txt \
    -outfmt "6 std qlen slen" \
    -perc_identity 99

echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Second blast Runtime: $runtime sec" >> $log

cat blastn_2.txt \
    | awk -F '\t' '{if($3 >= 97 && ($4/$14) > 0.97) print $0}' \
    > blastn_2_filtered.txt

echo -e "\nRunning dereplicating blast_2...\n    Starts: $(date)">>$log
start=$(date +%s.%N)
${SCRIPTS}/derep_hitsets.py blastn_2_filtered.txt \
    > blastn_2_derep_IDs.txt
echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Dereplicate Second blast Runtime: $runtime sec" >> $log

echo -e "\nRunning extracting dereplicated blastn_2...\n    Starts: $(date)">>$log
start=$(date +%s.%N)
java -jar ${RDPHOME}/ReadSeq.jar \
    select-seqs blastn_2_derep_IDs.txt \
    blastn_2_derep.fasta fasta Y refset.fasta
echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Extracted dereplicated blast_2 hits Runtime: $runtime sec" >> $log

echo -e "\nRunning making blastDB with dereplicated blastn_2 hits ...\n    Starts: $(date)">>$log
start=$(date +%s.%N)
$MAKEBLASTDB \
    -dbtype nucl \
    -in blastn_2_derep.fasta \
    -out blastn_2_derep
echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Making blastDB with dereplicated blastn_2 hits Runtime: $runtime sec" >> $log

echo -e "\nRunning blastn with all asv seqs against dereplicated blastn_2 hits ...\n    Starts: $(date)">>$log
start=$(date +%s.%N)
$BLASTN \
    -num_threads 4 \
    -db blastn_2_derep \
    -max_target_seqs 2000 \
    -query asv_uniq.fasta \
    -out blastn_3.txt \
    -outfmt "6 std qlen slen" \
    -perc_identity 97 \
    -qcov_hsp_perc 0.97
echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Done blastn with all asv seqs against dereplicated blastn_2 hits Runtime: $runtime sec" >> $log

cat blastn_3.txt \
    | awk -F '\t' '{if($3 >= 97 && ($4/$13) > 0.97) print $0}' \
    > blastn_3_filtered.txt

#yank blastn results and create pickle dictionary file for each of all samples
${SCRIPTS}/pickle_blastn_by_sample.py asv_count.csv blastn_3_filtered.txt asv.uc

#Find the top 1 closest (K1) sequence as the template for each ASV
echo -e "\nRunning RDP Seqmatch with all asv PEs  against RDP11.5  ...\n    Starts: $(date)">>$log
start=$(date +%s.%N)
java -jar ${RDPHOME}/SequenceMatch.jar \
    seqmatch \
    -k 1 \
    -s 0.9 \
    -o asv_PE_K1.txt \
    ${RDP_SEQMATCH}/ \
    asv_PE.fasta
echo "    Ends: $(date)">>$log
end=$(date +%s.%N)
runtime=$(python -c "print(${end} - ${start})")
echo "    Done RDP Seqmatch with all asv PEs against RDP11.5  Runtime: $runtime sec" >> $log

#pickle the seqmatch K1 hits
${SCRIPTS}/pickle_k1.py asv_PE_K1.txt
