#!/usr/bin/bash

mkdir -p fastq

list=("normal_rep1" "normal_rep2" "normal_rep3")

for i in ${list[@]}; do 

	perl CIRI_simulator.pl -O ${i} -G arm2.gtf -C 10 -LC 10 -R 1 -LR 1 -L 100 -E 1 -D ./ -CHR1 1 -M 320 -M2 550 -PM 15 -S 70 -S2 70 -SE 0 -PSI 0
done

list1=("tumor_rep1" "tumor_rep2" "tumor_rep3")

for i in ${list1[@]}; do

	perl CIRI_simulator.pl -O ${i} -G arm1.gtf -C 10 -LC 10 -R 1 -LR 1 -L 100 -E 1 -D ./ -CHR1 1 -M 320 -M2 550 -PM 15 -S 70 -S2 70 -SE 0 -PSI 0

done

mv *.fq fastq/

gzip fastq/*.*

rm fastq/*rep2*

cp fastq/normal_rep1_1.fq.gz fastq/normal_rep2_1.fq.gz
cp fastq/normal_rep1_2.fq.gz fastq/normal_rep2_2.fq.gz

cp fastq/tumor_rep1_1.fq.gz fastq/tumor_rep2_1.fq.gz
cp fastq/tumor_rep1_2.fq.gz fastq/tumor_rep2_2.fq.gz

