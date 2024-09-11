Data were generated using the following Makefile:

the `data/genomics/sarscov2/wgsim` directory contains 4 bams WT, DEL, DUP, INV.
Fastqs for those reads were generated using wgsim.

The coverage of those bams looks like this:

```
MT192765.1 (29.8Kbp)
>    24.1 │              ▂                           ▂                                   ▇     ▅                                     ▁  ▅  █                 █         ▁                             │ Number of reads: 8000
>    21.4 │    ▁ █       █ ▇     █  ▇        ▅  ▆    █ ▁▅ ▁▅                          ▁▆ █ ▃   █                    █              ▁ █ ▇█▃▆█         █    ▂  █ █▂  ▄▄  █                             │ 
>    18.7 │    █▂█▃▂    ▂█ █▅ ▂  █▃ █▅    ▂  █▆ █▁▁▅ █ ██ ██      █ ▅▁    ▁  ▃    ▃ █▄██▃█ █▅ ▇█▆▂ ▃       ▂        █  ▆ ▁▇  ▆  ▁  █▁█▇█████  ▂   █ ▂█ ▇▃▂█▇▅█ ██ ▅██▂ █▅ ▇▆█ ▁▁ ▁ ▂█ █     ▆▅▂▄▆█   │ Covered bases:   29.8Kbp
>    16.0 │   ██████▃▃▇▆██▃██▅█▁▆██▅██▃▄  █ ▄██ ████ █ █████▂▁▇▃▂ █▅██▅▆  █ ▇█ ▇ ▄█ ██████▅██▆████▃█ ▃ ▂▅  █      ▄ █▄ █▅██ ▂█ ▃█ ▃█████████▇███ ▇█ ██▇███████▄██▂████ ██████▄██ █ ██▄█▃▆  ▃██████   │ Percent covered: 99.92%
>    13.4 │   ██████████████████████████▅▇█▃██████████▇█████████████████▆▇█▂██▄█ ███████████████████▄█▂██ ██▁▇▆▃▅▇█▅██▇███████▇██▅████████████████████████████████████▅█████████▃████████▂▁███████▁  │ Mean coverage:   18.8x
>    10.7 │ ▁ ██████████████████████████████████████████████████████████████████▆██████████████████████████████████████████████████████████████████████████████████████████████████████████████████  │ Mean baseQ:      17
>     8.0 │ █▆██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▆│ Mean mapQ:       60
>     5.3 │▁█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ 
>     2.7 │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Histo bin width: 160bp
>     0.0 │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Histo max cov:   26.725
          1        1.6K      3.2K      4.8K      6.4K      8.0K      9.6K     11.2K     12.8K     14.4K     16.0K     17.6K     19.2K     20.8K     22.4K     24.0K     25.6K     27.2K           29.8K  

MT192765.1 (29.8Kbp)
>    28.7 │                                  ▃                                                                                                               █                                       │ Number of reads: 7974
>    25.5 │                               ▃  █                                           ▃                                           ▇                       █  ▁             ▃                      │ 
>    22.3 │    ▂        ▁▃ ▂  ▅  █  ▆  ▂  █ ▅██ █▂   ▆ ▄▆ ▁                            ▁▁█▂▂▁  ▂ ▃     ▂                           ▄ █▄▅ ▇               ▁▆ ▄█  █ ▄▅ ▂  ▆  ▄  █▁ ▁    █ ▃   ▁        │ Covered bases:   29.8Kbp
>    19.1 │    █ ▅▄  ▁  ██▆█▅ █  █▂▂█▄▄█ ▃█ ███▂██ ▅ █ ██ █▄                           ██████▇██ █   ▄ █ ▃▄       ▁ ▂ █▁ ██  ▇  ▃▇▁█▄███▄█▆▂▃  ▃ ▁▆ ▃▁▆▁ ██▄██▂██▁██▃█▁ █  █▃ ██ █▆  ██ █   █▃▂▂█    │ Percent covered: 99.96%
>    15.9 │   ▁█▁██  █ ▇█████▃█ ▅██████████▇██████▇█▇███████                          ▅█████████ █▅▃▄█▅█▇██▅  ▆ ▇▂█▂█ ██ ██▅▄█▃▅█████████████▃▆█ ███████ ██████████████▄██▁█████▇██▅███▅█▂ ▆█████▁   │ Mean coverage:   18.7x
>    12.8 │   █████▇▇█▆█████████████████████████████████████                         ▆██████████▆███████████▆███████████████████████████████████▄█████████████████████████████████████████▄███████   │ Mean baseQ:      17
>     9.6 │ ▂▅██████████████████████████████████████████████  ▂   ▆               ▂ ▁█████████████████████████████████████████████████████████████████████████████████████████████████████████████▄▆ │ Mean mapQ:       60
>     6.4 │▄████████████████████████████████████████████████▂▇█▂▆▃█▃█▆▇▅ ▇▇▁▆▇▂█▁▆█▆████████████████████████████████████████████████████████████████████████████████████████████████████████████████▆│ 
>     3.2 │█████████████████████████████████████████████████████████████▇████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Histo bin width: 160bp
>     0.0 │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Histo max cov:   31.9
          1        1.6K      3.2K      4.8K      6.4K      8.0K      9.6K     11.2K     12.8K     14.4K     16.0K     17.6K     19.2K     20.8K     22.4K     24.0K     25.6K     27.2K           29.8K  

MT192765.1 (29.8Kbp)
>    26.5 │                                                                                                                                          █                                               │ Number of reads: 7975
>    23.5 │                         ▃ ▃                                                    ▂                                         █ █   ▄         █                 ▂                             │ 
>    20.6 │    ▃  ▁      █ ▆     ▄ ▃█ █  ▄▅ ▆▆▄ ▁    ▃  █         ▂ ▃▁                   ▅ █ ▂     ▃ ▂              █ ▁   █ ▁     ▄ ▂█▂█▇▆██         █  ▁    ▁  ▄ ▅▇▂▁ █    ▃ ▂    ▂  ▁          ▆   │ Covered bases:   29.8Kbp
>    17.6 │    █ ▄█ ▅ ▁▇▄█▁█  ▆ ▃█▁██ █▄▆██▇███▂█▅ ▇ █ ▃██▄▅▅ ▂█  █ ██ ▁▃  ▁  ▃▂  ▅  █ ▄▄█▁█▁█▆▇█▃ █ █ ▂   ▄▂     ▅▂█ █▆▇▆█▁█▄▄▆  █▁████████▄▇▅   ▇  █▃▅█ ▆▄ █ ▂█ ████ █  ▇▃█ ██ █ █  █        ▅ █   │ Percent covered: 99.95%
>    14.7 │   ██▃██ █▆██████▅▃█ █████▄████████████▇█▅█▅█████████▃▃█ ██ ██ ▁█████ ▃█▁▂█▁███████████▁█▆█▄█▆▃▂██▄▇▄▅▃███▆██████████▂▁█████████████▆ ▆█▅▆████▁██▃█▇██▁████▃█▆▁████████ ██▃█ ▇▇▃▃▃▁██▄█   │ Mean coverage:   18.7x
>    11.8 │   █████████████████▅██████████████████████████████████████▅██▄██████▄█████████████████████████████████████████████████████████████████████████████████████████████████▆████▅██████████▂  │ Mean baseQ:      17
>     8.8 │ ▃▅█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▂▆│ Mean mapQ:       60
>     5.9 │▃█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ 
>     2.9 │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Histo bin width: 160bp
>     0.0 │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Histo max cov:   29.413
          1        1.6K      3.2K      4.8K      6.4K      8.0K      9.6K     11.2K     12.8K     14.4K     16.0K     17.6K     19.2K     20.8K     22.4K     24.0K     25.6K     27.2K           29.8K  

MT192765.1 (29.8Kbp)
>    46.9 │                                                                                                     ▂▆ ▆▅ █▄▁▂▃▂ ▂  ▁                                                                    │ Number of reads: 7979
>    41.7 │                                                                                                    ▄██▇██▅██████▇█▂ █                                                                    │ 
>    36.5 │                                                                                                    ████████████████▃█                                                                    │ Covered bases:   29.8Kbp
>    31.3 │                                                                                                    ███████████████████                                                                   │ Percent covered: 99.95%
>    26.1 │                                                                                                    ███████████████████                                                                   │ Mean coverage:   18.7x
>    20.9 │                                                                                    ▄              ████████████████████ ▁       ▂                                                         │ Mean baseQ:      17
>    15.6 │    ▃ ▄ ▁   ▄ ▇▁▆▅▂▁▁ ▄  ▇ ▂      ▁▁ ▂    ▂  ▄  ▂      ▁ ▂  ▁       ▄ ▄▃ ▅▅ ▅▆▃▁▅▁▁▂█▇▁▆▅▂▂ ▁      ████████████████████ █ ▂▃▄▃▅▄█▂ ▄▃ ▃▂ ▃▅▁▆▁▂▂  ▆ ▁  ▁▅     ▁▂▁█▂ ▆▂  ▂        ▃  ▆▂    │ Mean mapQ:       60
>    10.4 │   ▅███▆█▅▃▅█▅█████████▆▂█▆█▃▆▆▆▅▅██▅█▆█▄▆█ ▆█▁██▅▇▆█▅▂█▆█▁██▃██▃██▇█▅██ ██▆███████████████▄█▅ ▃▆▄█████████████████████ ██████████▅██ ██▅██████████▇██▆██▅█▁▆▇█████▆██▅▄█▄▂▆▃▄▂▅▃█▆▇██▄   │ 
>     5.2 │ ▅▆████████████████████████████████████████▇██████████████████████████████████████████████████▆████████████████████████▆█████████████▇██████████████████████████████████████████████████▃ │ Histo bin width: 160bp
>     0.0 │▇████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████▇│ Histo max cov:   52.138
          1        1.6K      3.2K      4.8K      6.4K      8.0K      9.6K     11.2K     12.8K     14.4K     16.0K     17.6K     19.2K     20.8K     22.4K     24.0K     25.6K     27.2K           29.8K  
```




```make
OUTDIR=out

define run

$(OUTDIR)/SAMPLE$(1).R2.fq: $(OUTDIR)/SAMPLE$(1).R1.fq
	touch $$@

$(OUTDIR)/SAMPLE$(1).R1.fq: $(OUTDIR)/genome$(1).fa $(OUTDIR)/genome$(1).fa.fai
	rm -f $(OUTDIR)/R1$(1).fq $(OUTDIR)/R2$(1).fq
	/LAB-DATA/BiRD/users/lindenbaum-p/packages/samtools/misc/wgsim -N 4000 $$< $(OUTDIR)/SAMPLE$(1).R1.fq $(OUTDIR)/SAMPLE$(1).R2.fq

$(OUTDIR)/SAMPLE$(1).delly.bcf: $(OUTDIR)/delly $(OUTDIR)/SAMPLE$(1).bam $(OUTDIR)/genome_WT.fa $(OUTDIR)/genome_WT.fa.fai
	$(OUTDIR)/delly call -q 0 -r 0  -c 1 -g $(OUTDIR)/genome_WT.fa $(OUTDIR)/SAMPLE$(1).bam |\
		bcftools view -O b -o $$@ && bcftools index -f $$@ 

$(OUTDIR)/SAMPLE$(1).delly.gt.bcf: $(OUTDIR)/delly $(OUTDIR)/SAMPLE$(1).bam $(OUTDIR)/delly.sites.bcf  $(OUTDIR)/genome_WT.fa 
	$(OUTDIR)/delly call -q 0 -r 0  -c 1 -v $(OUTDIR)/delly.sites.bcf  -g $(OUTDIR)/genome_WT.fa $(OUTDIR)/SAMPLE$(1).bam |\
		bcftools view -O b -o $$@ && bcftools index -f $$@ 


$(OUTDIR)/SAMPLE$(1).cov: $(OUTDIR)/SAMPLE$(1).bam
	samtools coverage --plot-depth $$<

endef

SUFF=WT DEL INV DUP

all: coverage $(OUTDIR)/archive.zip
	ls -lah $(OUTDIR)/archive.zip

$(OUTDIR)/archive.zip: $(OUTDIR)/delly.filter.bcf $(OUTDIR)/bcftools.annot.bcf  $(OUTDIR)/MT192765.1.gtf.gz $(OUTDIR)/hapcaller.annot.bcf Makefile
	rm -f $@
	zip -j $@ \
		Makefile \
		$(OUTDIR)/delly.filter.bcf \
		$(OUTDIR)/delly.filter.bcf.csi \
		$(OUTDIR)/hapcaller.annot.bcf \
		$(OUTDIR)/hapcaller.annot.bcf.csi \
		$(OUTDIR)/hapcaller.call.vcf.gz \
		$(OUTDIR)/hapcaller.call.vcf.gz.tbi \
		$(OUTDIR)/bcftools.annot.bcf \
		$(OUTDIR)/bcftools.annot.bcf.csi \
		$(OUTDIR)/bcftools.call.bcf \
		$(OUTDIR)/bcftools.call.bcf.csi \
		$(OUTDIR)/MT192765.1.gff3.gz \
		$(OUTDIR)/MT192765.1.gff3.gz.tbi \
		$(OUTDIR)/MT192765.1.gtf.gz \
		$(OUTDIR)/MT192765.1.gtf.gz.tbi \
		$(OUTDIR)/delly.merge.bcf \
		$(OUTDIR)/delly.merge.bcf.csi \
		$(OUTDIR)/delly.sites.bcf \
		$(OUTDIR)/delly.sites.bcf.csi \
		$(addsuffix .bam,$(addprefix $(OUTDIR)/SAMPLE_, $(SUFF))) \
		$(addsuffix .bam.bai,$(addprefix $(OUTDIR)/SAMPLE_, $(SUFF))) \
		$(addsuffix .delly.bcf,$(addprefix $(OUTDIR)/SAMPLE_, $(SUFF))) \
		$(addsuffix .delly.bcf.csi,$(addprefix $(OUTDIR)/SAMPLE_, $(SUFF))) \
		$(addsuffix .delly.gt.bcf,$(addprefix $(OUTDIR)/SAMPLE_, $(SUFF))) \
		$(addsuffix .delly.gt.bcf.csi,$(addprefix $(OUTDIR)/SAMPLE_, $(SUFF)))

$(OUTDIR)/bcftools.annot.bcf : $(OUTDIR)/bcftools.call.bcf  $(OUTDIR)/MT192765.1.gff3.gz  $(OUTDIR)/genome_WT.fa
	bcftools csq --force -v 2 --local-csq -g $(OUTDIR)/MT192765.1.gff3.gz -f  $(OUTDIR)/genome_WT.fa  $(OUTDIR)/bcftools.call.bcf -O b -o $@ && bcftools index -f $@

$(OUTDIR)/bcftools.call.bcf : $(addsuffix .bam,$(addprefix $(OUTDIR)/SAMPLE_, $(SUFF)))  $(OUTDIR)/genome_WT.fa
	bcftools mpileup -a 'FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/QS,FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/ADF,INFO/ADR,INFO/SCR'  --fasta-ref $(OUTDIR)/genome_WT.fa  $(filter %.bam,$^) |\
		bcftools call -mv -a 'INFO/PV4,FORMAT/GQ,FORMAT/GP' -Ob -o $@ && bcftools index -f $@

$(OUTDIR)/hapcaller.annot.bcf : $(OUTDIR)/hapcaller.call.vcf.gz  $(OUTDIR)/MT192765.1.gff3.gz  $(OUTDIR)/genome_WT.fa
	bcftools csq --force -v 2 --local-csq -g $(OUTDIR)/MT192765.1.gff3.gz -f  $(OUTDIR)/genome_WT.fa  $(OUTDIR)/hapcaller.call.vcf.gz -O b -o $@ && bcftools index -f $@

$(OUTDIR)/hapcaller.call.vcf.gz :  $(addsuffix .bam,$(addprefix $(OUTDIR)/SAMPLE_, $(SUFF)))  $(OUTDIR)/genome_WT.fa $(OUTDIR)/genome_WT.dict 
	gatk HaplotypeCaller --min-base-quality-score 0   -stand-call-conf  1 --base-quality-score-threshold 6 -R $(OUTDIR)/genome_WT.fa $(addprefix -I ,$(filter %.bam,$^) ) -O $@

coverage:  $(addsuffix .cov,$(addprefix $(OUTDIR)/SAMPLE_, $(SUFF)))

$(OUTDIR)/delly.filter.bcf:$(OUTDIR)/delly.merge.bcf
	$(OUTDIR)/delly filter --tag -f germline -o $@ $< && bcftools index -f $@

$(OUTDIR)/delly.merge.bcf : $(addsuffix .delly.gt.bcf,$(addprefix $(OUTDIR)/SAMPLE_, $(SUFF)))
	bcftools merge -m id -O b -o $@ $^ && \
		bcftools index -f $@

$(OUTDIR)/delly.sites.bcf :  $(addsuffix .delly.bcf,$(addprefix $(OUTDIR)/SAMPLE_, $(SUFF)))
	$(OUTDIR)/delly merge -o $@ $^

$(OUTDIR)/delly:
	mkdir -p $(dir $@)
	wget -O $@ "https://github.com/dellytools/delly/releases/download/v1.2.6/delly_v1.2.6_linux_x86_64bit"
	chmod +x $@

$(OUTDIR)/SAMPLE_DUP.bam: $(OUTDIR)/genome_WT.fa.bwt $(OUTDIR)/SAMPLE_DUP.R1.fq $(OUTDIR)/SAMPLE_DUP.R2.fq  $(OUTDIR)/SAMPLE_WT.R1.fq $(OUTDIR)/SAMPLE_WT.R2.fq
	paste <(cat $(OUTDIR)/SAMPLE_DUP.R1.fq $(OUTDIR)/SAMPLE_WT.R1.fq | paste - - - - ) <(cat $(OUTDIR)/SAMPLE_DUP.R2.fq $(OUTDIR)/SAMPLE_WT.R2.fq| paste - - - - ) |\
		shuf |\
		awk '(rand()<0.5)' |\
		tr "\t" "\n" |\
		bwa mem -p  -R '@RG\tID:SAMPLE_DUP\tSM:SAMPLE_DUP\tPL:ILLUMINA' $(basename $< ) -|\
                samtools sort -T $(OUTDIR)/tmp  -O BAM -o $@ && samtools index $@

$(OUTDIR)/SAMPLE_DEL.bam: $(OUTDIR)/genome_WT.fa.bwt $(OUTDIR)/SAMPLE_DEL.R1.fq $(OUTDIR)/SAMPLE_DEL.R2.fq  $(OUTDIR)/SAMPLE_WT.R1.fq $(OUTDIR)/SAMPLE_WT.R2.fq
	paste <(cat $(OUTDIR)/SAMPLE_DEL.R1.fq $(OUTDIR)/SAMPLE_WT.R1.fq | paste - - - - ) <(cat $(OUTDIR)/SAMPLE_DEL.R2.fq $(OUTDIR)/SAMPLE_WT.R2.fq| paste - - - - ) |\
		shuf |\
		awk '(rand()<0.5)' |\
		tr "\t" "\n" |\
		bwa mem -p  -R '@RG\tID:SAMPLE_DEL\tSM:SAMPLE_DEL\tPL:ILLUMINA' $(basename $< ) -|\
                samtools sort -T $(OUTDIR)/tmp  -O BAM -o $@ && samtools index $@

$(OUTDIR)/SAMPLE_INV.bam: $(OUTDIR)/genome_WT.fa.bwt $(OUTDIR)/SAMPLE_INV.R1.fq $(OUTDIR)/SAMPLE_INV.R2.fq  $(OUTDIR)/SAMPLE_WT.R1.fq $(OUTDIR)/SAMPLE_WT.R2.fq
	paste <(cat $(OUTDIR)/SAMPLE_INV.R1.fq $(OUTDIR)/SAMPLE_WT.R1.fq | paste - - - - ) <(cat $(OUTDIR)/SAMPLE_INV.R2.fq $(OUTDIR)/SAMPLE_WT.R2.fq| paste - - - - ) |\
		shuf |\
		awk '(rand()<0.5)' |\
		tr "\t" "\n" |\
		bwa mem -p  -R '@RG\tID:SAMPLE_INV\tSM:SAMPLE_INV\tPL:ILLUMINA' $(basename $< ) -|\
                samtools sort -T $(OUTDIR)/tmp  -O BAM -o $@ && samtools index $@

$(OUTDIR)/SAMPLE_WT.bam: $(OUTDIR)/genome_WT.fa.bwt $(OUTDIR)/SAMPLE_WT.R1.fq $(OUTDIR)/SAMPLE_WT.R2.fq  $(OUTDIR)/genome_WT.fa.fai
	bwa mem  -R '@RG\tID:SAMPLE_WT\tSM:SAMPLE_WT\tPL:ILLUMINA' $(basename $< ) $(word 2,$^) $(word 3,$^) |\
		samtools sort -T $(OUTDIR)/tmp  -O BAM -o $@ && samtools index $@

$(eval $(call run,_DEL))
$(eval $(call run,_DUP))
$(eval $(call run,_WT))
$(eval $(call run,_INV))

$(OUTDIR)/genome_WT.fa.bwt: $(OUTDIR)/genome_WT.fa
	bwa index $<

$(OUTDIR)/genome_DUP.fa.fai: $(OUTDIR)/genome_DUP.fa
	samtools faidx $<

$(OUTDIR)/genome_DUP.fa: $(OUTDIR)/genome_WT.fa
	awk '(NR<200)' $< > $@
	awk '(NR>200 && NR<240)' $< >> $@
	awk '(NR>200 && NR<240)' $< >> $@
	awk '(NR>200 && NR<240)' $< >> $@
	awk '(NR>200 && NR<240)' $< >> $@
	awk '(NR>200 && NR<240)' $< >> $@
	awk '(NR>200)' $< >> $@

$(OUTDIR)/genome_DEL.fa.fai: $(OUTDIR)/genome_DEL.fa
	samtools faidx $<

$(OUTDIR)/genome_DEL.fa: $(OUTDIR)/genome_WT.fa
	awk '(NR<100)' $< > $@
	awk '(NR>150)' $< >> $@

$(OUTDIR)/genome_INV.fa.fai: $(OUTDIR)/genome_INV.fa
	samtools faidx $<

$(OUTDIR)/genome_INV.fa: $(OUTDIR)/genome_WT.fa
	awk '(NR<350)' $< > $@
	awk '(NR>=350 && NR<360)' $< | paste -s -d "" | tr "ATGC" "TACG" | rev | fold -w 80 >> $@
	awk '(NR>=360)' $< >> $@

$(OUTDIR)/genome_WT.fa.fai: $(OUTDIR)/genome_WT.fa
	wget -O $@ "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta.fai"

$(OUTDIR)/genome_WT.dict : $(OUTDIR)/genome_WT.fa
	wget -O $@ "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.dict"

$(OUTDIR)/genome_WT.fa:
	mkdir -p $(dir $@)
	wget -O $@ "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta"

$(OUTDIR)/MT192765.1.gtf.gz : MT192765.1.gff3
	agat_convert_sp_gff2gtf.pl --gff $< -o $(addsuffix .tmp,$@)
	LC_ALL=C sort -t $$'\t' -k1,1 -k4,4n  $(addsuffix .tmp,$@) | bgzip > $@ && tabix -p gff -f $@
	rm  $(addsuffix .tmp,$@)
	
$(OUTDIR)/MT192765.1.gff3.gz : MT192765.1.gff3
	LC_ALL=C sort -t $$'\t' -k1,1 -k4,4n $< | bgzip > $@ && tabix -p gff -f $@

workdir:
	@echo $(OUTDIR)
clean:
	rm -rvf "$(OUTDIR)/"

```



