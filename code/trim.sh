#!/bin/bash

# remove 5` gene primers and degenerate heterogeneity spacers from forward and reverse paired-end sequences
# using cutadapt v1.18 https://github.com/marcelm/cutadapt and SeqPurge 2019_05 (part of https://github.com/imgag/ngs-bits)


#input - folder in output/demux
libPath="output/demux/"

#make folders for trimming output
#these folders were changed to match denoise.R code 2/23/22 - may not work next time around so check
outf1="output/trim/ITS1/fungi"
mkdir -p $outf1
outf2="output/trim/ITS2/fungi"
mkdir -p $outf2
outb="output/trim/16S/bact"
mkdir -p $outb

#make directory for scratch data
scratch="output/scratch/"
mkdir -p $scratch 

# define gene primers for trimming
fwd_primer_f1="CTHGGTCATTTAGAGGAASTAA" #ITS1F-KYO1
fwd_primer_rc_f1="TTASTTCCTCTAAATGACCDAG"
rev_primer_f1="TTYRCTRCGTTCTTCATC" #ITS2-KYO2
rev_primer_rc_f1="GATGAAGAACGYAGYRAA"

fwd_primer_f2="TCCTCCGCTTATTGATATGC" #ITS4
fwd_primer_rc_f2="GCATATCAATAAGCGGAGGA"
rev_primer_f2="CAHCGATGAAGAACRYAG" #ITS3_kyo1
rev_primer_rc_f2="CTRYGTTCTTCATCGDTG"

fwd_primer_b="GTGYCAGCMGCCGCGGTAA" #515f
fwd_primer_rc_b="TTACCGCGGCKGCTGRCAC" 
rev_primer_b="GGACTACNVGGGTWTCTAAT" #806r
rev_primer_rc_b="ATTAGAWACCCBNGTAGTCC"

#Make file for trimming summary output
log="output/trim/summary.tab"
echo -e "Sample\tRaw\tf1.trim1\tf1.trim2\tf2.trim1\tf2.trim2\tb.trim1\tb.trim2" > ${log}

#Loop over samples

for fwd in $( find $libPath -name "*R1.fastq.gz" | grep -v "undetermined"); do
    #fwd=$( find $libPath -name "P01_01_A.R1.fastq.gz")
    rev=`echo $fwd | sed 's/R1.fastq.gz/R2.fastq.gz/'`
    samp=`echo $fwd | awk -F'/' '{print $3}'| awk -F'.' '{print $1}'`
    
    ############
    ### ITS1 ###
    ############

        ## suggested practice ##
        #1. run cutadapt without the --overlap (the length of the primer), --minimum-length, and --maximum-length flags on a test sample
        #2. observe the output for each sample, noting the amount of error, minimum, and maximum read lengths
        #3. use these to determine the appropriateness of each relevant flag

        #Ideally, a priori primer length and amplicon length (301 in this case) are all that is needed to provide values for these flags, but this may need to be adjusted
        #After a run through with a lowish maximum accepted error rate [-e] and [-gt] 0 at the "trim overhang" step, check the summary file to see if:

        #(a) [-e] needs to be increased to a higher value, which is better for larger libraries where it would be difficult to do (b)
        #(b) [-gt] needs to be increased to the great observed erroneous sequence occurrence
        #           e.g. samples expected to only have ITS2 are found to have trace (~25) amounts of trimmed ITS1 sequences, so -[gt] set to 25

    
    # trim primers in paired-end mode
    trim1f1R1="${scratch}${samp}.f1.R1.trim1.fq.gz"
    trim1f1R2="${scratch}${samp}.f1.R2.trim1.fq.gz"
    
    # important cutadapt setting: overlap = length of shorter primer, min/max = set to global expectations, e = error rate
    cutadapt --quiet -g $fwd_primer_f1 -G $rev_primer_f1 --discard-untrimmed --overlap 18 -e 0.20 --minimum-length 220 --maximum-length 230 -o $trim1f1R1 -p $trim1f1R2 $fwd $rev
    
    #Count seqs after first trim step (the "@M01498" should match the first 5+ characters of your fastq headers)
    if [[ -f $trim1f1R1 && -f $trim1f1R2 ]]; then
        f1_mid=`gzip -cd $trim1f1R1 | grep -c '^@M01498'`
    else
        f1_mid=0
    fi

    #trim any remaining primer overhang
    if [ $f1_mid -gt 25 ]; then
        trim2f1R1="${scratch}${samp}.f1.R1.trim2.fq.gz"
        trim2f1R2="${scratch}${samp}.f1.R2.trim2.fq.gz"
        SeqPurge -in1 $trim1f1R1 -in2 $trim1f1R2 -out1 $trim2f1R1 -out2 $trim2f1R2 -a1 $rev_primer_rc_f1 -a2 $fwd_primer_rc_f1 -qcut 0 -ncut 0 -min_len 100 -summary ${scratch}${samp}.f1.log.txt
    fi

    #Count seqs after final step
    if [[ -f $trim2f1R1 && -f $trim2f1R2 ]]; then
        f1_final=`gzip -cd $trim2f1R1 | grep -c '^@M01498'`
    else
        f1_final=0
    fi
    #Keep only samples that made it through the processing (also truncate reads)
    if [ $f1_final -gt 25 ]; then
        outf1R1="${outf1}${samp}.f1.R1.fq.gz"
        outf1R2="${outf1}${samp}.f1.R2.fq.gz"
        cutadapt --quiet -g XX -G XX --length 220 -o $outf1R1 -p $outf1R2 $trim2f1R1 $trim2f1R2
    fi

    ############
    ### ITS2 ###
    ############
    # trim primers in paired-end mode
    trim1f2R1="${scratch}${samp}.f2.R1.trim1.fq.gz"
    trim1f2R2="${scratch}${samp}.f2.R2.trim1.fq.gz"


    cutadapt --quiet -g $fwd_primer_f2 -G $rev_primer_f2 --discard-untrimmed --overlap 18 -e 0.20 --minimum-length 224 --maximum-length 230 -o $trim1f2R1 -p $trim1f2R2 $fwd $rev

    #Count seqs after first trim step (the "@M01498" should match the first 5+ characters of your fastq headers)
    if [[ -f $trim1f2R1 && -f $trim1f2R2 ]]; then
        f2_mid=`gzip -cd $trim1f2R1 | grep -c '^@M01498'`
    else
        f2_mid=0
    fi

    #trim any remaining primer overhang
    if [ $f2_mid -gt 25 ]; then
        trim2f2R1="${scratch}${samp}.f2.R1.trim2.fq.gz"
        trim2f2R2="${scratch}${samp}.f2.R2.trim2.fq.gz"
        SeqPurge -in1 $trim1f2R1 -in2 $trim1f2R2 -out1 $trim2f2R1 -out2 $trim2f2R2 -a1 $rev_primer_rc_f2 -a2 $fwd_primer_rc_f2 -qcut 0 -ncut 0 -min_len 100 -summary ${scratch}${samp}.f2.log.txt
    fi

    #Count seqs after final step
    if [[ -f $trim2f2R1 && -f $trim2f2R2 ]]; then
        f2_final=`gzip -cd $trim2f2R1 | grep -c '^@M01498'`
    else
        f2_final=0
    fi
    #Keep only samples that made it through the processing (also truncate reads)
    if [ $f2_final -gt 25 ]; then
        outf2R1="${outf2}${samp}.f2.R1.fq.gz"
        outf2R2="${outf2}${samp}.f2.R2.fq.gz"
        cutadapt --quiet -g XX -G XX --length 224 -o $outf2R1 -p $outf2R2 $trim2f2R1 $trim2f2R2
    fi

    ############
    ### BACT ###
    ############

 # trim primers in paired-end mode
    trim1bR1="${scratch}${samp}.b.R1.trim1.fq.gz"
    trim1bR2="${scratch}${samp}.b.R2.trim1.fq.gz"


    cutadapt --quiet -g $fwd_primer_b -G $rev_primer_b --discard-untrimmed --overlap 19 -e 0.20 --minimum-length 223 --maximum-length 230 -o $trim1bR1 -p $trim1bR2 $fwd $rev

    #Count seqs after first trim step (the "@M01498" should match the first 5+ characters of your fastq headers)
    if [[ -f $trim1bR1 && -f $trim1bR2 ]]; then
        b_mid=`gzip -cd $trim1bR1 | grep -c '^@M01498'`
    else
        b_mid=0
    fi

    #trim any remaining primer overhang
    if [ $b_mid -gt 25 ]; then
        trim2bR1="${scratch}${samp}.b.R1.trim2.fq.gz"
        trim2bR2="${scratch}${samp}.b.R2.trim2.fq.gz"
        SeqPurge -in1 $trim1bR1 -in2 $trim1bR2 -out1 $trim2bR1 -out2 $trim2bR2 -a1 $rev_primer_rc_b -a2 $fwd_primer_rc_b -qcut 0 -ncut 0 -min_len 100 -summary ${scratch}${samp}.b.log.txt
    fi

    #Count seqs after final step
    if [[ -f $trim2bR1 && -f $trim2bR2 ]]; then
        b_final=`gzip -cd $trim2bR1 | grep -c '^@M01498'`
    else
        b_final=0
    fi
    #Keep only samples that made it through the processing (also truncate reads)
    if [ $b_final -gt 25 ]; then
        outbR1="${outb}${samp}.b.R1.fq.gz"
        outbR2="${outb}${samp}.b.R2.fq.gz"
        cutadapt --quiet -g XX -G XX --length 223 -o $outbR1 -p $outbR2 $trim2bR1 $trim2bR2
    fi

###end
    #make summary and append to file
    RAW=`gzip -cd $fwd | grep -c '^@M01498'`
    echo -e $samp"\t"$RAW"\t"$f1_mid"\t"$f1_final"\t"$f2_mid"\t"$f2_final"\t"$b_mid"\t"$b_final | cat >> ${log}

    #remove tmp files
    rm ${scratch}${samp}.*
done

rm -r ${scratch}


