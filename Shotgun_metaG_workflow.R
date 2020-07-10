######################Sequence Quality Control Workflow####################

#First Created: April 18th 2018
#Author: RMJ
#Last Edited: April 19th 2018

#This code is run in a linux environment utilizing the tools found in the bbtools toolkit: https://jgi.doe.gov/data-and-tools/bbtools/

#This code has three primary steps 
#1. Trimming
#2. Filtering
#3. Deduplicating

#The code is natively run as a one line script that pipes the output into the next command. For the purposes of this explanation, the parts will be addressed separately. Please see line [16] for the code in its entirety

for bob in *R1_001.fastq.gz; 
do bbduk.sh in=$bob out=trim/$bob.trim.fastq qtrim=rl trimq=10 &>>trim/$bob.trim.stats;
bbduk.sh in=trim/$bob.trim.fastq out=filtered/$bob.MinQ10.trim.fastq maq=10 &>>filtered/$bob.filtered.stats;
dedupe.sh in=filtered/$bob.MinQ10.trim.fastq out=deduped/$bob.dd.MinQ10.trim.fastq ac=f &>>deduped/$bob.deduped.stats; 
done

###bob explained

#this code is a for loop. This allows the code to sequentially run a file through the process and then automatically do the next file, until all files are complete. In this case i have told the code to look for files that end with 'R1_001.fastq.gz' the * indicates that anything can come before the words I specified but it MUST end in 'R1_001.fastq.gz'. The code then makes a list of all the files in the current directory that satisfy this requirement. It will then replace every instance of 'bob' in the code with a filename in this list and workthrough them sequentially until all files have been run. 

#the choice of bob is meaningless. it just has to be a character or string that is unique, not a command, and does not occur elsewhere in the code except for where we want it. For example, calling our variable 'trim' instead of 'bob' would be a poor choice because trim is a command and appears elsewhere in the code so the computer would get confused. 

####### 1.  Trimming ###########

#This first step will perform quality trimming. The program will assess the quality of the sequence and look for regions or instances of low quality (<q10) and trim them. If there are too many bases below q10 across the sequence, such that the sequence length after trimming is <=10 the sequence is removed. This will also look for and remove adapters based on internal reference.

#It usually prints to the terminal screen the output but instead the '&>>' tells it to create a stats file with the name of the input file

for bob in *R1_001.fastq.gz; 

do bbduk.sh in=$bob out=trim/$bob.trim.fastq qtrim=rl trimq=10 ref=/bbmap/resources/adapters.fa &>>trim/$bob.trim.stats;

# in=input sequences
# out=output file
# qtrim=rl (right and left)
# trimq=10 (phred 10)
# ref=adapter ref

########## 2. Filtering ###############

#This next part of the code takes the trimmed sequences and then further filters them based on the average sequence quality. The sequence must have an average quality of >10 in order to pass.

#Similar to the trim code, it would normally output to the terminal screen however, we tell it to store the stats as a file with the name of the input sequences 

bbduk.sh in=trim/$bob.trim.fastq out=filtered/$bob.MinQ10.trim.fastq maq=10 &>>filtered/$bob.filtered.stats; 

# in=input sequences
# out=output file
# maq=average quality score threshold

########### 3. Deduplicating ################

#to reduce the computational load and potentially over amplified reads, we deduplicate exact sequence matches to one representative.

#This code takes the output from the filtered and trimmed files and deduplicates them

#Again it will save the stats as a file with the name of the input file as well

dedupe.sh in=filtered/$bob.MinQ10.trim.fastq out=deduped/$bob.dd.MinQ10.trim.fastq ac=f &>>deduped/$bob.deduped.stats; 

# in=input sequences
# out=output file
# ac=containment absorption (if true would allow for mismatches)

done;

###################################### Diamond Protein Translation and Alignment #####################

#Diamond leverages BLAST and blast databases to align translated DNA sequences to protein sequences within the database allowing us to find a taxonomic and functional match. This is a looped code that functions in a similar way to the QC code but instead of 'bob' it uses 'pups' 

#again the word is meaningless and only stands in to ensure that there is no ambiguity within the code. Here is the code in its entirety and explained below

#Diamond ref: Buchfink, B., Xie, C. and Huson, D.H., 2015. Fast and sensitive protein alignment using DIAMOND. Nature methods, 12(1), p.59.

for pups in *.fastq.gz.dd.MinQ10.trim.fastq; 
do diamond blastx -d ~/Diamond_Tools/BLAST_DB/nr_4_23_2018.dmnd -q $pups -o $pups.daa -f 100 --sallseqid --top 5 -b 8 -e .0001; 
done;

#diamond = engages the diamond program
# blastx = runs it in blastx mode i.e. translate DNA into a query first
# -q = query file, the file containing sequences you want matches for
# -o = output file, where the output of the diamond process will be stored
# -f = determines the output format. 100 codes it to DAA (diamond proprietary format) which can be read by MEGAN
# --sallseqid = include the subject (database sequence) id
# --top = only return alignments whose score is at most 5% lower than the top alignment score. 
# -b = How big to make the chunks of file to be processed as a control for processing time and memory usage. measured in GB x6 i.e. -b 8 = 24GB
# -e = the maximum expected value to report an alignment. It will not report an alignment with an evalue greater than this value. Just because it is reported does not mean that it is the best match. 

#### Defaults not mentioned in code #########

#minimum open reading frame of 20 Amino acids
#scoring matrix=BLOSUM62 
#11/1 gap open/gap extend penalty







