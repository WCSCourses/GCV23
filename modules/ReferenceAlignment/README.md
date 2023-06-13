# Reference Alignment Practical

## Course Details:

* [GCV2023: Genomics and Clinical Virology 2023](https://coursesandconferences.wellcomeconnectingscience.org/event/genomics-and-clinical-virology-20230611/)
* 11th-16th June 2023
* Wellcome Genome Campus, UK
* [https://github.com/WCSCourses/GCV23](https://github.com/WCSCourses/GCV23)

## Contact Details:

[Dr. Richard Orton](https://www.gla.ac.uk/researchinstitutes/iii/staff/richardorton/)  
[Medical Research Council– University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/)  
464 Bearsden Road  
Glasgow  
G61 1QH

E-mail: Richard.Orton@glasgow.ac.uk

## Content:

This practical is associated with a VirtualBox image containing the tools and data sets pre-installed, and an accompnaying lecture on Reference Alignment of High-Throughoput Sequencing (HTS) reads to a reference a sequence.

* [0: Overview](#0-overview)
* [1: Setup](#0-setup)
	+ [1.1: Basic read statistics](#11-basic-read-statistics)
* [2: Read Alignment](#2-read-alignment)
	+ [2.1: Indexing the reference sequence](#21-indexing-the-reference-sequence)
	+ [2.2: Aligning the reads to the reference](#22-aligning-the-reads-to-the-reference])
	+ [2.3: Convert the SAM to a BAM file](#23-convert-the-SAM-to-a-BAM-file)
	+ [2.4: Basic alignment statistics](#24-basic-alignment-statistics)
* [3: Alignment on your own](#3-alignment-on-your-own)
* [4: Extra data](#4-extra-data)
* [5: Assembly Visualisation and Statistics Practical](#5-assembly-visualisation-and-statistics-practical)
	+ [5.1: Setup](#51-setip)
	+ [5.2: Summary Statistics - weeSAM](#52-summary-statistics-weeSAM)
	+ [5.3: Coverage plot on your own](#53-coverage-plot-on-your-own)
	+ [5.4: Visualisation with Tablet](#54-visualisation-with-tablet)
 
# 0: Overview

**YOU DO NOT NEED TO ENTER THE COMMANDS IN THIS OVERVIEW SECTION!**

In this practical, we will be aligning paired end reads to a reference sequence. Commands that you need to enter into the Terminal window (i.e. the command line) are presented in a box and with a different font, like this:

```
ls
```

Sometimes a command is long and doesn’t fit on a single line on the screen (and screen sizes vary), but it should still be entered as one single line on the computer terminal. 

```
bwa mem -t 14 my_reference_file.fasta my_read_file_1.fastq my_read_file_2.fastq > my_output_file.sam
```

A few Linux tips to remember:

1.	Use the **Tab button** to automatically complete filenames – especially long ones
2.	Use the **Up Arrow** to scroll through your previous commands, it enables you to easily re-run or re-use/change/correct old commands
3.	**Case matters**, the following file names are all different:

```
Myfile.txt
MyFile.txt
MYFILE.txt
myfile.txt
my file.txt
my_file.txt
```

Also watch out for number 1s being confused with lowercase letter L’s, and capital O’s being confused with zeroes

```
l = lower case letter L
1 = number one
O = capital letter O
0 = zero
```

# 1: Setup

In this session, we will be using a set of Illumina paired end reads which were simulated from a hepatitis c virus (HCV) genome; these simulated reads were created using ART (Huang et al., 2012: [10.1093/bioinformatics/btr708](10.1093/bioinformatics/btr708)). The goal now is to align these reads to a reference genome sequence, with an ultimate goal of creating a consensus sequence for mutation anlysis.

To start off, we will need to move into the correct folder:

```
cd ~/Richard/HCV/
```

***Command breakdown:***

* **cd** = change directory 

Next, list the contents of the directory so you can see the files we will be working with:

```
ls
```

You should see the FASTQ paired end read files:

**hcv\_sim\_R1.fq**  
**hcv\_sim\_R2.fq**

And also two FASTA reference sequence files:

**1b\_hcv\_ref.fasta**  
**1a\_hcv\_ref.fasta**

We will be aligning the paired end reads to the two reference sequences in turn. The reference sequences represent the 1a and 1b subtypes of HCV, and we will use the alignment results to determine which subtype the sample contains, and to also highlight the importance of selecting an appropriate reference.
 
## 1.1: Basic read statistics

We will first use a tool called prinseq to count the number of reads in each file. As these are paired end reads, there should be one read from each read pair in each file – and hence the same number of reads in each file. We will also use prinseq to output statistics on the read lengths

```
prinseq-lite.pl -stats_info -stats_len -fastq hcv_sim_R1.fq -fastq2 hcv_sim_R2.fq
```

***Command breakdown:***

1.	**prinseq-lite.pl** is the name of the program
2.	**-stats\_info** tells prinseq to output basic stats on the reads (number of reads and bases)
3.	**-stats\_len** tells prinseq to output basic stats on read lengths (min, max, mean etc)
4.	**-fastq hcv\_sim\_R1.fq** the name of the 1st FASTQ file
5.	**-fastq2 hcv\_sim\_R2.fq** the name of the 2nd FASTQ file in the pair

### Common Issue
* A common issue here is not entering the prinseq command on one line in the terminal - you should only use the enter key at the end of the command to execute it.
* Another common issue is typos - check the command carefully if you get an error - it is likely you have a typo

***
### Questions
**Question 1** – How many reads and bases are in the read files 1 and 2?

**Question 2** – What is the average (mean) length of the reads? 
***

The statistics are split into those for the first FASTQ file of the read pair (e.g. stats\_info, stats\_len, etc) and those for the second FASTQ file of the read pair (e.g. stats\_info2, stats\_len2, etc), and should look a like this

```
stats_info	bases		48000000
stats_info	reads		320000
stats_info2	bases		48000000
stats_info2	reads		320000
stats_len	max		150
stats_len	mean		150.00
stats_len	median		150
stats_len	min		150
stats_len	mode		150
stats_len	modeval		320000
stats_len	range		1
stats_len	stddev		0.00
stats_len2	max		150
stats_len2	mean		150.00
stats_len2	median		150
stats_len2	min		150
stats_len2	mode		150
stats_len2	modeval		320000
stats_len2	range		1
stats_len2	stddev		0.00 
```

Paired read files should always have the same number of lines/reads (the ordering of the reads in each file is also critical), so if your two paired files have a different number of lines, something has gone wrong (e.g. filtering/trimming went wrong and corrupted the output, or maybe files from different samples are being used). 
 
# 2: Read Alignment

There are many tools available to align reads onto a reference sequence: bwa, bowtie2, minimap2, bbMap, to name but a few.

We will be using [BWA](http://bio-bwa.sourceforge.net) to align our paired end reads to a reference sequence and output a [SAM (Sequence Alignment Map)](https://samtools.github.io/hts-specs/SAMv1.pdf) file. The SAM file contains the result of each read’s alignment to the given reference sequence. 

## 2.1: Indexing the reference sequence

First, we need to create a BWA index of the reference sequence. Tools such as BWA need to index the sequence first to create a fast lookup (or index) of short sequence seeds within the reference sequence. This enables the tools to rapidly align millions of reads:

```
bwa index 1b_hcv_ref.fasta
```

If you list (ls) the contents of the Refs directory, you should see the BWA index files, they will all have the prefix 1b\_hcv\_ref.fasta, and will have extensions such as **.amb**, **.ann**, **.bwt**, **.pac**, and **.sa**.

```
ls
```

## 2.2: Aligning the reads to the reference

Next, we want to align our reads to the reference sequence using the BWA mem algorithm:

```
bwa mem -t 4 1b_hcv_ref.fasta hcv_sim_R1.fq hcv_sim_R2.fq > 1b.sam
```

***Command breakdown:***

1. **bwa** = the name of the program we are executing
2. **mem** = the BWA algorithm to use (recommended for illumina reads > 70nt)
3. **-t 4** = use 4 computer threads
4. **1b\_hcv\_ref.fasta** = the name (and location) of the reference genome to align to
5. **hcv\_sim\_R1.fq** = the name of read file 1
6. **hcv\_sim\_R2.fq** = the name of read file 2
7. **>** = direct the output into a file
8. **1b.sam** = the name of the output SAM file to create 
 
Typically, a SAM file contains a single line for each read alignment in the data set, and this line stores a single alignment result (reference name, alignment location, CIGAR string, the read sequence itself, read quality, etc). 

Overall, this will create an output file called 1b.sam in the current directory, which contains the results (in SAM format) of aligning all our reads to the reference sequence 1b\_hcv\_ref.fasta.

When bwa has finished (and your prompt comes back), check that the SAM file has been created.

```
ls
```

There should now be a file called 1b.sam in the directory.

### Common issue
A common mistake is not waiting for your previous command to finish, and entering the next command into the terminal before the prompt has returned. You need to wait until the **manager@GCV2023** command prompt returns before entering the next command.

## 2.3: Convert the SAM to a BAM file

Typically, a SAM file contains a single line for each read in the data set, and this line stores the alignment result of each read (reference name, alignment location, CIGAR string, the read sequence itself, quality, etc).

SAM files are in a text format (which you can open and view if you like: head 1b.sam), but can take up a lot of disk storage space. It is good practice to convert your SAM files to BAM (Binary Alignment Map) files, which are compressed binary versions of the same data, and can be sorted and indexed easily to make searches faster. We will use samtools to convert our SAM to BAM, and sort and index the BAM file:

```
samtools sort 1b.sam -o 1b.bam
```

```
samtools index 1b.bam
```

***Command breakdown:***

1.	The first command tells samtools to **sort** the SAM file, and to also output (**-o**)the sorted data in BAM format to a file called **1b.bam**
2.	We then use samtools to **index** the BAM file 1b.bam (indexing [which relies on sorted data] enables faster searches downstream).


There should now be two new files in the directory called: 

**1b.bam** (the BAM file)  
**1b.bam.bai** (the BAM index file) 

Now let’s list (ls) the contents of the directory to check we have our new files, and also check out their sizes:

```
ls -lh
```

***Command breakdown:***
-l tells the list command to give the output in a long list format, whilst the h tells it to provide file sizes in a human readable format, this is the 5th column, which will have the size of each file in a format such as 2.5M (M for megabytes) or 9.5G (G for gigabytes).

***

### Questions
**Question 3** – How big is the SAM file compared to the BAM file?

***

**NB:** If your SAM file is 0B (i.e. 0 bytes, empty) then something went wrong with the bwa alignment step, so restart from there. If you SAM file is fine (i.e. >0), but your BAM file is 0B (i.e. empty), then something went wrong with your SAM to BAM conversion so re-do that step. 

We don’t need our original SAM file anymore (as we have the BAM file now) so we remove (rm) the SAM file 1b.sam:

```
rm 1b.sam
```

## 2.4: Basic alignment statistics

[Samtools](https://samtools.github.io) is one of the key pieces of software in analyses involving High Throughput Sequencing (HTS) data, it has a wide range of functions and combines easily with related tools such as vcftools (for calling variants). One common thing to check is how many reads have aligned to the reference (called mapped), and how many did not (called unmapped). Samtools can report this for us easily, utilising the aligner SAM flags you learnt about in the previous session.

Reminder: the 2nd column in the SAM file contains the flag for the read alignment. If the flag includes the number 4 flag in its makeup then the read is unmapped, if it doesn’t include the number 4 in it's makeup
then it is mapped.

### Number of mapped read alignments:
```
samtools view -c -F4 1b.bam
```

***Command breakdown***

1.	**samtools view** = to view the file 1b.bam
2.	**–c** = count the read alignments
3.	**–F4** = skip read alignments that contain the unmapped Flag 4 

### Number of unmapped reads
```
samtools view -c -f4 1b.bam
```

This time we use –f4, only include read alignments that do have the unmapped flag 4

***
### Questions

**Question 4** – how many reads are mapped to the 1b_hcv_ref.fasta genome?

**Question 5** – how many reads are unmapped?
***

Technically, the above command gives the number of mapped read alignments not reads. A read could be mapped equally well to multiple positions (one will be called the primary alignment, and others secondary alignments [sam flag 256]), or a read could be split into two parts (e.g. spliced) with one part being the primary alignment and the others supplementary [sam flag 2048]

So to get the true number of mapped reads you need to count only the alignments that do not have flags 4 (unmapped), 256 (not primary), and 2048 (supplementary) = 4 + 256 + 2048 = 2308

### Number of mapped read reads
```
samtools view -c -F2308 1b.bam
```
or

```
samtools view -c -F4 -F256 -F2048 1b.bam
```

For small RNA viruses, secondary and supplementary alignments tend to be rare. 

# 3: Alignment on your own

You now need to use bwa to align the reads to the 1a_hcv_ref.fasta reference sequence – later in the visualisation and summary statistics section we will be comparing the 1a vs 1b alignment results.

You need to work out the commands yourself based on the previous commands for the 1b_hcv_ref.fasta reference. 

Here is a reminder of the commands you used for 1b HCV which you will need to adapt. 

**NB:** Essentially, you will want to change the reference name in the bwa command, and all of the SAM/BAM filenames in the bwa and samtools commands from 1b to 1a.

```
bwa index 1b_hcv_ref.fasta
```
```
bwa mem -t 4 1b_hcv_ref.fasta hcv_sim_R1.fq hcv_sim_R2.fq > 1b.sam
```
```
samtools sort 1b.sam -o 1b.bam
```
```
samtools index 1b.bam
```
```
rm 1b.sam
```
```
samtools view -c -f4 1b.bam
```
```
samtools view -c -F2308 1b.bam
```

***
### Questions

**Question 6** – how many reads are mapped to the 1a_hcv_ref.fasta genome?

**Question 7** – how many reads are unmapped?

**Question 8** – which reference assembly has the most mapped reads: 1a\_hcv or 1b\_hcv? Therefore, which reference sequence is better (1a or 1b)?
***

# 4: Extra Data

If you are looking for something extra to do, there are additional data sets located in the folder:

### ~/Training/Richard/Ebola/

You will find a set of (gzipped) FASTQ paired end read files, and a reference FASTA sequence to align them to.

The reads are from a patient from the ebola epidemic in West Africa 2014 {Gire et al, 2014} [https://www.ncbi.nlm.nih.gov/pubmed/25214632](https://www.ncbi.nlm.nih.gov/pubmed/25214632)

The reference ebola sequence is from a 2007 outbreak in Democratic Republic of Congo. 

Try aligning the reads to the reference yourself.

### ~/Training/Richard/Noisey/

This is a real HCV sample, but the read quality is quite poor making it quite noisey. Again, two HCV ref sequences are supplied (HCV_1a and HCV_1B). Align the paired end reads to each reference and determine what subtype the sample is by comparing mapping and coverage statistics.

### ~/Training/Richard/Mystery/

This is a mystery sample, combine all the given references sequences into one file using the “cat” command, align the reads to that combined reference and then determine what the virus in the sample is.
 
# 5: Assembly Visualisation and Statistics Practical

In this practical, we will be checking our reference assembly from the previous session. We will use tools to generate summary statistics of the depth and breadth of the coverage across the genome, coverage plots, and visualisation of our assembly using tools such as Tablet and weeSAM. Later sessions of the course will cover how to call the consensus sequence and variants.

## 5.1: Setup

In the previous session, you should have tanoti aligned the paired reads onto two different HCV genomes (types 1a and 1b).

This should have resulted in two BAM files in your HCV folder, lets check:

```
cd ~/Richard/HCV/
```

```
ls
```

You should see (amongst others):

**1a.bam**  
**1a.bam.bai**  
**1b.bam**  
**1b.bam.bai**  

Along with the two reference sequences:

**1a\_hcv\_ref.fasta**  
**1b\_hcv\_ref.fasta**  

We need all these files to proceed, so if you don’t have them – ask for help and we can copy across pre-computed versions.

## 5.2: Summary Statistics - weeSAM

We previously used samtools to count the number of mapped and unmapped reads (using the command samtools view -c -F4 1a.bam), which suggested that HCV 1a was a better reference sequence for our sample based on a greater number of mapped reads, but let’s explore this is more detail using a tool called weeSAM: https://github.com/centre-for-virus-research/weeSAM

weeSAM analyses a SAM or BAM file, generates a graphical coverage plot, and reports a range of summary statistics such as:

* **Ref_Name**: The identifier of the reference.
* **Ref_Len**: The length in bases of each reference.
* **Mapped\_Reads**: Number of reads mapped to each reference.
* **Breadth**: The number of sites in the genome covered by reads.
* **%\_Covered**: The percent of sites in the genome which have coverage.
* **Min\_Depth**: Minimum read depth observed.
* **Max\_Depth**: Max read depth observed.
* **Avg\_Depth**: Mean read depth observed.
* **Std\_Dev**: Standard deviation of the mean (Avg_Depth).
* **Above\_0.2_Depth**: Percentage of sites which have greater than 0.2 * Avg_Depth.
* **Above\_1_Depth**: Percentage of sites which are above Avg_Depth.
* **Above\_1.8_Depth**: Percentage of sites which have greater than 1.8 * Avg_Depth.
* **Variation\_Coefficient**: The mean of Std_Dev of the mean.

The Average Depth (Avg_Depth) is perhaps the most important field, along with Breadth which will tell you how much of the genome is covered by aligned reads. But the fields such as StdDev and Above_0.2_Depth can give an indication of the variability in the coverage across the genome.

Let’s run weeSAM on our samples:

```
weeSAM --bam 1b.bam --html 1b
```

An explanation of this command is:

1.	**weeSAM**: the name of the program we are using
2.	**--bam**: flag to signify input bam file
3.	**1b.bam**: the name of our bam file to analyse
4.	**--html**: flag to signify output html file
5.	**1b**: the name to label the output html file

If you list the contents of the directory you should see that the folder 1b_html_results has been created:

```
ls
```

Inside this folder is a HTML file that we can view in a web browser (like Firefox or Chrome), the HTML file has the summary statistics and coverage plot so lets take a look and open the html file: 

```
firefox 1b_html_results/1b.html
```

You should see something like this:

*INSERT WEESAM PICTURE HERE*

***
### Questions
**Question 9** – what is the average depth of coverage across the 1b HCV reference genome?
***

Now let’s view the coverage plot by clicking on the hyperlink (blue and underlined) in the Ref_Name column, you should see a coverage plot similar to this:

*INSERT WEESAM COVERAGE PICTURE HERE*

The x-axis represents the genome position, whilst the y-axis represents the Depth of Coverage at each genome position. 

**NB:** The reference sequence filename is 1b_hcv_ref.fasta, but the actual name of the sequence itself is 1|b|EU781827, you can open up the file yourself to check this if you want (head –n1 1b_hcv_ref.fasta).

Although you do expect variation in coverage across the genome, the numerous regions of near zero coverage suggest that the HCV 1b reference is not ideal, and the aligner has struggled to effectively map reads onto it in this regions – presumably because the reference is too divergent from the viral population in the sample at these regions. 

Close FireFox and the weeSAM display before proceeding! 

### Common issue
A common issue here is due to the fact that we have launched firefox from the terminal (wihtout running it background - see advanced linux commands). In order to get our command promot nack (the manager@GCV2023) we need to close firefox window down.

## 5.3: Coverage plot on your own

Your task now is to run weeSAM on the 1a.bam file. So you will need to adapt the previous weeSAM command, a reminder of it is:

```
weeSAM --bam 1b.bam --html 1b
```

***
### Questions
**Question 10** – what is the average depth of coverage across the HCV 1a reference genome?

**Question 11** – how does the coverage plot of HCV 1a compare to HCV 1b? Do you think it is better?
***

## 5.4. Visualisation with Tablet

[Tablet](https://ics.hutton.ac.uk/tablet/) is a tool for the visualisation of next generation sequence assemblies and alignments. It goes beyond simple coverage plots, and allows you to scroll across the genome, zoom into errors of interests, highlight mutations to the reference, and investigate the assembly.

Tablet requires three files:

1.	A bam file, e.g. 1a.bam
2.	A bam index file, e.g. 1a.bam.bai
3.	A reference sequence file: e.g. 1a\_hcv\_ref.fasta

To launch Tablet, type:

```
tablet
```

**NB:** You will not be able to use this command line for other things until you have closed down tablet – but you can open another command line window if you want to leave tablet open and do other things.

You should see the Tablet graphical user interface:

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet.png)

**NB:** Sometimes a small popup window also appears, giving information on how to correctly cite Tablet, with a brief countdown timer.

We want to load in our read alignment from the HCV 1a genome. So **Click** on the **Open Assembly** button on the top menu bar.

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet2.png)

This will launch the Open Assembly window, **Click** **Browse** and then **navigate** to your **~/Richard/HCV/** folder and **Select** the **1a.bam** file for **Primary Assembly**. Afterward, **Click** **Browse** and **select** the **1a\_hcv\_ref.fasta** file for **Reference/Consensus File**, before **Clicking** **Open**.

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet3.png)

After loading you should see the message **-select a contig to begin visualisation-** along with a list of contigs in the left hand panel. In our analysis, we have used a single sequence (the HCV 1a reference sequence), so our contig list only has one entry (the contig NC\_004102.1\_HCV1a), **click** on this entry.

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet4.png)

**NB:** Although our reference file is called 1a\_hcv\_ref.fasta, the actual sequence itself inside the file is called NC\_004102.1\_HCV1a (you can check for yourself if you want to: head –n1 1a\_hcv\_ref.fasta)

**NB:** We only have one contig as our reference sequence only consisted of one sequence (the HCV genome). However, you can align reads to a reference containing multiple sequences, such as the human genome consisting of multiple separate chromosome sequences, or a segmented virus such as influenza consisting of multiple separate segment sequences, or all of the contigs generated from a metagenomics data set. 

Tablet should now load the entire BAM file for visualisation. You can use the **scrollbars** to move across the genome and view all the reads aligned against the reference.

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet5.png)

### Read Display
In the read display, As, Cs, GS and Ts are represented with different colour blocks, Variants are highlighted with Red Text and a different shading, Deletions are represented with Red Asterisks, whilst the location of Insertions is highlighted with red boxes.

**NB:** Like insertions, soft clipping at the end of the reads are also highlighted with red boxes. Soft clipping is where the aligner has decided to discount a portion of the read (at the read’s beginning or end) to improve the alignment.

You can easily jump about the BAM alignment by **clicking** within the **Coverage Overview** window, and the read display will be updated to show this region.

### Variants
One of the (many) useful features of Tablet is the ability to highlight variants. **Slide** the **Variants Slider** all the way to the right hand side to highlight variants. If you now scroll along the genome, you should be able to easily spot consensus level mutations (as virtually every read at a position will have a mutation) and also spot minority variants.

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet6.png)

**NB:** Minority variants could be real viral mutations from the viral population or be errors introduced by RT-PCR or the sequencer itself.

***
### Questions

**Question 12:** Can you find a genome position that has a consensus level mutation?
Hint: hold the mouse over a mutation and the genome location will be reported above the read display in red text
***
 
### Colours Schemes

Tablet also has a few other colour schemes for visualisation, accessed through the “Colour Schemes” tab at the top. Try a few out, perhaps the most commonly used schemes are:

1.	Nucleotide: this is the default one: As, Cs, Gs, and Ts represented with different colours
2.	Direction: reads aligned in the forward direction are highlighted in light blue, whilst those in the reverse direction are highlighted in dark blue
3.	Variants: represents As, Cs, Gs, and Ts with grey, and highlights any mutations with red.

### Exit

Remember that you need to close Tablet down in order to get your command line back.

Either click on the red cross in the top left hand corner, or click the Tablet icon (red circle) (located above Open Assembly) and select Exit Tablet.



 



