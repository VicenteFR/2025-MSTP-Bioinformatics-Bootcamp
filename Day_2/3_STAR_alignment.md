# <div align="center"><b><span style="color:yellow">STAR</span>: Aligning Reads to the Genome (Here We Go!)</b></div>

![mario_star](../static/Day_2/mario_star.png)

# 1) Generate <span style="color:yellow">STAR</span> genome index

Open the STAR user [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). We will go through this briefly together to get an understanding of how to read documentation.

## 1a) Download human fasta file

Open UCSC genome [browser](https://genome.ucsc.edu/). The link to the specific annotations we will use is provided below, but first take a look through the website to see all the available annotations and features. We will briefly go through this together. We will use UCSC to download the **chromosome fasta files** that are needed to build the STAR index. Use the wget command followed by a copy of the web link address to download the files to TSCC. The annotations are located [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/). Scroll to the bottom of the page and get the link for **chromFa.tar.gz**. Once you have made a folder on your TSCC account, move into it so your annotations will land in the proper place during downloading.

```bash
cd ~/scratch
mkdir annotations
cd annotations
mkdir hg19
cd ~/scratch/annotations/hg19/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
```

This will download a zipped file that you will proceed to unzip with:

```bash
tar -xvf chromFa.tar.gz
```

Unfortunately this downloads everything and we only want the chr#.fa files. Let's check how many files we actually everye:

```bash
ls | wc -l
```
```plaintext
94
```

We will remove the unneeded files in this directory with `rm`. To remove more than one file at a time you have to use the `-r` flag (recursive). Remember you can use the star character to remove all things that contain common characters. For example:

```bash
rm -r *random*
rm -r *chrUn*
rm -r *hap*
```

Once the folder is clean and contains only one fasta file per chromosome (and the original tar.gz file), you can merge them all together using `cat` and assign the output to a new file called allchrom.fa using `>`. This is the chromosome fasta file that you will need to use to generate the genome index.

```bash
cat *.fa > allchrom.fa
```

*NOTE - the > character saves the result of your command to a new file.*

## 1b) Download gtf annotation from gencode

We will use gencode release (19) for genome build GRCh37 (hg19). We want the gtf file of the comprehensive gencode annotation for chromosomes. 

```bash
cd ~/scratch/annotations/hg19/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
```

Unzip the compressed "wgot" file with `gunzip`:

```bash
gunzip gencode.v19.annotation.gtf.gz
```

*Note the different compression formats. This file was gzipped so it needs to unzipped with gunzip. The fasta files were tar compressed and required tar -xvf to unzip*.

## 1c) Generating Genome Indices:

Before we begin scriptifying things we need to make a **<span style="color:red">genomeDir</span>** folder so STAR doesn't error out. Why do you think this is (The documentation, specifically section 2.1 under --genomeDir may be helpful)? 

We can easily construct this from our current path with the following:

```bash
mkdir star_own
```

STAR requires a lot of processing power, we are going to submit this command as a job to the cluster. Remember that handy fake submission script we made (at the end of Introduction_to_bash_scripting.ipynb)? Let's use it here by copying and updating the necessary parameters:

`cd ~/scripts` <-- If this causes an error it likely means either that you did not make a scripts folder from your home directory, or you made your scripts folder somewhere else. To make a scripts folder from your home directory do the following before rerunning the above command: 

```bash
cd
mkdir scripts
```

Now copy our fake script into that directory with a new, meaningful name such as star_generate_index.sh

```bash
cp ~/fake_script.sh ~/scripts/star_generate_index.sh
```

For this script, we will use a time limit of 5 hours, 1 node, and 16 processors.

Use the STAR manual to decide how your STAR command should look like. Once you have decided what your STAR command should look like, add it to your **star_generate_index.sh** script below the PBS flags.

In case you are lost here is what I did:

```bash
mkdir ~/scratch/annotations/hg19/star
```

```bash
#!/bin/bash
#SBATCH --job-name=star_generate_index   # Specify a name for the job
#SBATCH --output=star_generate_index.out # Standard output file
#SBATCH --error=star_generate_index.err # Standard error file
#SBATCH --nodes=1                  # Request 1 node
#SBATCH --ntasks=16                # Request 16 tasks (processes)
#SBATCH --mem=32G                   # Request 32 GB of memory
#SBATCH --time=5:00:00             # Set a time limit of 5 hours
#SBATCH --partition=hotel          # Specify the partition name
#SBATCH --qos=hotel                # Specify the quality of service
#SBATCH --account=htl191           # Specify the account

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ~/scratch/annotations/hg19/star_own --genomeFastaFiles ~/scratch/annotations/hg19/allchrom.fa --sjdbGTFfile ~/scratch/annotations/hg19/gencode.v19.annotation.gtf --sjdbOverhang 49 --outFileNamePrefix ~/scratch/annotations/hg19/star_own
```

#### Submit your script to the cluster:

First make sure the conda environment with STAR in it is activated:

```bash
conda activate 2025-mstp-bootcamp
sbatch star_generate_index.sh
```

This will take a **long** while to complete (**ETA**: 5 hrs). If successful your out file, *star_generate_index.out*, should look something like this: <br><br>

```bash
	/tscc/nfs/home/vfajardorosas/opt/miniconda3/envs/2025-mstp-bootcamp/bin/STAR-avx2 --runThreadN 16 --runMode genomeGenerate --genomeDir /tscc/nfs/home/vfajardorosas/scratch/annotations/hg19/star --genomeFastaFiles /tscc/nfs/home/vfajardorosas/scratch/annotations/hg19/allchrom.fa --sjdbGTFfile /tscc/nfs/home/vfajardorosas/scratch/annotations/hg19/gencode.v19.annotation.gtf --sjdbOverhang 49 --outFileNamePrefix /tscc/nfs/home/vfajardorosas/scratch/annotations/hg19/star
	STAR version: 2.7.11b   compiled: 2024-11-25T09:14:51+0000 :/opt/conda/conda-bld/star_1732525954305/work/source
Aug 14 18:36:20 ..... started STAR run
Aug 14 18:36:20 ... starting to generate Genome files
Aug 14 18:37:00 ..... processing annotations GTF
Aug 14 18:37:22 ... starting to sort Suffix Array. This may take a long time...
Aug 14 18:37:33 ... sorting Suffix Array chunks and saving them to disk...
Aug 14 20:13:57 ... loading chunks from disk, packing SA...
Aug 14 20:16:13 ... finished generating suffix array
Aug 14 20:16:13 ... generating Suffix Array index
Aug 14 20:19:00 ... completed Suffix Array index
Aug 14 20:19:00 ..... inserting junctions into the genome indices
Aug 14 20:24:02 ... writing Genome to disk ...
Aug 14 20:24:06 ... writing Suffix Array to disk ...
Aug 14 20:24:44 ... writing SAindex to disk
Aug 14 20:24:47 ..... finished successfully
```

## 1d) Check the status of your job

```bash
squeue -u [your_username]
```

Take a look at the status (The column labeled S). `PD` means your job is pending to run and has not started yet. `R` means your job is running (you will see the time updated according to how long it has been running). `C` means your job is complete.

## 1e) Make a copy of an already existing STAR genome index

Since the genome index genration takes a long while, we can use a version of it that has been previously generated by making a copy to our annotations directory. **Note that this is possible because the STAR genome index is deterministic and not stochastic**, which means that `genomeGenerate` will always provide the exact same output given the same input.<br/>
To do this, let's type the following commands:

```bash
cd ~/scratch/annotations/hg19
cp -R /tscc/nfs/home/hkcarter/star .
```

## 1f) Dealing with Errors

The error file you have designated will output any errors or warnings associated with the job.

If you have an error with your job, see if you can understand the error and try go back into your script to correct it/them. Common errors are misspelling or incorrect paths to data.

*Example error*

`EXITING: FATAL INPUT ERROR: unrecoginzed parameter name "sjdbGTFFile" in input "Command-Line-Initial"
SOLUTION: use correct parameter name (check the manual)`

```bash
Jul 21 14:19:02 ...... FATAL ERROR, exiting
```

In this case, you would go back and change the typo. The second **F** should not be capitalized.

# 2) Map reads to the genome

Again let's make a folder for the output of our alignments:

```bash
mkdir ~/scratch/star_alignment
```

Once your genome index job is complete, you can move onto the next step of mapping your reads to the genome. Again you can copy your `fake_script.sh` file and make the necessary changes for this particular job submission.

```bash
cp ~/fake_script.sh ~/scripts/star_align.sh
```

Using the STAR manual, try to write out the command for mapping reads. If you don't like a challenge though I'll be magnanimous and elucidate the esoteric thus allowing you to embrace your inner-astronomer and go STAR-gazing:

```bash
#!/bin/bash
#SBATCH --job-name=star_align   # Specify a name for the job
#SBATCH --output=star_align.out # Standard output file
#SBATCH --error=star_align.err # Standard error file
#SBATCH --nodes=1                  # Request 1 node
#SBATCH --ntasks=16                 # Request 16 tasks (processes)
#SBATCH --mem=32G                   # Request 16 GB of memory
#SBATCH --time=3:00:00             # Set a time limit of 3 hours
#SBATCH --partition=hotel          # Specify the partition name
#SBATCH --qos=hotel                # Specify the quality of service
#SBATCH --account=htl191           # Specify the account

STAR --runThreadN 16 --genomeDir ~/scratch/annotations/hg19/star --readFilesIn ~/scratch/fastq/DMSO_1_ATCACG.combined.fastq  --outFileNamePrefix ~/scratch/star_alignment/DMSO_1_ATCACG

STAR --runThreadN 16 --genomeDir ~/scratch/annotations/hg19/star --readFilesIn ~/scratch/fastq/DMSO_2_CGATGT.combined.fastq --outFileNamePrefix ~/scratch/star_alignment/DMSO_2_CGATGT

STAR --runThreadN 16 --genomeDir ~/scratch/annotations/hg19/star --readFilesIn ~/scratch/fastq/DTP_1_CAGATC.combined.fastq --outFileNamePrefix ~/scratch/star_alignment/DTP_1_CAGATC

STAR --runThreadN 16 --genomeDir ~/scratch/annotations/hg19/star --readFilesIn ~/scratch/fastq/DTP_2_CCGTCC.combined.fastq --outFileNamePrefix ~/scratch/star_alignment/DTP_2_CCGTCC

STAR --runThreadN 16 --genomeDir ~/scratch/annotations/hg19/star --readFilesIn ~/scratch/fastq/DTP_3_GTGAAA.combined.fastq --outFileNamePrefix ~/scratch/star_alignment/DTP_3_GTGAAA
```

# DONE!
When finished, it's time to do the script submission thing. So go ahead and unleash the `sbatch star_align.sh`.<br>
Congratulations! You are now a full-fledged bio-stronomer (or rock-STAR if you prefer)!🌟

 ---
