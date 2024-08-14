# <div align="center"><b>Return of the owl: How many reads map to a gene?</b></div>

![owl](../static/Day_3/owl.png)

Now that we have our reads aligned to the genome, we can count the number of reads that map to each gene. One way to do this is with the **featureCounts** tool.


Try to find the user documentation for featureCounts. Googling "user manual featurecounts" should help. If not you might find [this](http://manpages.org/featurecounts) helpful. You can also find the documentation in a larger package called subread. [Subread user guide](https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf)

To begin as usual, let's check that you have featureCounts installed properly:

`which featureCounts` --> This should yield `~/miniconda3/bin/featureCounts`

If featureCounts is not installed (though it should be - see [installations](https://github.com/jvtalwar/2022-MSTP-Bioinformatics-Bootcamp/blob/main/Day_0_Setup/Installations/Installations.ipynb) section 4.4), install it now with conda:

`conda install -c bioconda subread`

Refer to the subread manual to find the command you want to run. Scroll down in the manual to section 6: Read summarization. This package also contains an aligner, but as you hopefully recall we used STAR instead. So we will skip that part and only use the part of the package important for quantification of reads mapping to genes.

Note the arguments/flags you can use. From the documentation can you figure out the right combination of flags and inputs to count 'em all? Try discussing the flags/parameters amongst yourselves if you are stumped and see if that can lead to a breakthrough. 

**Hint:** Notice you can count multiple files at the same time! Use this to your advantage to "count" all the sam/bam files at once. Try not to look below until after discussion.

# Creating a SLURM script

Put this output, which is a form of processed data somewhere meaningful (i.e., where YOU know it is)! Some suggestions include making a new folder called `feature_counts` (I'll be doing that), or you can simply include it in another folder with `feature_counts` in the filename. You have complete discretion here, but make sure you understand your setup/organizational system.

Try to write your own SLURM script before looking below!

# DONE!
---

# My SLURM script

Here is my completed script below:

`cp ~/fake_script.sh ~/scripts/featureCounts.sh`

```bash
#!/bin/bash
#SBATCH --job-name=featureCounts   # Specify a name for the job
#SBATCH --output=featureCounts.out # Standard output file
#SBATCH --nodes=1                  # Request 1 node
#SBATCH --ntasks=2                 # Request 2 tasks (processes)
#SBATCH --mem=8G                   # Request 8 GB of memory
#SBATCH --time=1:00:00             # Set a time limit of 1 hour
#SBATCH --partition=hotel          # Specify the partition name
#SBATCH --qos=hotel                # Specify the quality of service
#SBATCH --account=htl191           # Specify the account

featureCounts -a ~/scratch/annotations/hg19/gencode.v19.annotation.gtf -G ~/scratch/annotations/hg19/allchrom.fa -o ~/scratch/feature_counts/hangauer.results.counts ~/scratch/star_alignment/*sorted.bam
```
