# <div align="center"><b>Postprocessing with `samtools`</b></div>

Our STAR alignment script should have generated a bunch of `.sam` (or sequence alignment map) files. These files are "human-readable" alignment files. However, this "human-readable-ness" comes at a huge space/memory cost. To save memory/space we need to generate a compressed (non human-readable) version of these alignment files (you kind conceptualize these as zipped alignment maps). These compressed `.sam` files are known as `.bam` files

Since `.bam` (BAM-BAM!) files are binary, they can only be read by the computer. [`samtools`](https://www.htslib.org/doc/samtools.html) is a great tool that lets us view the contents of bamfiles and perform various manipulations on them. I suggest moseying on over to the samtools documentation and giving it a once over.

Check that you have samtools installed. You can do this with either `which samtools` or `samtools --help`.

Assuming that you do have it installed let's review the basic usage. Simply, commands in samtools follow the below format:


```bash
samtools <command> [options]
```

We will be using 3 samtools commands some of the samtools commands (using the `.sam` files you generated previously).

# 1) `samtools view`

The `view` command can be used for a few different things, but the one we care about now is converting `.sam` files to `.bam` files. The basic usage is as follows:


```bash
samtools view -S -b interesting_file.sam interesting_file.bam
```

# 2) `samtools sort`

We also need to be able to sort a bam file by position and save it to a new file with the extension `.sorted.bam.`:

```bash
samtools sort -@ 8 -o interesting_file.sorted.bam interesting_file.bam
```

# 3) `samtools index`

Lastyl, we need create an index for the sorted bam file. You can think of a index as a table of contents for your `.bam` file. Like a textbook, if we want to read-up on certain topics, we would skim the table of contents to find where we need to jump to.

```bash
samtools index interesting_file.sorted.bam
```

This will create a file called `interesting_file.sorted.bam.bai` which is the index file.

# 4) Scripting `samtools`

Let's put everything together in an automagical script.

`cp ~/fake_script.sh ~/scripts/samtools.sh`

```bash
#!/bin/bash
#SBATCH --job-name=samtools   # Specify a name for the job
#SBATCH --output=samtools.out # Standard output file
#SBATCH --nodes=1                  # Request 1 node
#SBATCH --ntasks=12                 # Request 12 tasks (processes)
#SBATCH --mem=16G                   # Request 16 GB of memory
#SBATCH --time=1:00:00             # Set a time limit of 1 hour
#SBATCH --partition=hotel          # Specify the partition name
#SBATCH --qos=hotel                # Specify the quality of service
#SBATCH --account=htl191           # Specify the account

cd ~/scratch/star_alignment

for x in DMSO_1_ATCACGAligned.out DMSO_2_CGATGTAligned.out DTP_1_CAGATCAligned.out DTP_2_CCGTCCAligned.out DTP_3_GTGAAAAligned.out

do

echo "Beginning $x"

samtools view -S -b $x.sam > $x.bam

samtools sort -@ 8 -o $x.sorted.bam $x.bam

samtools index $x.sorted.bam

done
```

 Notice the -@ 8 flag above for sort? The sorting takes 8 processors, so we need to submit a job requesting at least this amount of resources. Keep in mind that you can include two commands in the same script. Just put one below the other and your second one will run after the first one is finished (yes you can run things in parallel, and if we had a ton of time-intensive files, we certainly would do this, but for the sake of simplicity we are just going to focus on the iterative process for now).

And run with:

```bash
sbatch samtools.sh
```

# <div align="center"><b>Visualizing Reads with IGV</b></div>


Check out the [IGV](http://software.broadinstitute.org/software/igv/) website. Go to Downloads page and follow the instructions based on your operating system. 
 - **Note:** Using the web app version of IGV seems to fail, so we recommend downloading the Windows/Mac version and using that.


In order to view alignments, you need to upload the bam files to an external server/computer (i.e., not TSCC) for viewing. One option is download the bam and the indexed bai files to your desktop and load them from there. But since the files are big and take time to process I have put them to an external server for you to view (without having to download them). There are 10 files in total, two for each condition (2 parental, 3 persister) corresponding to the bam and the index files for each.

**After IGV finishes installing, open IGV:**

Select your genome with **Genomes** - *Load Genome From Server...* Choose hg19.

Upload the bam files from the external server with - **File** - *Load from URL...*<br><br>
The URL links are:<br>
<br>
DMSO-1:<br>
https://mstp-bootcamp-2020.s3-us-west-1.amazonaws.com/DMSO_1_ATCACGAligned.out.sorted.bam<br>
https://mstp-bootcamp-2020.s3-us-west-1.amazonaws.com/DMSO_1_ATCACGAligned.out.sorted.bam.bai<br>
<br>
DMSO-2:<br>
https://mstp-bootcamp-2020.s3-us-west-1.amazonaws.com/DMSO_2_CGATGTAligned.out.sorted.bam<br>
https://mstp-bootcamp-2020.s3-us-west-1.amazonaws.com/DMSO_2_CGATGTAligned.out.sorted.bam.bai<br>
<br>
DTP-1:<br>
https://mstp-bootcamp-2020.s3-us-west-1.amazonaws.com/DTP_1_CAGATCAligned.out.sorted.bam<br>
https://mstp-bootcamp-2020.s3-us-west-1.amazonaws.com/DTP_1_CAGATCAligned.out.sorted.bam.bai<br>
<br>
DTP-2:<br>
https://mstp-bootcamp-2020.s3-us-west-1.amazonaws.com/DTP_2_CCGTCCAligned.out.sorted.bam<br>
https://mstp-bootcamp-2020.s3-us-west-1.amazonaws.com/DTP_2_CCGTCCAligned.out.sorted.bam.bai<br>
<br>
DTP-3:<br>
https://mstp-bootcamp-2020.s3-us-west-1.amazonaws.com/DTP_3_GTGAAAAligned.out.sorted.bam<br>
https://mstp-bootcamp-2020.s3-us-west-1.amazonaws.com/DTP_3_GTGAAAAligned.out.sorted.bam.bai<br>

After you have uploaded the files, IGV will likely ask you to zoom in before visualization. To start type TP53 in the search bar at the top and you should see the text replaced with the actual visualization.

Play around by viewing different genes or chromosome locations. Can you see genes that clearly have fewer reads in the parental vs persister datasets? What about differences in called variants at specific positions? We'll come back to gene analyses later on after running differential expression.

When you are ready to quit IGV, you can save the session with _File - Save session_. Next time you open IGV you can open your saved session without having to reload the BAM files.

# DONE!

---
