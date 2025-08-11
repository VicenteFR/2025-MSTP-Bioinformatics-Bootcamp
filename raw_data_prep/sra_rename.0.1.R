# ---> Rename SRA data on the basis of their SRR IDs.
# Bioproject from paper titled: Drug-tolerant persister cancer cells are vulnerable to GPX4 inhibition
# RunInfo and Accession List retrieved from this page: https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=335583

# ---> Arguments.
# General files.
bio.prj <- 'PRJNA335583'
# Path and file definition.
sra.info.file <- paste0('/tscc/nfs/home/vfajardorosas/raw_data/SraRunInfo.csv')
geo.info.file <- paste0('/tscc/nfs/home/vfajardorosas/raw_data/GEO-SRA_', bio.prj, '.csv')
gen.reports.path <- paste0('/tscc/lustre/ddn/scratch/vfajardorosas/raw_data/', bio.prj)
reports.path <- paste0(gen.reports.path, '/fastqs')

# ---> Load files
sra.info <- read.csv(file=sra.info.file, stringsAsFactors=FALSE)
geo.info <- read.csv(file=geo.info.file, stringsAsFactors=FALSE)

# ---> Main program
# Define new file names.
tmp.data.1 <- list.files(path=reports.path, full.names=TRUE)
tmp.data.1 <- data.frame(
    srr.id=sub(x=basename(tmp.data.1), pattern='\\.fastq$', replacement=''),
    fastq.file=tmp.data.1
)
tmp.cols <- c(`srr.id`='Run', `gsm.id`='SampleName')
tmp.data.2 <- sra.info[, tmp.cols]
colnames(tmp.data.2) <- names(tmp.cols)
tmp.data <- merge(
    x=tmp.data.1, y=tmp.data.2,
    by='srr.id', all=TRUE
)
tmp.data <- merge(
    x=tmp.data, y=geo.info,
    by='gsm.id', all=TRUE
)
tmp.data <- tmp.data[order(tmp.data$srr.id), ]
tmp.data$gsm.occurence <- ave(seq_along(tmp.data$gsm.id), tmp.data$gsm.id, FUN=seq_along)
tmp.data$new.fastq.file <- paste0(
    reports.path, '/',
    tmp.data$final.name,
    '_L00', tmp.data$gsm.occurence,
    '_R1.fastq'
)
# Rename files.
for(tmp.idx in 1:nrow(tmp.data)){
    old.file.name <- tmp.data[tmp.idx, 'fastq.file']
    new.file.name <- tmp.data[tmp.idx, 'new.fastq.file']
    file.rename(from=old.file.name, to=new.file.name)
}

sessionInfo()