# NGS-Pipeline
Next-Generation Sequencing (NGS) pipeline designed for identifying genetic variants through a comparative analysis of the mt1 sequence and the reference genome hq19.

This project was done as part of the UCSD extension course on processing large genomic datasets. https://extendedstudies.ucsd.edu/courses-and-programs/processing-actionable-data-in-genomics

The pipeline encompasses the following key steps:
- Sets up directories and file paths on the local machine.
- Connects to Google Cloud and downloads reference genome and sequencing files.
- Indexes the reference genome using bwa index command.
- Aligns sequencing files to the reference genome using bwa mem.
- Sorts the aligned file using samtools sort and converts .sam to .bam.
- Calls variants using bcftools and outputs a .vcf file.
- Performs a quality check by counting the number of variants found.
- Handles errors at each step and prints corresponding error messages.
- Requires user input for file names to start the pipeline.
- Includes three classes: NGSPipeline, GCloud, and FilePath, handling different tasks.
