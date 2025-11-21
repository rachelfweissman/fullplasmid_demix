### fullplasmid_demix 
This tool allows you to pool many plasmids into one full plasmid sequencing run, separate out individual plasmid sequences from one .fastq, and create a consensus for each plasmid in the mixed pool.

This script demixes pooled plasmid sequencing files by searching for unique sequences for each plasmid in the pooled .fastq files, assembles each plasmid, and generates consensus sequences.

## Requirements
<ul>
    <li>minimap2</li>
    <li>samtools</li>
    <li>flye</li>
    <li>racon</li>
    <li>medaka (optional)</li>
    <li>circlator (optional)</li>
</ul>

## Set up Environment
1. Create new conda environment from demixEnv2.yml file

        conda env create -f demixEnv2.yml
2. Activate env

        conda activate demixEnv2

## Set up files
The script requires an excel file with the following columns: 
<ul>
    <li>Plasmid_name - name of single plasmid in pool, MUST match reference fasta file name</li>
    <li>Unique_Sequences - unique sequences for each plasmid in the pool, separated by commas</li>
    <li>Colony_ID - unique identifier if sequencing the same plasmid across multiple pools</li>
    <li>Fastq - name of fastq file of pooled reads in fastq_dir</li>
</ul>
    
## Command-line Parameters
<ul>
    <li>-e, --excel_file: path to excel file with plasmid names and unique sequences</li>
    <li>-r, --ref_dir: (optional) path to directory with reference plasmid fasta files</li>
    <li>-f, --fastq_dir: path to directory with pooled fastq files</li>
    <li>-o, --output_dir: path to output directory for results</li>
    <li>-k, --keep_temp: keep temporary files (intermediate alignments, assemblies)</li>
    <li>-t, --threads: maximum number of threads to use for external tools (default: 4)</li>
    <li>-q, --quick_method: Use quick alignment method (minimap2+samtools consensus) instead of full assembly with racon. Default is false</li>
</ul>

## Output
For each plasmid in the excel file, the script will create a subfolder in the output directory with the following files:
<ul>
    <li>{Plasmid_name}_{Colony_ID}_reads.fasta - fasta file with reads containing unique sequences for the plasmid</li>
    <li>{Plasmid_name}_{Colony_ID}_sorted.bam - sorted bam file of aligned reads to reference plasmid</li>
    <li>{Plasmid_name}_{Colony_ID}_consensus.fa - consensus sequence of aligned reads</li>
    <li>{Plasmid_name}_{Colony_ID}_coverage.html - coverage plot of aligned reads</li>
    <li>{Plasmid_name}_{Colony_ID}_readLengths.html - histogram of read lengths</li>
</ul>

## Usage Examples

### Run with Python Script
Basic usage:
```
python fullPlasmidSeq_demix_RFW.py -e example/plasmid_uniqueseq.xlsx -f example/raw_reads -o example/output_directory
```

With reference directory:
```
python fullPlasmidSeq_demix_RFW.py -e example/plasmid_uniqueseq.xlsx -r reference_plasmids -f raw_reads -o output_directory
```

Keep temporary files:
```
python fullPlasmidSeq_demix_RFW.py -e example/plasmid_uniqueseq.xlsx -f raw_reads -o output_directory -k
```

Specify thread count for external tools:
```
python fullPlasmidSeq_demix_RFW.py -e example/plasmid_uniqueseq.xlsx -f raw_reads -o output_directory -t 8
```

## Performance Optimizations

The script includes several optimizations to improve processing speed:

### Resource Management
- Limits thread usage for each external bioinformatics tool (minimap2, flye, racon, etc.)
- Controls memory usage by setting environment variables for numerical libraries

### Efficient Processing
- Groups plasmids by FASTQ file for minimal file loading
- Uses optimized parameters for Flye assembly (--meta option)
- Implements early returns in sequence searches
- Sets timeouts for external commands to prevent indefinite hangs

### I/O Optimization
- Indexes each FASTQ file only once
- Processes reads sequentially by group
- Cleans up intermediate files to save disk space
