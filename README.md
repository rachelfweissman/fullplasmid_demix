### full plasmid sequencing demix
This tool allows you to pool many plasmids into one full plasmid sequencing run, separate out individual plasmid sequences from one or more `.fastq` files, and create a consensus for each plasmid in the mixed pool.

## Desktop App (recommended for most users)

A click-to-run desktop app is available for **macOS and Linux** — no Python or conda setup required for end users. The bioinformatics tools are installed automatically the first time you launch the app.

**[Download the latest release](https://github.com/rachelfweissman/fullplasmid_demix/releases/latest)**

See [GUI_README.md](GUI_README.md) for installation steps, including the one-time `xattr -cr` command needed on macOS to bypass the "damaged" Gatekeeper warning that affects all unsigned apps.

The command-line instructions below are for users who prefer to script the pipeline or run it on a server.

---

## Requirements
- minimap2
- samtools
- flye
- racon
- medaka *(optional, used in full assembly mode)*
- circlator *(optional, used in full assembly mode)*

## Set up Environment
1. Create new conda environment from demixEnv2.yml file

        conda env create -f demixEnv2.yml
2. Activate env

        conda activate demixEnv2

## Set up files
The script requires an Excel file with the following columns:

| Column | Description |
|---|---|
| `Plasmid_Name` | Name of a single plasmid in the pool. Must match the reference fasta filename (without extension). |
| `Unique_Sequences` | One or more sequences unique to this plasmid, comma-separated. Reads containing all listed sequences are assigned to this plasmid. |
| `Colony_ID` | Unique identifier for this sample (e.g. if the same plasmid appears across multiple pools). |
| `Fastq` | Name of the pooled fastq file (without `.fastq` extension) in `fastq_dir`. |

## Command-line Parameters

| Flag | Description |
|---|---|
| `-e`, `--excel_file` | Path to Excel file with plasmid names and unique sequences *(required)* |
| `-f`, `--fastq_dir` | Path to directory with pooled fastq files *(required)* |
| `-o`, `--output_dir` | Path to output directory for results *(required)* |
| `-r`, `--ref_dir` | Path to directory with reference plasmid fasta files *(optional)* |
| `-q`, `--quick_method` | Use quick alignment method (minimap2 + samtools consensus) instead of full de novo assembly. Requires a reference directory. |
| `-t`, `--threads` | Max threads per external tool (default: 4). Parallel jobs are scaled automatically based on available CPUs. |
| `-k`, `--keep_temp` | Keep intermediate files (alignments, assemblies, Racon/Medaka outputs). |

## Output
Each plasmid gets its own subfolder `{output_dir}/{Plasmid_Name}_{Colony_ID}/` containing:

- `*_reads.fasta` — reads assigned to this plasmid
- `*_sorted.bam` / `*.bam.bai` — reads aligned to the consensus
- `*_consensus.fa` — final consensus sequence
- `*_coverage.html` — interactive coverage plot

A shared `{output_dir}/log/` subfolder contains:
- `fullPlasmidSeq_demix.log` — full run log
- `*_readLengths.html` — read length histogram for each plasmid

Unassigned reads are written to `{output_dir}/unused_reads/`.

## How It Works

Processing runs in two phases:

**Phase 1 — Demixing:** Each FASTQ file is loaded and scanned in a single pass. Each read is tested against the plasmids in order and assigned to the first one whose unique sequences all appear in the read. Different FASTQ files are processed in parallel.

**Phase 2 — Consensus:** All per-plasmid alignment and consensus jobs run in parallel, with the number of concurrent jobs scaled to keep total CPU usage bounded (`cpu_count / threads_per_job`).

### Assembly modes

**Quick mode** (`-q`): Maps reads directly to the provided reference using minimap2, then calls consensus with `samtools consensus`. Fast, but requires a reference fasta.

**Full assembly mode** (default): Runs de novo assembly with Flye, polishes with two rounds of Racon, and optionally polishes further with Medaka. Does not require a reference.

## Usage Examples

Basic usage (full assembly):
```
python fullPlasmidSeq_demix_RFW.py -e plasmid_uniqueseq.xlsx -f raw_reads -o output_dir
```

Quick mode with reference:
```
python fullPlasmidSeq_demix_RFW.py -e plasmid_uniqueseq.xlsx -r reference_plasmids -f raw_reads -o output_dir -q
```

Increase threads and keep intermediate files:
```
python fullPlasmidSeq_demix_RFW.py -e plasmid_uniqueseq.xlsx -f raw_reads -o output_dir -t 8 -k
```