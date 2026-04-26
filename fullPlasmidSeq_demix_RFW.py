import os
import subprocess
import tempfile
import warnings
import logging
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import plotly.graph_objects as go
import argparse
import shutil
import time
import glob
import concurrent.futures

"""
This script demixes pooled plasmid sequencing files by searching for unique sequences for each plasmid in the pooled .fastq files. 

The script requires an excel file with the following columns: 
    Plasmid_name - name of single plasmid in pool, MUST match reference fasta file name
    Unique_Sequences - unique sequences for each plasmid in the pool, separated by commas
    Colony_ID - unique identifier if sequencing the same plasmid across multiple pools
    Fastq - name of fastq file of pooled reads in fastq_dir
    
Parameters:
    input_excel: path to excel file with plasmid names and unique sequences for each plasmid as described above
    ref_dir: path to directory with reference plasmid fasta files
    fastq_dir: path to directory with pooled fastq files
    output_dir: path to output directory for results
    
Output:
    For each plasmid in the excel file, the script will create a subfolder in the output directory with the following files:
        {Plasmid_name}_{Colony_ID}_reads.fasta - fasta file with reads containing unique sequences for the plasmid
        {Plasmid_name}_{Colony_ID}_sorted.bam - sorted bam file of aligned reads to reference plasmid
        {Plasmid_name}_{Colony_ID}_consensus.fa - consensus sequence of aligned reads
        {Plasmid_name}_{Colony_ID}_coverage.html - coverage plot of aligned reads
        {Plasmid_name}_{Colony_ID}_readLengths.html - histogram of read lengths
"""

# Configure logging (file handler is added in main() once output_dir is known)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()]
)

# Set thread limits for external tools to optimize resource usage
# Most bioinformatics tools use multiple threads by default which can cause resource contention
MAX_THREADS = 4  # Adjust based on your system's capabilities
os.environ["OMP_NUM_THREADS"] = str(MAX_THREADS)
os.environ["OPENBLAS_NUM_THREADS"] = str(MAX_THREADS)

def read_sequence_sets_from_excel(excel_file):
    logging.info(f"Reading sequence sets from Excel file: {excel_file}")
    required_columns = {'Plasmid_Name', 'Unique_Sequences', 'Colony_ID', 'Fastq'}
    df = pd.read_excel(excel_file)
    df['Unique_Sequences'] = df['Unique_Sequences'].str.upper()
    if not required_columns.issubset(df.columns):
        missing = required_columns - set(df.columns)
        logging.error(f"Excel file is missing required columns: {', '.join(missing)}")
        raise ValueError(f"Excel file is missing required columns: {', '.join(missing)}")
    return df

def read_fasta(file_path):
    logging.info(f"Reading FASTA file: {file_path}")
    return SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))

def write_fasta(sequences, file_path):
    logging.info(f"Writing FASTA file: {file_path}")
    with open(file_path, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def write_fastq(sequences, file_path):
    logging.info(f"Writing FASTQ file: {file_path}")
    with open(file_path, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fastq")

# Optimized sequence search that returns early when match is found
def find_sequences(index_file, sequence_sets):
    logging.info("Finding sequences matching unique identifiers.")
    matched_reads = []
    unmatched_reads = list(index_file.values())
    
    # If no sequences to find, return early
    if not sequence_sets or not unmatched_reads:
        return [], unmatched_reads
        
    for seq_set in sequence_sets:
        sequences = [seq.strip() for seq in seq_set.split(',') if seq.strip()]
        temp_matched = []
        temp_unmatched = []
        
        # If no sequences in this set, skip
        if not sequences:
            continue
            
        for record in unmatched_reads:
            # Convert to string once for multiple checks
            record_seq = str(record.seq)
            if all(seq in record_seq for seq in sequences):
                temp_matched.append(record)
            else:
                temp_unmatched.append(record)
        
        matched_reads.extend(temp_matched)
        unmatched_reads = temp_unmatched
        
    logging.info(f"Found {len(matched_reads)} matched reads.")
    return matched_reads, unmatched_reads

def run_command(command):
    logging.info(f"Running command: {command}")
    try:
        # Use a timeout to prevent commands from hanging indefinitely
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE, text=True, timeout=7200)  # 2-hour timeout
        if result.stderr:
            warnings.warn(f"Command '{command}' had the following warnings/errors:\n{result.stderr}")
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"Command '{command}' failed with error:\n{e.stderr}")
        raise
    except subprocess.TimeoutExpired:
        logging.error(f"Command '{command}' timed out after 7200 seconds")
        raise
def align_sequences_quick(input_fasta, reference, output_dir, outputname, log_dir=None):
    if not os.path.exists(reference):
        logging.warning(f"Reference FASTA file {reference} not found. Skipping alignment for {outputname}.")
        return

    makeSam_cmd = f'minimap2 -t {MAX_THREADS} -ax map-ont {reference} {input_fasta} 2>/dev/null| samtools view -b -F 0x900 | samtools sort -o {os.path.join(output_dir, outputname)}_sorted.bam'
    index_cmd = f'samtools index {os.path.join(output_dir, outputname)}_sorted.bam'
    samtools_consensus_cmd = f"samtools consensus --config r10.4_sup --output {os.path.join(output_dir, outputname)}_consensus.fa {os.path.join(output_dir, outputname)}_sorted.bam"
    run_command(makeSam_cmd)
    run_command(index_cmd)
    run_command(samtools_consensus_cmd)
    make_plots(os.path.join(output_dir, f"{outputname}_sorted.bam"), output_dir, outputname, log_dir=log_dir)

def align_sequences(input_fasta, output_dir, outputname, keep_temp=False, log_dir=None):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: De novo assembly with Flye (circular-aware for plasmids)
    flye_outdir = os.path.join(output_dir, f"{outputname}_flye")
    
    # Optimized parameters for Flye to improve speed
    # --meta speeds up processing but may reduce assembly quality slightly
    # --threads limits CPU utilization to prevent resource contention
    flye_cmd = f"flye --nano-raw {input_fasta} --out-dir {flye_outdir} --plasmids --meta --threads {MAX_THREADS}"
    run_command(flye_cmd)
    
    # Get the assembly file from Flye
    assembly_file = os.path.join(flye_outdir, "assembly.fasta")
    
    # Step 2: Polish with Racon (2 rounds)
    # First round - use optimized parameters for minimap2
    temp_align_file = os.path.join(output_dir, f"{outputname}_align.paf")
    minimap_cmd = f"minimap2 -t {MAX_THREADS} -x ava-ont {assembly_file} {input_fasta} > {temp_align_file}"
    run_command(minimap_cmd)
    
    racon_output1 = os.path.join(output_dir, f"{outputname}_racon1.fasta")
    # Limit Racon threads
    racon_cmd1 = f"racon -t {MAX_THREADS} {input_fasta} {temp_align_file} {assembly_file} > {racon_output1}"
    run_command(racon_cmd1)
    
    # Second round of Racon
    temp_align_file2 = os.path.join(output_dir, f"{outputname}_align2.paf")
    minimap_cmd2 = f"minimap2 -t {MAX_THREADS} -x ava-ont {racon_output1} {input_fasta} > {temp_align_file2}"
    run_command(minimap_cmd2)
    
    racon_output2 = os.path.join(output_dir, f"{outputname}_racon2.fasta")
    racon_cmd2 = f"racon -t {MAX_THREADS} {input_fasta} {temp_align_file2} {racon_output1} > {racon_output2}"
    run_command(racon_cmd2)
    
    # Step 3: Final polishing with Medaka
    medaka_outdir = os.path.join(output_dir, f"{outputname}_medaka")
    consensus_output = os.path.join(output_dir, f"{outputname}_consensus.fa")
    
    # Limit threads for Medaka
    medaka_cmd = f"medaka_consensus -i {input_fasta} -d {racon_output2} -o {medaka_outdir} -m r941_min_hac_g507 -t {MAX_THREADS}"
    
    try:
        run_command(medaka_cmd)
        # Copy final consensus to expected location with expected name
        final_consensus = os.path.join(medaka_outdir, "consensus.fasta")
        shutil.copy(final_consensus, consensus_output)
        # Delete medaka directory after copying the consensus file if not keeping temp files
        if not keep_temp and os.path.exists(medaka_outdir) and os.path.isdir(medaka_outdir):
            shutil.rmtree(medaka_outdir)
    except Exception as e:
        logging.warning(f"Medaka failed: {e}. Using Racon output as consensus.")
        # Use Racon output as consensus
        shutil.copy(racon_output2, consensus_output)
    
    # Optional: Use Circlator to fix start point of circular sequence
    circlator_cmd = f"circlator fixstart {consensus_output} {os.path.join(output_dir, f'{outputname}_circular.fa')}"
    try:
        run_command(circlator_cmd)
    except Exception as e:
        print(f"Warning: Circlator failed, but continuing. Error: {e}")
    
    # Generate a sorted BAM file for plotting (using the final consensus as reference)
    # Limit threads for minimap2 and samtools
    makeSam_cmd = f'minimap2 -t {MAX_THREADS} -ax map-ont {consensus_output} {input_fasta} 2>/dev/null | samtools view -@ {MAX_THREADS} -b -F 0x900 | samtools sort -@ {MAX_THREADS} -o {os.path.join(output_dir, outputname)}_sorted.bam'
    index_cmd = f'samtools index {os.path.join(output_dir, outputname)}_sorted.bam'
    
    run_command(makeSam_cmd)
    run_command(index_cmd)
    
    # Generate plots using the alignment to the final consensus
    make_plots(os.path.join(output_dir, f"{outputname}_sorted.bam"), output_dir, outputname, log_dir=log_dir)
    
    # Clean up intermediate files if not keeping temp files
    if not keep_temp:
        logging.info(f"Cleaning up intermediate files for {outputname}")
        # Delete alignment files
        if os.path.exists(temp_align_file):
            os.remove(temp_align_file)
        if os.path.exists(temp_align_file2):
            os.remove(temp_align_file2)
        
        # Delete intermediate Racon file
        if os.path.exists(racon_output1) and racon_output1 != racon_output2:
            os.remove(racon_output1)
        
        # Delete flye directory
        if os.path.exists(flye_outdir) and os.path.isdir(flye_outdir):
            shutil.rmtree(flye_outdir)
    
    return consensus_output

def make_plots(bamfile, output_dir, outputname, log_dir=None):
    logging.info(f"Generating plots for {outputname}.")
    lengths_dir = log_dir if log_dir else output_dir
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        coverage_file = temp_file.name
    with tempfile.NamedTemporaryFile(delete=False, suffix='.txt') as tmp_len:
        lengths_file = tmp_len.name
    try:
        coverage_cmd = f'samtools depth -@ {MAX_THREADS} -a {bamfile} > {coverage_file}'
        run_command(coverage_cmd)
        coverage_data = pd.read_csv(coverage_file, sep='\t', header=None, names=['name', 'position', 'coverage'])
        fig = go.Figure()
        for plasmidname in coverage_data['name'].unique():
            chrom_data = coverage_data[coverage_data['name'] == plasmidname]
            fig.add_trace(go.Scatter(
                x=chrom_data['position'],
                y=chrom_data['coverage'],
                mode='lines',
                name=plasmidname
            ))
        fig.update_layout(
            title='Coverage Plot',
            xaxis_title='Position',
            yaxis_title='Coverage Depth',
            template='plotly_white',
            showlegend=True,
            hovermode='x unified',
        )
        fig.write_html(f"{os.path.join(output_dir, outputname)}_coverage.html")

        lengths_cmd = f"samtools view -@ {MAX_THREADS} {bamfile} | awk '{{print length($10)}}' > {lengths_file}"
        run_command(lengths_cmd)
        read_lengths = pd.read_csv(lengths_file, header=None, names=['length'])
        fig = go.Figure()
        fig.add_trace(go.Histogram(
            x=read_lengths['length'],
            nbinsx=50,
            marker_color='blue'
        ))
        fig.update_layout(
            title='Read Length Distribution',
            xaxis_title='Read Length (bp)',
            yaxis_title='Count',
            template='plotly_white',
            showlegend=False,
        )
        fig.write_html(f"{os.path.join(lengths_dir, outputname)}_readLengths.html")
    finally:
        os.remove(coverage_file)
        if os.path.exists(lengths_file):
            os.remove(lengths_file)

def _demix_fastq(fastq_name, group, fastq_dir, output_dir, ref_dir):
    """Load one FASTQ file and assign reads to plasmids in a single pass.

    Pre-computes str(record.seq) once per read so substring searches never
    repeat that conversion. Uses SeqIO.parse (faster than SeqIO.index for
    bulk loading) and touches each read exactly once regardless of how many
    plasmids are in the pool.

    Returns (tasks, fastq_name, unmatched_records).
    """
    fastq_path = os.path.join(fastq_dir, f"{fastq_name}.fastq")
    logging.info(f"Loading FASTQ: {fastq_path}")

    # Parse all reads; pre-compute string sequence once per record
    reads = [(r, str(r.seq)) for r in SeqIO.parse(fastq_path, "fastq")]
    logging.info(f"Loaded {len(reads)} reads from {fastq_name}")

    # Build ordered list of (plasmid_name, colony_id, required_sequences)
    plasmid_order = [
        (row['Plasmid_Name'], row['Colony_ID'],
         [s.strip() for s in str(row['Unique_Sequences']).split(',') if s.strip()])
        for _, row in group.iterrows()
    ]

    # Single pass: assign each read to the first plasmid whose sequences all match.
    # Using a for/else so we only append to unmatched when no plasmid claimed the read.
    buckets = {(pname, cid): [] for pname, cid, _ in plasmid_order}
    unmatched = []
    for record, record_seq in reads:
        for pname, cid, seqs in plasmid_order:
            if all(s in record_seq for s in seqs):
                buckets[(pname, cid)].append(record)
                break
        else:
            unmatched.append(record)

    # Write per-plasmid FASTA files and build consensus task descriptors
    tasks = []
    for pname, cid, _ in plasmid_order:
        sequences = buckets[(pname, cid)]
        if sequences:
            group_output_dir = os.path.join(output_dir, f"{pname}_{cid}")
            os.makedirs(group_output_dir, exist_ok=True)
            group_fasta = os.path.join(group_output_dir, f"{pname}_{cid}_reads.fasta")
            reference_plasmid = os.path.join(ref_dir, f"{pname}.fa") if ref_dir else None
            write_fasta(sequences, group_fasta)
            tasks.append({
                'group_fasta': group_fasta,
                'reference_plasmid': reference_plasmid,
                'group_output_dir': group_output_dir,
                'outputname': f"{pname}_{cid}",
                'plasmid_name': pname,
                'colony_id': cid,
            })
        else:
            logging.warning(f"No sequences found for {pname}, Colony {cid}")

    logging.info(f"Demixed {fastq_name}: {len(tasks)} plasmids matched, {len(unmatched)} unused reads")
    return tasks, fastq_name, unmatched


def _run_consensus(task, quick_method, keep_temp, log_dir=None):
    """Run alignment and consensus for a single plasmid task (called in parallel)."""
    group_fasta = task['group_fasta']
    reference_plasmid = task['reference_plasmid']
    group_output_dir = task['group_output_dir']
    outputname = task['outputname']
    plasmid_name = task['plasmid_name']
    colony_id = task['colony_id']
    try:
        if reference_plasmid and quick_method:
            logging.info(f"Quick mapping to reference for {outputname}")
            align_sequences_quick(group_fasta, reference_plasmid, group_output_dir, outputname, log_dir=log_dir)
        else:
            if reference_plasmid:
                logging.info(f"Reference available but quick_method not set; running de novo assembly for {outputname}")
            else:
                logging.info(f"No reference provided; running de novo assembly for {outputname}")
            align_sequences(group_fasta, group_output_dir, outputname, keep_temp=keep_temp, log_dir=log_dir)
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to process {plasmid_name}, Colony {colony_id}. Error: {e}")


def main(excel_file, ref_dir, fastq_dir, output_dir, keep_temp=False, threads=4, quick_method=False):
    """
    Main processing function.

    Parameters:
      - excel_file: path to excel file
      - ref_dir: path to reference fasta directory (can be None or '')
      - fastq_dir: path to fastq directory
      - output_dir: path to write outputs
      - keep_temp: bool, whether to keep temporary files
      - threads: int, max threads for external tools
      - quick_method: bool, use quick alignment/consensus pipeline instead of full assembly
    """
    global MAX_THREADS
    MAX_THREADS = threads
    os.environ["OMP_NUM_THREADS"] = str(MAX_THREADS)
    os.environ["OPENBLAS_NUM_THREADS"] = str(MAX_THREADS)

    log_dir = os.path.join(output_dir, "log")
    os.makedirs(log_dir, exist_ok=True)
    log_handler = logging.FileHandler(os.path.join(log_dir, "fullPlasmidSeq_demix.log"), mode="w")
    log_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logging.getLogger().addHandler(log_handler)

    logging.info(f"Options: keep_temp={keep_temp}, threads={threads}, quick_method={quick_method}")
    start_time = time.time()
    logging.info("Starting main process.")

    df = read_sequence_sets_from_excel(excel_file)
    unused_reads_dir = os.path.join(output_dir, "unused_reads")
    os.makedirs(unused_reads_dir, exist_ok=True)

    # --- Phase 1: Demixing (parallel across FASTQ files) ---
    # Each FASTQ file is independent, so load and demix them concurrently.
    # Within each file, _demix_fastq does a single pass through all reads,
    # assigning each read to the first plasmid whose unique sequences match.
    logging.info("Phase 1: Demixing — finding unique reads per plasmid.")

    fastq_groups = list(df.groupby('Fastq'))
    n_fastq_workers = min(len(fastq_groups), max(1, os.cpu_count() or 1))

    plasmid_tasks = []
    remaining_reads = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_fastq_workers) as executor:
        futures = {
            executor.submit(_demix_fastq, fastq_name, group, fastq_dir, output_dir, ref_dir): fastq_name
            for fastq_name, group in fastq_groups
        }
        for future in concurrent.futures.as_completed(futures):
            fastq_name = futures[future]
            try:
                tasks, fname, unused = future.result()
                plasmid_tasks.extend(tasks)
                remaining_reads[fname] = unused
            except Exception as e:
                logging.error(f"Failed to demix {fastq_name}: {e}")
                raise

    # Write unused reads after all demixing is complete
    for fastq_name, unused_reads in remaining_reads.items():
        if unused_reads:
            unused_reads_fastq = os.path.join(unused_reads_dir, f"{fastq_name}_unused.fastq")
            write_fastq(unused_reads, unused_reads_fastq)

    logging.info(f"Phase 1 complete. {len(plasmid_tasks)} plasmids to process.")

    # --- Phase 2: Consensus in parallel ---
    # Each plasmid's alignment/consensus is independent; run them concurrently.
    # Parallel workers = CPUs / per-tool threads, so total thread usage stays bounded.
    n_workers = max(1, (os.cpu_count() or 1) // MAX_THREADS)
    logging.info(f"Phase 2: Running consensus for {len(plasmid_tasks)} plasmids ({n_workers} parallel workers).")

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(_run_consensus, task, quick_method, keep_temp, log_dir): task['outputname']
            for task in plasmid_tasks
        }
        for future in concurrent.futures.as_completed(futures):
            name = futures[future]
            exc = future.exception()
            if exc:
                logging.error(f"Consensus task for {name} raised an exception: {exc}")
            else:
                logging.info(f"Consensus complete for {name}")

    elapsed = time.time() - start_time
    logging.info(f"Main process completed in {elapsed:.1f}s.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process plasmid sequences")
    parser.add_argument("-e", "--excel_file", required=True, help="Path to the Excel file with plasmid names and unique sequences")
    parser.add_argument("-r", "--ref_dir", default="", help="Directory containing reference plasmid fasta files (optional)")
    parser.add_argument("-f", "--fastq_dir", required=True, help="Directory containing pooled fastq files")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to store output files")
    parser.add_argument("-k", "--keep_temp", action="store_true", help="Keep temporary/intermediate files")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Maximum number of threads to use for external tools (default: 4)")
    parser.add_argument("-q", "--quick_method", action="store_true", help="Use quick alignment & consensus method instead of full assembly")
    args = parser.parse_args()

    # Call main using named arguments (avoids positional mismatches)
    main(
        excel_file=args.excel_file,
        ref_dir=args.ref_dir or "",
        fastq_dir=args.fastq_dir,
        output_dir=args.output_dir,
        keep_temp=args.keep_temp,
        threads=args.threads,
        quick_method=args.quick_method
    )