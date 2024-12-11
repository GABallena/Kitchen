import os
import subprocess
import argparse
import time

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Trimming, cleaning, and quality checking script.")
    parser.add_argument("--fastqc_output_dir", required=True, help="Directory for FastQC output.")
    parser.add_argument("--adapter_file", required=True, help="Adapter sequences file.")
    parser.add_argument("--sample", required=True, help="Sample name.")
    parser.add_argument("--qc_report", required=True, help="Output QC report file.")
    parser.add_argument("--retries", type=int, default=3, help="Number of retries for each step.")
    return parser.parse_args()

def run_command(command, description, retries=3):
    """Run a system command with retries."""
    for attempt in range(1, retries + 1):
        try:
            print(f"Attempt {attempt} for {description}...")
            subprocess.run(command, shell=True, check=True)
            print(f"Success: {description}")
            return
        except subprocess.CalledProcessError as e:
            print(f"Error during {description}: {e}")
            if attempt < retries:
                print(f"Retrying...")
                time.sleep(2)
            else:
                raise RuntimeError(f"Failed {description} after {retries} attempts.") from e

def run_trimmomatic(args, seq):
    """Run Trimmomatic for trimming reads."""
    command = (
        f"trimmomatic PE -threads 4 {seq['input_R1']} {seq['input_R2']} "
        f"{seq['output_paired1']} {seq['unpaired1']} {seq['output_paired2']} {seq['unpaired2']} "
        f"ILLUMINACLIP:{args.adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36"
    )
    run_command(command, f"Trimmomatic for {seq['name']}", retries=args.retries)

def run_fastqc(trimmed_file, output_dir, retries=3):
    """Run FastQC on a given file."""
    command = f"fastqc --threads 4 --outdir {output_dir} {trimmed_file}"
    run_command(command, f"FastQC for {trimmed_file}", retries=retries)

def parse_fastqc_summary(fastqc_summary_file):
    """Parse FastQC summary file to check for quality metrics."""
    required_metrics = ["Per base sequence quality", "Adapter Content"]
    passed_metrics = {metric: False for metric in required_metrics}

    try:
        with open(fastqc_summary_file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) > 1:
                    status, metric = fields[0], fields[1]
                    if metric in required_metrics:
                        passed_metrics[metric] = status == "PASS"
    except FileNotFoundError:
        print(f"Error: Summary file {fastqc_summary_file} not found.")
        return False

    if all(passed_metrics.values()):
        print(f"FastQC passed for {fastqc_summary_file}")
        return True
    else:
        print(f"FastQC failed for {fastqc_summary_file}. Details:")
        for metric, passed in passed_metrics.items():
            if not passed:
                print(f"- {metric}: FAILED")
        return False

def process_sequences(sequences, args, max_retries):
    """
    Process a list of sequences. Retry failed sequences up to max_retries.
    """
    attempt = 0

    while attempt < max_retries and sequences:
        attempt += 1
        print(f"Processing attempt {attempt} for {len(sequences)} sequences...")

        remaining_sequences = []

        for seq in sequences:
            try:
                print(f"Processing sequence: {seq['name']}")

                # Step 1: Run Trimmomatic
                run_trimmomatic(args, seq)

                # Step 2: Run FastQC for paired reads
                run_fastqc(seq['output_paired1'], args.fastqc_output_dir, retries=args.retries)
                run_fastqc(seq['output_paired2'], args.fastqc_output_dir, retries=args.retries)

                # Step 3: Verify FastQC metrics
                summary_file_1 = os.path.join(
                    args.fastqc_output_dir, os.path.basename(seq['output_paired1']) + "_fastqc/summary.txt"
                )
                summary_file_2 = os.path.join(
                    args.fastqc_output_dir, os.path.basename(seq['output_paired2']) + "_fastqc/summary.txt"
                )

                if not (parse_fastqc_summary(summary_file_1) and parse_fastqc_summary(summary_file_2)):
                    raise RuntimeError(f"Quality check failed for sequence {seq['name']}.")

                print(f"Sequence {seq['name']} processed successfully.")

            except Exception as e:
                print(f"Error processing sequence {seq['name']}: {e}")
                remaining_sequences.append(seq)

        sequences = remaining_sequences

    if sequences:
        print(f"Failed to process {len(sequences)} sequences after {max_retries} attempts.")
    else:
        print("All sequences processed successfully.")

def main():
    """
    Main function for trimming, cleaning, and QC checks.
    """
    args = parse_arguments()

    # Example sequence structure
    sequences = [
        {
            "name": f"Sequence_{i}",
            "input_R1": f"input_R1_seq_{i}.fastq",
            "input_R2": f"input_R2_seq_{i}.fastq",
            "output_paired1": f"paired_R1_seq_{i}.fastq",
            "output_paired2": f"paired_R2_seq_{i}.fastq",
            "unpaired1": f"unpaired_R1_seq_{i}.fastq",
            "unpaired2": f"unpaired_R2_seq_{i}.fastq",
        }
        for i in range(1, 11)  # Adjust to your actual number of sequences
    ]

    max_retries = 3
    process_sequences(sequences, args, max_retries)

if __name__ == "__main__":
    main()
