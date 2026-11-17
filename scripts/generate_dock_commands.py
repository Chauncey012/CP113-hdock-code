#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import logging
from datetime import datetime

# ==============================================================================
# Variable Definition Section
# ==============================================================================
# Default log file name
LOG_FILE_NAME = "docking_prep.log"
# Default output file for commands
DEFAULT_COMMAND_FILE = "dock_commands.sh"
# Default file to record completed jobs
DEFAULT_DONE_FILE = "completed_jobs.log"

# ==============================================================================
# Main Logic
# ==============================================================================
def setup_logging(log_file):
    """Sets up the logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        filename=log_file,
        filemode='w'
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

def get_completed_jobs(done_file):
    """Reads the list of completed ligand files."""
    if not os.path.exists(done_file):
        return set()
    with open(done_file, 'r') as f:
        return {line.strip() for line in f}

def create_docking_commands(receptor_file, ligand_dir, output_dir, command_file, done_file):
    """
    Generates HDOCK commands in a sorted order based on subdirectory names.
    """
    start_time = datetime.now()
    setup_logging(LOG_FILE_NAME)
    logging.info("Script started.")
    logging.info(f"Receptor: {receptor_file}")
    logging.info(f"Ligand Directory: {ligand_dir}")
    logging.info(f"Output Directory: {output_dir}")
    logging.info("Assuming 'hdock' and 'createpl' are in the system's PATH.")

    receptor_abs_path = os.path.abspath(receptor_file)
    ligand_abs_dir = os.path.abspath(ligand_dir)
    output_abs_dir = os.path.abspath(output_dir)
    done_file_abs_path = os.path.abspath(done_file)
    
    hdock_executable = "hdock"
    createpl_executable = "createpl"

    os.makedirs(output_abs_dir, exist_ok=True)

    completed_jobs = get_completed_jobs(done_file_abs_path)
    logging.info(f"Found {len(completed_jobs)} completed jobs in '{done_file}'.")

    commands_to_run = []
    total_ligands_found = 0
    skipped_jobs = 0

    # --- MODIFICATION START: Sort directories before processing ---
    try:
        # Get all entries in the ligand directory and filter for directories
        subdirs = [d for d in os.listdir(ligand_abs_dir) if os.path.isdir(os.path.join(ligand_abs_dir, d))]
        subdirs.sort() # Sort the directory list alphabetically/numerically
        logging.info(f"Found and sorted {len(subdirs)} subdirectories. Processing in order...")
    except FileNotFoundError:
        logging.error(f"Ligand directory not found: {ligand_abs_dir}")
        return

    # Loop through the sorted subdirectories
    for subdir_name in subdirs:
        current_subdir_path = os.path.join(ligand_abs_dir, subdir_name)
        
        # Now get all PDB files within this specific subdirectory
        for filename in os.listdir(current_subdir_path):
            if filename.endswith(".pdb"):
                total_ligands_found += 1
                
                if filename in completed_jobs:
                    skipped_jobs += 1
                    continue

                # --- Path Construction ---
                ligand_abs_path = os.path.join(current_subdir_path, filename)
                output_subdir = os.path.join(output_abs_dir, subdir_name)
                
                receptor_rel_path = os.path.relpath(receptor_abs_path, start=output_subdir)
                ligand_rel_path = os.path.relpath(ligand_abs_path, start=output_subdir)

                # --- Job and Output Naming ---
                ligand_name = os.path.splitext(filename)[0]
                receptor_name = os.path.splitext(os.path.basename(receptor_abs_path))[0]
                job_id = f"{receptor_name}_{ligand_name}"
                
                output_hdock_file_local = f"{job_id}.out"
                output_pdb_file_local = f"{job_id}.pdb"

                # --- Build the Command ---
                command = (
                    f"mkdir -p {output_subdir} && "
                    f"cd {output_subdir} && "
                    f"{hdock_executable} {receptor_rel_path} {ligand_rel_path} -out {output_hdock_file_local} && "
                    f"{createpl_executable} {output_hdock_file_local} {output_pdb_file_local} -nmax 1 -complex && "
                    f"echo '{filename}' >> {done_file_abs_path}"
                )
                commands_to_run.append(command)
    # --- MODIFICATION END ---

    logging.info(f"Total ligands found: {total_ligands_found}")
    logging.info(f"Skipped completed jobs: {skipped_jobs}")
    logging.info(f"Generated {len(commands_to_run)} new commands.")

    with open(command_file, 'w') as f:
        for cmd in commands_to_run:
            f.write(cmd + '\n')
    
    logging.info(f"All commands have been written to '{command_file}' in sorted order.")
    end_time = datetime.now()
    logging.info(f"Script finished. Total running time: {end_time - start_time}.")
    
    print(f"\n✅ Generation complete. Commands are sorted by directory.")
    print(f"📄 Log file saved to: {LOG_FILE_NAME}")
    print(f"🚀 Command file for parallel execution saved to: {command_file}")
    print(f"To run the docking jobs, use the following command:")
    print(f"   parallel -j [number_of_cores] < {command_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate HDOCK docking commands in a sorted order for parallel execution.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Example:
  python3 %(prog)s \\
    -r MC1R.pdb \\
    -l /path/to/uniprot \\
    -o /path/to/uniprot_results

Make sure 'hdock' and 'createpl' are in your system's PATH.
"""
    )
    parser.add_argument("-r", "--receptor_file", required=True, help="Path to the receptor PDB file.")
    parser.add_argument("-l", "--ligand_dir", required=True, help="Path to the top-level directory containing ligand subdirectories.")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the top-level directory where results will be saved.")
    parser.add_argument("-cf", "--command_file", default=DEFAULT_COMMAND_FILE, help=f"Output file name for the shell commands (default: {DEFAULT_COMMAND_FILE}).")
    parser.add_argument("-df", "--done_file", default=DEFAULT_DONE_FILE, help=f"File to track completed jobs (default: {DEFAULT_DONE_FILE}).")

    args = parser.parse_args()
    create_docking_commands(args.receptor_file, args.ligand_dir, args.output_dir, args.command_file, args.done_file)

if __name__ == "__main__":
    main()
