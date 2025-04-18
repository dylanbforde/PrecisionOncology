#!/usr/bin/env python
"""
This script runs all data preparation scripts to prepare the data for unsupervised learning.
"""

import os
import subprocess
import time
import pathlib

def main():
    # Get the current script directory
    current_dir = pathlib.Path(__file__).parent.absolute()
    
    scripts = [
        "prepare_mutations.py",
        "prepare_cna.py",
        "prepare_mrna.py",
        "prepare_ancestry.py"
    ]
    
    for script in scripts:
        script_path = os.path.join(current_dir, script)
        
        # Make sure script is executable
        os.chmod(script_path, 0o755)
        
        print(f"\n{'='*50}")
        print(f"Running {script}...")
        print(f"{'='*50}\n")
        
        start_time = time.time()
        
        # Run the script
        subprocess.run(['python', script_path], check=True)
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        print(f"\nCompleted {script} in {elapsed_time:.2f} seconds")
    
    print("\nAll data preparation complete!")

if __name__ == "__main__":
    main()