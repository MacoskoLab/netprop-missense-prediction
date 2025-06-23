"""
Pipeline run management utilities
"""

import os
import shutil
from datetime import datetime


def setup_logging():
    """Setup logging to file with timestamp"""
    timestamp = datetime.now().strftime("%Y-%m-%dT%H-%M-%S%z")
    log_file = f"logs/snakemake.{timestamp}.log"
    os.makedirs("logs", exist_ok=True)
    return log_file


def log_message(message, log_file=None):
    """Log message to both stdout and log file"""
    print(message, flush=True)
    if log_file:
        with open(log_file, "a") as f:
            f.write(f"{datetime.now().isoformat()}: {message}\n")


def create_run_name(config):
    """Create a run name based on the run number and config."""
    run_name = config.get("run_name", None)
    if not run_name:
        run_name = "_"

    run_number = config.get("run_number", None)
    if not run_number:
        # List existing run directories
        existing_runs = [
            run_dir
            for run_dir in os.listdir("results")
            if os.path.isdir(os.path.join("results", run_dir))
        ]
        # Filter for directories with the same run name
        existing_runs = [
            run_dir for run_dir in existing_runs if run_dir.startswith(run_name)
        ]
        # Extract run numbers
        run_numbers = [int(run_dir.split("_")[-1]) for run_dir in existing_runs]
        run_number = max(run_numbers) + 1 if run_numbers else 1

    return run_name + "_" + str(run_number)


def copy_config_to_run(run):
    """Copy the config file to the run directory."""
    # Copy config file to run directory
    run_dir = f"results/{run}"
    os.makedirs(run_dir, exist_ok=True)
    shutil.copyfile("config/config.yml", f"{run_dir}/{run}_config.yml")
