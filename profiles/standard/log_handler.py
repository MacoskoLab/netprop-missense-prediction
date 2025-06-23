#!/usr/bin/env python3

"""
Snakemake log handler
Automatically logs workflow execution to timestamped files
"""

import os
import sys
from datetime import datetime


def main():
    # Create logs directory if it doesn't exist
    os.makedirs("../logs", exist_ok=True)

    # Generate timestamp for log file
    timestamp = datetime.now().strftime("%Y-%m-%dT%H-%M-%S%z")
    log_file = f"../logs/snakemake.{timestamp}.log"

    # Open log file for writing
    with open(log_file, "w") as f:
        # Write header
        f.write(f"# Snakemake execution log\n")
        f.write(f"# Started: {datetime.now().isoformat()}\n")
        f.write(f"# Command: {' '.join(sys.argv)}\n\n")

        # Read from stdin and write to both stdout and log file
        for line in sys.stdin:
            sys.stdout.write(line)
            f.write(line)
            f.flush()


if __name__ == "__main__":
    main()
