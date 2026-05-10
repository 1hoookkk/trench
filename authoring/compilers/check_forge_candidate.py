#!/usr/bin/env python3
import sys
import subprocess
from pathlib import Path

def main():
    # Simple wrapper to match the user's requested CLI signature:
    # --frame-map (ignored but accepted) --compiled --plot --recipe
    tools_dir = Path(__file__).parent
    validator = tools_dir / "validate_filter_candidate.py"
    
    cmd = [sys.executable, str(validator)]
    
    # Forward relevant args, stripping --frame-map which validate_filter_candidate doesn't take
    skip_next = False
    for i, arg in enumerate(sys.argv[1:]):
        if skip_next:
            skip_next = False
            continue
        if arg == "--frame-map":
            skip_next = True
            continue
        cmd.append(arg)
        
    result = subprocess.run(cmd)
    sys.exit(result.returncode)

if __name__ == "__main__":
    main()
