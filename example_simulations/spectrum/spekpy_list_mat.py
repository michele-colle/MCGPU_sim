import spekpy as sp
import sys
from contextlib import redirect_stdout

# --- Configuration ---
OUTPUT_FILE = 'supported_materials.txt'
# ---------------------

def save_supported_materials_to_file(filename):
    """Initializes SpekPy and redirects the output of s.show_matls() to a file."""
    
    # Initialize SpekPy instance
    s = sp.Spek()
    
    print(f"Attempting to save material list to '{filename}'...")

    try:
        # Open the file for writing ('w')
        with open(filename, 'w') as f:
            # Temporarily redirect standard output (print statements) to the file object 'f'
            with redirect_stdout(f):
                # This function call's output is now written to the file
                s.show_matls() 
        
        print(f"\nSUCCESS: The list of supported filtration materials has been saved to '{filename}'.")

    except Exception as e:
        print(f"\nERROR: Could not save materials list. Reason: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    save_supported_materials_to_file(OUTPUT_FILE)