import matplotlib.pyplot as plt
import numpy as np

def load_rayleigh_file(filepath):
    """
    Parses an MC-GPU style Rayleigh data file.
    Returns X and P as numpy arrays.
    """
    x_vals = []
    p_vals = []
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                # Skip comments and header lines
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.split()
                # A valid data line must have at least X and P (first two columns)
                if len(parts) >= 6:
                    x_vals.append(float(parts[0]))
                    p_vals.append(float(parts[1]))
                    
        return np.array(x_vals), np.array(p_vals)
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None, None

def plot_comparison(file1, file2):
    # Load data
    x1, p1 = load_rayleigh_file(file1)
    x2, p2 = load_rayleigh_file(file2)
    
    if x1 is None or x2 is None:
        return

    plt.figure(figsize=(10, 6))
    
    # Plot the Cumulative Probability (P) vs Momentum Transfer (X)
    plt.plot(x1, p1, 'o-', label=f'File 1: {file1}', markersize=4, alpha=0.8)
    plt.plot(x2, p2, 's-', label=f'File 2: {file2}', markersize=4, alpha=0.8)
    
    # In physics, the slope dP/dX is the actual Form Factor squared (F^2)
    plt.title('Comparison of Rayleigh Cumulative Probability (P) vs $X = (q/4\pi)^2$')
    plt.xlabel('Momentum Transfer Variable (X)')
    plt.ylabel('Cumulative Probability (P)')
    
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    
    # Annotate the "forward peaking" area
    plt.annotate('Small-angle scattering region', xy=(0.05, 0.2), xytext=(0.2, 0.4),
                 arrowprops=dict(facecolor='black', shrink=0.05))

    plt.tight_layout()
    
    # Save the plot
    output_plot = "rayleigh_comparison.png"
    plt.savefig(output_plot)
    print(f"Plot saved as: {output_plot}")
    plt.show()

if __name__ == "__main__":
    # Ensure you have generated these files with the previous script first
    file_a = "rayleigh_Al.txt"
    file_b = "rayleigh_AL_original.txt"  # Or another element you've generated
    
    import os
    if os.path.exists(file_a) and os.path.exists(file_b):
        plot_comparison(file_a, file_b)
    else:
        print(f"Missing data files. Please ensure {file_a} and {file_b} exist.")