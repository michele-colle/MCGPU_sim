import matplotlib.pyplot as plt
import numpy as np
from rayleigh_functions import *


def plot_comparison(file1, file2):
    # Load data
    x1, p1, a1, b1, itl1, itu1 = load_rayleigh_file(file1)
    x1smooth, dense_p1 = get_smooth_rita_curve(x1, p1, a1, b1, resolution=50)
    x2, p2, a2, b2, itl2, itu2 = load_rayleigh_file(file2)
    x2smooth, dense_p2 = get_smooth_rita_curve(x2, p2, a2, b2, resolution=50)
    
    if x1 is None or x2 is None:
        return

    plt.figure(1)
    
    # Plot the Cumulative Probability (P) vs Momentum Transfer (X)
    plt.plot(x1, p1, 'o', label=f'File 1: {file1}', markersize=4, alpha=0.8)
    plt.plot(x1smooth, dense_p1, '-', label=f'File 1 Smooth: {file1}', color='blue', alpha=0.7)
    plt.plot(x2, p2, 's', label=f'File 2: {file2}', markersize=4, alpha=0.8)
    plt.plot(x2smooth, dense_p2, '-', label=f'File 2 Smooth: {file2}', color='red', alpha=0.7)
    
    # In physics, the slope dP/dX is the actual Form Factor squared (F^2)
    plt.title('Comparison of Rayleigh Cumulative Probability (P) vs $X = (q/4\pi)^2$')
    plt.xlabel('Momentum Transfer Variable (X)')
    plt.ylabel('Cumulative Probability (P)')
    #plt.xlim(0,1e3)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    
    #plt.xlim(671,11351)
    #plt.ylim(0.9999990649435871,1.0000000607496007)
    plt.show()
    
    

if __name__ == "__main__":
    # Ensure you have generated these files with the previous script first
    file_a = "rayleigh_Al_rita_final.txt"
    file_b = "rayleigh_AL_original.txt"  # Or another element you've generated
    
    import os
    if os.path.exists(file_a) and os.path.exists(file_b):
        plot_comparison(file_a, file_b)
    else:
        print(f"Missing data files. Please ensure {file_a} and {file_b} exist.")