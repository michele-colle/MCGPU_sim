import matplotlib.pyplot as plt
from rayleigh_functions import *
import numpy as np



def plot_comparison(file1):
    # Load data
    x, p, a, b, itl, itu = load_rayleigh_file(file1)
    smooth_x_values, dense_p = get_smooth_rita_curve(x, p, a, b, resolution=50)

    # Generate a high-resolution probability grid (2000 points)
    ru_samples = np.linspace(0.99999, 1, 200000)**2

    # Evaluate the curve using our ITL/ITU lookup logic
    smooth_x2 = [evaluate_rita_gpu_emulator(r, x, p, a, b, itl, itu) for r in ru_samples]

    plt.figure(figsize=(10, 6))
    
    # Plot the Cumulative Probability (P) vs Momentum Transfer (X)
    plt.plot(x, p, 'o', label=f'File 1: {file1}', markersize=4, alpha=0.8)
    plt.plot(smooth_x_values, dense_p, '-', label='RITA Smooth Curve', color='red', alpha=0.7)
    plt.plot(smooth_x2, ru_samples, 'x-', label='RITA Smooth Curve', color='red', alpha=0.7)
    
    # In physics, the slope dP/dX is the actual Form Factor squared (F^2)
    plt.title('Comparison of Rayleigh Cumulative Probability (P) vs $X = (q/4\pi)^2$')
    plt.xlabel('Momentum Transfer Variable (X)')
    plt.ylabel('Cumulative Probability (P)')
    
    # plt.xlim(671,11351)
    # plt.ylim(0.9999990649435871,1.0000000607496007)
    
    #plt.xlim(8.172216057064375,8.172216057064375)
    #plt.ylim(0.9775596410445989,1.0037529366735076)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    

    plt.tight_layout()
    
    plt.show()

if __name__ == "__main__":
    # Ensure you have generated these files with the previous script first
    file_a = "rayleigh_AL_original.txt"
    file_a = "rayleigh_Al_rita_final.txt"
    
    import os
    if os.path.exists(file_a):
        plot_comparison(file_a)
    else:
        print(f"Missing data files. Please ensure {file_a} exist.")