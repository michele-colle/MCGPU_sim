import matplotlib.pyplot as plt
import numpy as np
import os

def load_mcgpu_material(filepath):
    """
    Parses ONLY the Mean Free Path section of an MC-GPU material file.
    Stops reading when it hits the next section header (starts with #).
    """
    data = {'energy': [], 'rayleigh': [], 'compton': [], 'photo': [], 'total': []}
    reading_data = False
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                stripped = line.strip()
                
                # If we encounter a '#' after we've already started reading data, STOP.
                if reading_data and stripped.startswith('#'):
                    break
                
                # Skip header comments at the very top
                if stripped.startswith('#') or not stripped:
                    continue
                
                parts = stripped.split()
                if len(parts) >= 5:
                    try:
                        data['energy'].append(float(parts[0]))
                        data['rayleigh'].append(float(parts[1]))
                        data['compton'].append(float(parts[2]))
                        data['photo'].append(float(parts[3]))
                        data['total'].append(float(parts[4]))
                        reading_data = True # Flag that we are now in the data block
                    except ValueError:
                        # In case of malformed lines or mid-table headers
                        continue
                        
        return {k: np.array(v) for k, v in data.items()}
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def plot_individual_canvases(file1_path, file2_path):
    data1 = load_mcgpu_material(file1_path)
    data2 = load_mcgpu_material(file2_path)
    
    if data1 is None or data2 is None:
        return

    # Configuration for the 4 plots
    plot_configs = [
        ('rayleigh', 'Rayleigh Mean Free Path', 'blue'),
        ('compton', 'Compton Mean Free Path', 'green'),
        ('photo', 'Photoelectric Mean Free Path', 'red'),
        ('total', 'Total Mean Free Path (with Pair Prod)', 'black')
    ]

    for key, title, color in plot_configs:
        # Create a brand new figure for each interaction
        plt.figure(figsize=(9, 6))
        
        plt.plot(data1['energy'], data1[key], label=f'File 1 ({os.path.basename(file1_path)})', 
                 color=color, linewidth=2, alpha=0.7)
        plt.plot(data2['energy'], data2[key], label=f'File 2 ({os.path.basename(file2_path)})', 
                 color='orange', linestyle='--', linewidth=2)
        
        plt.title(f'{title} vs Energy', fontsize=14)
        plt.xlabel('Energy (eV)', fontsize=12)
        plt.ylabel('Mean Free Path (cm)', fontsize=12)
        
        # Physics data is almost always best viewed on Log-Log scales
        plt.xscale('log')
        plt.yscale('log')
        
        plt.grid(True, which="both", ls="-", alpha=0.3)
        plt.legend()
        
        # Save each one separately
        filename = f"comparison_{key}.png"
        plt.savefig(filename)
        print(f"Created canvas for {key}: {filename}")

    # Display all windows at once
    plt.show()

if __name__ == "__main__":
    # Replace these with your actual filename paths
    f1 = "Al_5-120keV_xraylib.mcgpu"
    f2 = "Al__5-120keV.mcgpu"
    
    if os.path.exists(f1) and os.path.exists(f2):
        plot_individual_canvases(f1, f2)
    else:
        print("Error: Ensure both material files exist before running.")