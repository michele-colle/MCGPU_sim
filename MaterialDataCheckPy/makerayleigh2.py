import numpy as np
from scipy.integrate import cumulative_trapezoid

def generate_mcgpu_balanced_table(input_txt, output_txt, target_points=128):
    # 1. Load native G4 Data
    raw_data = []
    with open(input_txt, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                try:
                    raw_data.append([float(parts[0]), float(parts[1])])
                except ValueError: continue
    
    data = np.array(raw_data)
    
    # 2. Downsample to target_points if necessary (taking every Nth point)
    if len(data) > target_points:
        indices = np.linspace(0, len(data) - 1, target_points).astype(int)
        data = data[indices]
    
    q_g4 = data[:, 0]
    f_val = data[:, 1]
    x_raw = (q_g4 * 20.60744)**2
    f2_raw = f_val**2
    # normalization of f2
    f2_raw = f2_raw 
    n = len(x_raw)

    # 3. Probability Integration
    p_vals = cumulative_trapezoid(f2_raw, x_raw, initial=0)
    p_norm = p_vals / p_vals[-1]
    pdf = f2_raw / p_vals[-1]

    # --- 4. Calculate A and B Parameters (RITA Algorithm) ---
    # Formulas derived from the PenNuc PDF (Eq A.9a and A.9b)
    A = np.zeros(n)
    B = np.zeros(n)
    
    for i in range(n - 1):
        x_i, x_next = x_raw[i], x_raw[i+1]
        P_i, P_next = p_norm[i], p_norm[i+1]
        pdf_i, pdf_next = pdf[i], pdf[i+1]

        dx = x_next - x_i
        dP = P_next - P_i
        
        # Calculate B (Eq A.9a)
        # Note: If pdf is 0, set B=0 to avoid division by zero
        if pdf_i > 1e-30 and pdf_next > 1e-30:
            term = (dP / dx)**2 * (1.0 / (pdf_next * pdf_i))
            B[i] = 1.0 - term
        else:
            B[i] = 0.0
            
        # Calculate A (Eq A.9b)
        if pdf_i > 1e-30:
            A[i] = (dP / dx) * (1.0 / pdf_i) - B[i] - 1.0
        else:
            A[i] = 0.0

    ITL = np.zeros(n, dtype=int)
    ITU = np.zeros(n, dtype=int)
    
    N_intervals = n - 1
    
    for k in range(n):
        # The random number 'ru' for this bin falls in [k/(N-1), (k+1)/(N-1)]
        # However, since 'itn' is an integer truncation, 'ru' can be anywhere 
        # that maps to 'k'.
        # We need the widest possible range of X-indices that could contain 
        # a probability P in the range [k/(N-1), (k+1)/(N-1)].
        
        prob_lower = float(k) / N_intervals
        prob_upper = float(k + 1) / N_intervals
        
        # Find index i such that P[i] <= prob_lower
        # This is our Lower Bound
        # np.searchsorted finds the first index where P > value
        idx_low = np.searchsorted(p_norm, prob_lower, side='right') - 1
        idx_low = max(0, idx_low) # Clamp to 0
        
        # Find index j such that P[j] >= prob_upper
        # This is our Upper Bound
        idx_high = np.searchsorted(p_norm, prob_upper, side='left')
        idx_high = min(n - 1, idx_high) # Clamp to max index
        
        # MC-GPU uses 1-based indexing in the file
        ITL[k] = idx_low + 1
        ITU[k] = idx_high + 1

    # 6. Write File in exact MC-GPU format
    with open(output_txt, 'w') as f:
        f.write("#[RAYLEIGH SAMPLING DATA: BALANCED G4 PENELOPE]\n")
        f.write(f"    {n}\n")
        f.write("#[ X | P | A | B | ITL | ITU ]\n")
        for i in range(n):
            f.write(f" {x_raw[i]:16.10E} {p_norm[i]:16.10E} {A[i]:16.10E} "
                    f"{B[i]:16.10E} {int(ITL[i]):4} {int(ITU[i]):4}\n")

if __name__ == "__main__":
    generate_mcgpu_balanced_table("RayleighFF_G4_Al.txt", "rayleigh_Al_rita_final.txt", 128)