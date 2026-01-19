import numpy as np

def load_rayleigh_file(filepath):
    """
    Parses an MC-GPU style Rayleigh data file.
    Returns X and P as numpy arrays.
    """
    x_vals = []
    p_vals = []
    a_vals = []
    b_vals = []
    itl_vals = []
    itu_vals = []
    
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
                    a_vals.append(float(parts[2]))
                    b_vals.append(float(parts[3]))
                    itl_vals.append(int(parts[4]))
                    itu_vals.append(int(parts[5]))
                    
        return np.array(x_vals), np.array(p_vals), np.array(a_vals), np.array(b_vals), np.array(itl_vals), np.array(itu_vals)
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None, None, None, None, None, None
    
def get_smooth_rita_curve(x_co, p_co, aco, bco, resolution=10):
    """
    Evaluates the RITA formula at 'resolution' points between each tabulated bin
    to reconstruct the smooth curve used by the GPU.
    """
    smooth_x = []
    smooth_p = []
    
    # Iterate through each bin (except the last point)
    for i in range(len(x_co) - 1):
        d = p_co[i+1] - p_co[i]
        x_start = x_co[i]
        x_range = x_co[i+1] - x_co[i]
        
        # Sub-sample the probability 'rr' within this bin
        # rr goes from 0 to d
        for rr in np.linspace(0, d, resolution):
            if d > 1e-18:
                # The exact formula from the CUDA snippet:
                numerator = (aco[i] + 1.0 + bco[i]) * d * rr
                denominator = (d*d + (aco[i]*d + bco[i]*rr) * rr)
                
                # Calculate the interpolated x value
                xx = x_start + (numerator / denominator) * x_range
                smooth_x.append(xx)
            else:
                smooth_x.append(x_start)
            smooth_p.append(p_co[i] + rr)

    return np.array(smooth_x), np.array(smooth_p)

def evaluate_rita_gpu_emulator(ru, x_co, p_co, aco, bco, itl, itu):
    """
    Simulates the MC-GPU C++ kernel logic to evaluate the curve at probability 'ru'.
    """
    NP = len(x_co)
    # 1. Map probability to the lookup interval (itn)
    itn = int(ru * (NP - 1))
    if itn >= NP: itn = NP - 1
    
    # 2. Retrieve the search boundaries from the ITL/ITU tables
    # (Subtract 1 because the file uses 1-based Fortran indexing)
    i_idx = itl[itn]
    j_idx = itu[itn]
    
    # 3. Binary search within the narrowed range (Matches CUDA snippet)
    if (j_idx - i_idx) > 1:
        while (j_idx - i_idx) > 1:
            k = (i_idx + j_idx) // 2
            if ru > p_co[k - 1]:
                i_idx = k
            else:
                j_idx = k
    
    # 4. Final bin index (0-based)
    idx = i_idx - 1
    
    # 5. RITA formula (Rational Inverse Transform)
    rr = ru - p_co[idx]
    if rr > 1e-16:
        d = p_co[idx+1] - p_co[idx]
        numerator = (aco[idx] + 1.0 + bco[idx]) * d * rr
        denominator = (d*d + (aco[idx]*d + bco[idx]*rr) * rr)
        xx = x_co[idx] + (numerator / denominator) * (x_co[idx+1] - x_co[idx])
    else:
        xx = x_co[idx]
        
    return xx