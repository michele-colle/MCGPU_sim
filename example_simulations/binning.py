import numpy as np
import gzip
import os

# --- Parameters ---
ORIGINAL_SHAPE = (1140, 2415, 1740) # (Nz, Ny, Nx)
ORIGINAL_DTYPE = np.uint8
INPUT_FILE = 'phantom/Graff_scattered_22183101.raw.gz'
BINNING_FACTOR = 2

# --- Calculate New Dimensions (Width x Height x Length) ---
Nz, Ny, Nx = ORIGINAL_SHAPE
b = BINNING_FACTOR

# Calculate new size after truncation and binning
dim_z = (Nz - (Nz % b)) // b # Length (870)
dim_y = (Ny - (Ny % b)) // b # Height (1207)
dim_x = (Nx - (Nx % b)) // b # Width (570)

# Create the new filename string in the requested format: width x height x length
dimensions_string = f"{dim_x}x{dim_y}x{dim_z}"
OUTPUT_FILE = f'phantom/Graff_binned_median_{dimensions_string}.raw.gz'

print(f"New dimensions (width x height x length): {dimensions_string}")
print(f"New output filename: {OUTPUT_FILE}")

# --- Execution ---
print(f"Loading {INPUT_FILE}...")
try:
    with gzip.open(INPUT_FILE, 'rb') as f:
        data_flat = np.frombuffer(f.read(), dtype=ORIGINAL_DTYPE)
        original_phantom = data_flat.reshape(ORIGINAL_SHAPE, order='C')
    print("Loading complete.")
except Exception as e:
    print(f"ERROR: Could not load the input file. Ensure '{INPUT_FILE}' exists. Details: {e}")
    exit()


# --- Perform Vectorized Median Binning ---
print(f"Performing Median binning...")
truncated_phantom = original_phantom[:dim_z*b, :dim_y*b, :dim_x*b]
# ... (Previous code remains the same up to 'truncated_phantom')

# --- 1. Pre-allocate the Final Binned Array ---
# Create the final result array (870, 1207, 570) and initialize it
# This array will hold the results from each chunk.
binned_phantom = np.empty((dim_z, dim_y, dim_x), dtype=ORIGINAL_DTYPE)

# --- 2. Define Chunking Parameters ---
CHUNK_SIZE_Y = 100  # Process 100 slices of the Y-axis at a time
# --- Mode Calculation Function (Majority Vote) ---
def find_mode(arr):
    """Calculates the mode (most frequent value) of a 1D array."""
    # Find unique material IDs and their counts
    unique_values, counts = np.unique(arr, return_counts=True)
    
    # Return the unique value corresponding to the maximum count (the mode)
    # np.argmax handles ties by picking the first occurrence.
    return unique_values[np.argmax(counts)]
num_chunks = int(np.ceil(dim_y / CHUNK_SIZE_Y))

print(f"Processing in {num_chunks} chunks of size {CHUNK_SIZE_Y} along Y-axis to save memory...")

# --- 3. Process the Array in Chunks ---
for i in range(num_chunks):
    start_y = i * CHUNK_SIZE_Y
    end_y = min((i + 1) * CHUNK_SIZE_Y, dim_y)
    
    # Define the chunk boundary for the full data (remember to multiply by 'b')
    # Full data chunk dimensions: (dim_z*b, CHUNK_SIZE_Y*b, dim_x*b)
    start_y_full = start_y * b
    end_y_full = end_y * b
    
    # 3a. Extract the chunk (Memory-friendly view, not a copy yet)
    # The chunk is extracted from the original, already truncated array
    current_chunk = truncated_phantom[:, start_y_full:end_y_full, :]

    # --- Vectorized Binning for the Current Chunk ---
    
    # 3b. Reshape to 6D
    # (dim_z, 2, chunk_size_y, 2, dim_x, 2)
    chunk_size_y_actual = end_y - start_y
    reshaped_6d = current_chunk.reshape(dim_z, b, chunk_size_y_actual, b, dim_x, b)

    # 3c. Permute axes
    # (dim_z, chunk_size_y, dim_x, 2, 2, 2)
    transposed = reshaped_6d.transpose(0, 2, 4, 1, 3, 5)

    # 3d. Collapse block (Intermediate array created here, but small!)
    # (dim_z, chunk_size_y, dim_x, 8)
    collapsed = transposed.reshape(dim_z, chunk_size_y_actual, dim_x, b * b * b)

    # 3e. Calculate the median along the 8-element axis (AXIS=3 FIX APPLIED)
    # The temporary memory spike is now limited to a chunk size:
    # 870 * 100 * 570 * 8 * 8 bytes (float64) ~= 317 MB (manageable)
    median_chunk = np.max(collapsed, axis=3).astype(ORIGINAL_DTYPE)

    # 3f. Store the result in the pre-allocated array
    binned_phantom[:, start_y:end_y, :] = median_chunk
    
    # 3g. Explicitly delete temporary arrays to aggressively free memory
    del current_chunk, reshaped_6d, transposed, collapsed, median_chunk
    print(f"Chunk {i+1}/{num_chunks} processed and memory released.")

# --- Final Step: Save the New Phantom File ---
print(f"Binning complete. Final shape (Z x Y x X): {binned_phantom.shape}")
# ... (Save code goes here, using binned_phantom)

# --- Save the New Phantom File ---
print(f"Saving binned phantom to {OUTPUT_FILE}...")
os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
with gzip.open(OUTPUT_FILE, 'wb') as f:
    f.write(binned_phantom.tobytes(order='C'))
print(f"Successfully created: {OUTPUT_FILE}")