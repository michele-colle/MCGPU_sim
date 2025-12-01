import numpy as np
import gzip
import os

# --- Parameters ---
# Volume dimensions
VOLUME_WIDTH = 10.0    # cm (X dimension)
VOLUME_HEIGHT = 10.0   # cm (Y dimension)
VOLUME_LENGTH = 10.0   # cm (Z dimension)

# Cylinder dimensions and material
CYLINDER_RADIUS = 2.0  # cm
CYLINDER_HEIGHT = 5.0  # cm
CYLINDER_CENTER_X = 0.0   # cm (0,0,0 = volume center)
CYLINDER_CENTER_Y = 0.0   # cm
CYLINDER_CENTER_Z = 0.0   # cm

# Resolution
VOXEL_SIZE = 0.02    # cm (50 microns, matching typical MC-GPU resolution)

# Material IDs (from VICTRE defaults)
MATERIAL_ID = 7        # 7 = Muscle
BACKGROUND_ID = 0      # 0 = Air
OUTPUT_DTYPE = np.uint8

# --- Calculate Grid Dimensions ---
# Convert cm to number of voxels
num_voxels_x = int(np.ceil(VOLUME_WIDTH / VOXEL_SIZE))
num_voxels_y = int(np.ceil(VOLUME_HEIGHT / VOXEL_SIZE))
num_voxels_z = int(np.ceil(VOLUME_LENGTH / VOXEL_SIZE))

# Create output filename with dimensions
output_filename = f'phantom/cylinder_muscle_{num_voxels_x}x{num_voxels_y}x{num_voxels_z}.raw.gz'
os.makedirs(os.path.dirname(output_filename), exist_ok=True)

# Convert cylinder center from physical coordinates to voxel indices
# Origin (0,0,0) is at the center of the volume
center_voxel_x = num_voxels_x / 2.0 + CYLINDER_CENTER_X / VOXEL_SIZE
center_voxel_y = num_voxels_y / 2.0 + CYLINDER_CENTER_Y / VOXEL_SIZE
center_voxel_z = num_voxels_z / 2.0 + CYLINDER_CENTER_Z / VOXEL_SIZE

# Cylinder extends along Z axis
z_min_voxel = center_voxel_z - (CYLINDER_HEIGHT / (2.0 * VOXEL_SIZE))
z_max_voxel = center_voxel_z + (CYLINDER_HEIGHT / (2.0 * VOXEL_SIZE))
radius_voxels = CYLINDER_RADIUS / VOXEL_SIZE

print(f"Creating cylinder phantom:")
print(f"  Volume: {VOLUME_WIDTH} x {VOLUME_HEIGHT} x {VOLUME_LENGTH} cm")
print(f"  Volume grid: {num_voxels_x} x {num_voxels_y} x {num_voxels_z} voxels")
print(f"  Voxel size: {VOXEL_SIZE} cm")
print(f"  Cylinder center: ({CYLINDER_CENTER_X}, {CYLINDER_CENTER_Y}, {CYLINDER_CENTER_Z}) cm")
print(f"  Cylinder center (voxels): ({center_voxel_x:.1f}, {center_voxel_y:.1f}, {center_voxel_z:.1f})")
print(f"  Cylinder radius: {CYLINDER_RADIUS} cm ({radius_voxels:.1f} voxels)")
print(f"  Cylinder height: {CYLINDER_HEIGHT} cm (Z: {z_min_voxel:.1f} to {z_max_voxel:.1f} voxels)")
print(f"  Material: Muscle (ID={MATERIAL_ID})")
print(f"  Background: Air (ID={BACKGROUND_ID})")

# --- Create the Phantom Array ---
print("\nGenerating phantom geometry...")

# Initialize with background (air)
phantom = np.full((num_voxels_z, num_voxels_y, num_voxels_x), BACKGROUND_ID, dtype=OUTPUT_DTYPE)

# Create the cylinder by checking distance from axis and Z position
for iz in range(num_voxels_z):
    if iz % max(1, num_voxels_z // 10) == 0:  # Print progress every 10%
        print(f"  Processing Z-slice {iz}/{num_voxels_z} ({100*iz/num_voxels_z:.1f}%)")
    
    for iy in range(num_voxels_y):
        for ix in range(num_voxels_x):
            # Check if voxel is within cylinder height (Z dimension)
            if iz >= z_min_voxel and iz <= z_max_voxel:
                # Distance from cylinder axis in XY plane
                dx = ix - center_voxel_x + 0.5
                dy = iy - center_voxel_y + 0.5
                distance = np.sqrt(dx**2 + dy**2)
                
                # Fill cylinder if within radius
                if distance <= radius_voxels:
                    phantom[iz, iy, ix] = MATERIAL_ID

print(f"Phantom generated: {phantom.shape}")

# --- Calculate Material Statistics ---
num_muscle_voxels = np.sum(phantom == MATERIAL_ID)
num_air_voxels = np.sum(phantom == BACKGROUND_ID)
total_voxels = phantom.size

print(f"\nMaterial distribution:")
print(f"  Muscle voxels: {num_muscle_voxels} ({100*num_muscle_voxels/total_voxels:.1f}%)")
print(f"  Air voxels:    {num_air_voxels} ({100*num_air_voxels/total_voxels:.1f}%)")
print(f"  Total voxels:  {total_voxels}")

# --- Save the Phantom ---
print(f"\nSaving phantom to {output_filename}...")
with gzip.open(output_filename, 'wb') as f:
    f.write(phantom.tobytes(order='C'))

print(f"Successfully created: {output_filename}")

# --- Generate Corresponding MC-GPU Input File Snippet ---
input_snippet_filename = os.path.join(os.path.dirname(output_filename), 'input_cylinder_snippet.txt')
with open(input_snippet_filename, 'w') as f:
    f.write("# --- Cylinder Phantom Configuration Snippet ---\n")
    f.write("# Add these lines to your MC-GPU input file:\n\n")
    f.write("#[SECTION VOXELIZED GEOMETRY FILE v.2017-07-26]\n")
    f.write(f"{output_filename}     # VOXEL GEOMETRY FILE\n")
    f.write(f" 0.0    0.0    0.0              # OFFSET OF THE VOXEL GEOMETRY [cm]\n")
    f.write(f" {num_voxels_x}   {num_voxels_y}   {num_voxels_z}                 # NUMBER OF VOXELS: Nx Ny Nz\n")
    f.write(f" {VOXEL_SIZE} {VOXEL_SIZE} {VOXEL_SIZE}           # VOXEL SIZES [cm]\n")
    f.write(f" 1 1 1                          # SIZE OF LOW RESOLUTION VOXELS (binary tree) [powers of two]\n\n")
    f.write("#[SECTION MATERIAL FILE LIST v.2020-03-03]\n")
    f.write("material/air__5-120keV.mcgpu.gz                 density=0.0012   voxelId=0\n")
    f.write("material/muscle__5-120keV.mcgpu.gz              density=1.050    voxelId=7\n")

print(f"Input file snippet saved to: {input_snippet_filename}")
print("\nTo use this phantom in MC-GPU, copy the contents of the snippet file into your input file.")
