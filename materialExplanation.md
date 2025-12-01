# MC-GPU Material Management Guide

## Overview

MC-GPU v1.5b supports two modes for managing materials in simulations:
1. **Legacy Mode**: Uses hardcoded VICTRE defaults (old input files)
2. **New Mode**: Allows explicit density and voxel ID specification (modern input files)

---

## Two Input Modes

### Mode 1: Legacy Format (Without density= and voxelId=)

When using the old format like in `example_simulations/`:

```
#[SECTION MATERIAL FILE LIST v.2009-11-30]
material/air__5-120keV.mcgpu.gz                  #  1st MATERIAL FILE (.gz accepted)
material/adipose__5-120keV.mcgpu.gz              #  2nd MATERIAL FILE
material/skin__5-120keV.mcgpu.gz 
material/glandular__5-120keV.mcgpu.gz
```

The code automatically uses **hardcoded VICTRE default mappings** for voxel IDs and densities. The materials are assigned in order to specific voxel IDs.

### Mode 2: New Format (With density= and voxelId=)

```
#[SECTION MATERIAL FILE LIST v.2020-03-03]
material/air__5-120keV.mcgpu.gz                  density=0.0012   voxelId=0
material/adipose__5-120keV.mcgpu.gz              density=0.920    voxelId=1
material/mynewmaterial__5-120keV.mcgpu.gz        density=1.500    voxelId=100
```

When you include **ANY** `voxelId=` in the input file, the code **abandons the hardcoded defaults** and uses only what you explicitly specify.

---

## Default VICTRE Voxel-to-Material Mapping

The following hardcoded mappings are used when no `voxelId=` keywords are present:

| Voxel ID | Material Name | Material # | Default Density |
|----------|---------------|-----------|-----------------|
| 0 | Air | 1 | 0.00120 g/cm³ |
| 1 | Adipose/Fat | 2 | 0.920 g/cm³ |
| 2 | Skin | 3 | 1.090 g/cm³ |
| 29 | Glandular | 4 | 1.035 g/cm³ |
| 33 | Nipple | 5 | 1.090 g/cm³ |
| 40 | Muscle | 7 | 1.050 g/cm³ |
| 50 | Compression Paddle | 11 | 1.060 g/cm³ |
| 65 | Tungsten | 14 | 19.30 g/cm³ |
| 66 | Selenium | 15 | 4.50 g/cm³ |
| 88 | Ligament | 6 | 1.120 g/cm³ |
| 95 | Terminal Duct Lobular Unit | 10 | 1.050 g/cm³ |
| 125 | Duct | 8 | 1.050 g/cm³ |
| 150 | Artery | 9 | 1.000 g/cm³ |
| 200 | Mass/Signal | 12 | 1.060 g/cm³ |
| 225 | Vein | 9 | 1.000 g/cm³ |
| 250 | Microcalcification | 13 | 1.781 g/cm³ |

---

## How to Add a New Material

### Step 1: Obtain or Create a Material File

Material files must be in `.mcgpu` format (can be gzipped as `.mcgpu.gz`). Available materials in `example_simulations/material/`:

- `air__5-120keV.mcgpu.gz`
- `adipose__5-120keV.mcgpu.gz`
- `glandular__5-120keV.mcgpu.gz`
- `muscle__5-120keV.mcgpu.gz`
- `blood__5-120keV.mcgpu.gz`
- `skin__5-120keV.mcgpu.gz`
- `connective_Woodard__5-120keV.mcgpu.gz`
- `polystyrene__5-120keV.mcgpu.gz`
- `CalciumOxalate__5-120keV.mcgpu.gz`
- `W__5-120keV.mcgpu.gz` (Tungsten)
- `Se__5-120keV.mcgpu.gz` (Selenium)

#### Creating Custom Materials

**Option A: Use Existing Materials with Different Densities** (Easiest)
Simply use an existing material file with a custom density value in your input file. This is the fastest approach for material variations.

**Option B: Create from PENELOPE Database**
1. The utility `MC-GPU_create_material_data.f` (FORTRAN) converts PENELOPE 2006 `.mat` files to MC-GPU format
2. Found in the original MC-GPU repository: https://github.com/DIDSR/MCGPU
3. Compile with: `gfortran MC-GPU_create_material_data.f -o create_material_data`
4. Run: `./create_material_data < penelope_material.mat > output.mcgpu`
5. Gzip the output: `gzip output.mcgpu`

**Option C: Obtain Pre-made Files**
Download PENELOPE 2006 material databases from: http://www.oecd-nea.org/science/pubs/2009/nea6416-penelope.pdf

### Step 2: Add the Material to Your Input File

In the `[SECTION MATERIAL FILE LIST]` of your input `.in` file, add a new line:

```
material/mynewmaterial__5-120keV.mcgpu    density=1.050    voxelId=100
```

**Format:**
```
material_file.mcgpu [density=value] [voxelId=comma_separated_list]
```

**Parameters:**
- **`material_file.mcgpu`** - Path to your material properties file (relative or absolute; `.gz` compression supported)
- **`density=value`** (optional) - Density in g/cm³. If omitted, uses the nominal density from the material file
- **`voxelId=list`** (optional) - Comma-separated list of voxel IDs that use this material (e.g., `voxelId=100` or `voxelId=50,75,100`)

### Step 3: Update Your Voxel Geometry

Your voxel geometry file must assign the corresponding voxel ID to locations where you want the new material. The voxel file contains material ID and density for each voxel position.

---

## How to Avoid Overwriting Pre-assigned Materials

### Key Rule

When you include **ANY** `voxelId=` in the input file, the code **abandons ALL hardcoded defaults** and uses only what you explicitly specify.

### Safe Voxel ID Ranges

The following voxel IDs are **NOT used** by default and are safe to use for new materials:

- **3-28** (26 IDs)
- **30-32** (3 IDs)
- **34-39** (6 IDs)
- **41-49** (9 IDs)
- **51-64** (14 IDs)
- **67-87** (21 IDs)
- **89-94** (6 IDs)
- **96-124** (29 IDs)
- **126-149** (24 IDs)
- **151-199** (49 IDs)
- **201-224** (24 IDs)
- **226-249** (24 IDs)
- **251-255** (5 IDs)

### Best Practice

**Use IDs above 200** to avoid any conflict with the original VICTRE phantom voxel values:
- 201, 202, 203, ... up to 224

---

## Examples

### Example 1: Adding a Single New Material (New Mode)

```ini
#[SECTION MATERIAL FILE LIST v.2020-03-03]
material/air__5-120keV.mcgpu.gz                  density=0.0012   voxelId=0
material/adipose__5-120keV.mcgpu.gz              density=0.920    voxelId=1
material/skin__5-120keV.mcgpu.gz                 density=1.090    voxelId=2
material/glandular__5-120keV.mcgpu.gz            density=1.035    voxelId=29
material/muscle__5-120keV.mcgpu.gz               density=1.050    voxelId=40
material/mynewmaterial__5-120keV.mcgpu.gz        density=1.500    voxelId=201  ← New material
```

### Example 2: Same Material with Different Densities

```ini
material/muscle__5-120keV.mcgpu.gz               density=1.050    voxelId=40
material/muscle__5-120keV.mcgpu.gz               density=1.080    voxelId=202  ← Same file, different density
```

### Example 3: Material Assigned to Multiple Voxel IDs

```ini
material/blood__5-120keV.mcgpu.gz                density=1.00     voxelId=150,225  ← Multiple IDs
```

### Example 4: Using Nominal Density from Material File

```ini
material/W__5-120keV.mcgpu                                         voxelId=65   ← No density specified
material/Se__5-120keV.mcgpu                                        voxelId=66   ← Uses nominal from file
```

---

## Important Notes

### Constraints

1. **Maximum Materials**: The code supports up to `MAX_MATERIALS=15` materials (defined in `MC-GPU_v1.5b.h`)
   - This is **NOT a hard limit** - it can be increased by modifying one line in the code
   - See **Increasing Material Limit** section below
2. **Energy Range**: Material files must cover the full energy range of your x-ray spectrum
3. **Voxel ID Range**: Valid values are 0-255 (unsigned char)
4. **Compatibility**: If energy spectrum is outside material file range, the code will report an error

### Memory Considerations

- Material data is stored in GPU constant memory (~64-96 KB limit on most GPUs)
- Rayleigh/Compton interaction tables also consume GPU memory
- **Practical limit: 15-25 materials** (tested and working)
- **Beyond 30 materials**: Risky - may exceed GPU constant memory
- More materials = slightly slower material lookups (negligible performance impact)

### Backward Compatibility

- If you don't specify `density=` and `voxelId=`, the code defaults to VICTRE project values
- Old input files work without modification (100% backward compatible)

---

## Increasing the Material Limit

The 15-material limit is **NOT hard-coded in a critical way** - it can be modified.

### How to Increase MAX_MATERIALS

**Step 1: Edit the header file**

Open `MC-GPU_v1.5b.h` and find line 57:
```cpp
#define  MAX_MATERIALS      15
```

Change it to your desired value (e.g., 20, 25, 30):
```cpp
#define  MAX_MATERIALS      25    // Increase to 25 materials
```

**Step 2: Recompile the code**
```bash
make clean
make
```

### Practical Limits for Different GPU Memory

| Max Materials | Feasibility | Notes |
|---------------|-------------|-------|
| 15 | ✅ **Fully tested** | Current VICTRE setting, working |
| 20-25 | ✅ **Safe** | Should work on modern GPUs |
| 30 | ⚠️ **Risky** | May exceed constant memory on older GPUs |
| 50+ | ❌ **Likely failure** | GPU constant memory typically 64-96 KB |

### Why the Limit Exists

The `MAX_MATERIALS` constant affects several large data structures in GPU constant/global memory:
- **Density lookup table** (small)
- **Rayleigh interaction tables**: `128 * 25005 * MAX_MATERIALS` values (large)
- **Compton interaction tables**: `30 shells * MAX_MATERIALS` values (medium)

Modern NVIDIA GPUs (compute capability 3.5+) usually tolerate 20-25 materials without issues.

### Testing Your Limit

If you increase MAX_MATERIALS and get a CUDA compilation error like:
```
error: "identifier "density_LUT_CONST"" cannot have thread-local storage
```
or similar memory-related errors during runtime, your GPU has exceeded its constant memory limit. Reduce MAX_MATERIALS and recompile.

---

## Material File Format

Material files contain:
- **Material name**
- **Nominal density** (g/cm³)
- **Number of data values** in the interaction cross-section table
- **Mean free paths** for different interactions at various energies:
  - Rayleigh (elastic) scattering
  - Compton (inelastic) scattering
  - Photoelectric absorption
  - Pair production
  - TOTAL mean free path
- **Rayleigh interaction data** (atomic form factors from EPDL database)
- **Compton interaction data** (relativistic impulse model)

These are derived from PENELOPE 2006 database and cover the 5-120 keV energy range typically used in mammography.

---

## Troubleshooting

### Error: "Input voxelID values must be between 0 and 255"

**Solution**: Make sure all voxel IDs are in the range 0-255.

### Error: "Energy spectrum outside material file range"

**Solution**: Ensure your material file covers the full energy range of your x-ray spectrum. The example materials cover 5-120 keV.

### Materials not showing in simulation output

**Possible causes:**
1. Voxel geometry file doesn't contain voxels with the specified IDs
2. Density values are not set correctly
3. Material file path is incorrect or file doesn't exist

### Unexpected density or interaction behavior

**Solution**: Verify that:
- Material file path is correct (relative paths are relative to input file location)
- Density values are in g/cm³
- Voxel IDs in the input file match those in your voxel geometry file

---

## Quick Reference: Migration from Old to New Format

### Old Format (Legacy)
```
#[SECTION MATERIAL FILE LIST v.2009-11-30]
material/air__5-120keV.mcgpu.gz
material/adipose__5-120keV.mcgpu.gz
material/skin__5-120keV.mcgpu.gz
material/glandular__5-120keV.mcgpu.gz
```

### New Format (Explicit)
```
#[SECTION MATERIAL FILE LIST v.2020-03-03]
material/air__5-120keV.mcgpu.gz                  density=0.0012   voxelId=0
material/adipose__5-120keV.mcgpu.gz              density=0.920    voxelId=1
material/skin__5-120keV.mcgpu.gz                 density=1.090    voxelId=2
material/glandular__5-120keV.mcgpu.gz            density=1.035    voxelId=29
```

**Advantage of new format**: Full control over material assignment and density, no hardcoded dependencies.

---

## References

- MC-GPU Official Repository: https://github.com/DIDSR/VICTRE_MCGPU
- PENELOPE 2006 Database: http://www.oecd-nea.org/science/pubs/2009/nea6416-penelope.pdf
- Paper: Andreu Badal et al., "Mammography and breast tomosynthesis simulator for virtual clinical trials," Computer Physics Communications 261 (2021)
