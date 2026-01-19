## DETECTOR
posizione nel codice:

funzione:
```
__device__ inline void tally_image(float* energy, float3* position, float3* direction, signed char* scatter_state, unsigned long long int* image, struct source_struct* source_data_SHARED, struct detector_struct* detector_data_SHARED, int2* seed)     //!!detectorModel!!
```
conto:
```
  dist_detector += -detector_data_SHARED->scintillator_MFP*logf(ranecu(seed));   // Add distance to next interaction inside the detector material to the detector distance    //!!detectorModel!!
```
viene usato un unico valore all'energia media del fascio