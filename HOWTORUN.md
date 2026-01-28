per utilizzare il sw dopo aver compilato con make dalla cartella progetto:
```
make
cp MC-GPU_v1.5b example_simulations/ #credo funzioni
```
va lanciata la simulazione
```
cd example_simulations
./MC-GPU_v1.5b.x cbct_circular.in | tee cbct_circular.out
```

CREARE MATERIALI CON PENELOPE da sorgente
```
#compilo material.f di penelope e lo sposto in pendbase
cd ~/penelope/fsource
gfortran -O3 material.f -o material.x
mv material.x ../pendbase/

#compilo MC-GPU_create_material_data.f e sposto anche lui in pendbase (nota, il sorgente é stato adattato alla nuova versione di penelope.. dalla AI quindi paura vera)
gfortran -O3 MC-GPU_create_material_data.f penelope.f rita.f -o create_material.x
mv create_material.x ../pendbase/

cd ../pendbase/
./material.x
```
a questo punto il software chiede che materiale generare, la densitá e altra roba.
elenco materiali in pendbase/pdfiles/pdcompos.pen
poi crea un file, conviene dargli estensione .mat

materiali creati da me:
al.mat
pmma.mat 224  POLYMETHYL METHACRALATE (LUCITE, PERSPEX, PLEXIGLASS)  (224)
water.mat 278  WATER, LIQUID  --ICRU90  (278)
adipose_icrp.mat 103  ADIPOSE TISSUE (ICRP)
blood_icrp.mat 118  BLOOD (ICRP)  (118)
bone_cortical_icrp.mat 120  BONE, CORTICAL (ICRP)  (120) 
brain_icrp.mat 123  BRAIN (ICRP)  (123)
eye_lens_icrp.mat 156  EYE LENS (ICRP)  (156)
lung_icrp.mat 191  LUNG (ICRP)  (191)
muscle_skeletal_icrp.mat 202  MUSCLE, SKELETAL (ICRP)  (202)  
skin_icrp.mat 251  SKIN (ICRP)  (251)
testes_icrp.mat 259  TESTES (ICRP)  (259) 
tissue_soft_icrp.mat 262  TISSUE, SOFT (ICRP)  (262)

poi per creare effettivamente il file per mcgpu usa
`./create_material.x`
range: 
5000 120010
nbin
23002