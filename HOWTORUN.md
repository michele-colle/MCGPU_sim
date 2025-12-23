per utilizzare il sw dopo aver compilato con make dalla cartella progetto:
```
make .
cp MC-GPU_v1.5b example_simulations/ #credo funzioni
```
va lanciata la simulazione
```
cd example_simulations
./MC-GPU_v1.5b.x cbct_circular.in | tee cbct_circular.out
```