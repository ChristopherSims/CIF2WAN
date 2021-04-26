
# CIF2WAN v0.66
CIF2WAN
![arc_bulk](https://user-images.githubusercontent.com/7741705/112843580-b7d72800-9070-11eb-88be-43e407addfb9.png)

Edit input.in then run CIF2WAN.py

See https://arxiv.org/abs/2006.12647 for details

VASP: testing

Quantum espresso: well tested

Siesta: Well tested (Run on only 1 core, there is an error at wannerization in multicore use)

Abinit: well tested

Todo:

add GUI

Known issues:

LDA+U is not implemented fully in Siesta (beta) will not be added to the program

wannier90 convergence, wobbly bands. No way to fix, increase num_iter

Abinit SOC case, .mmn not output correctly, error at w90 interface
