# CIF2WAN v0.66
CIF2WAN

Edit input.Py then run main.py

See https://arxiv.org/abs/2006.12647 for details

VASP: testing

Quantum espresso: well tested

Siesta: Well tested

Abinit: well tested

Todo:

add GUI

add pure CIF interface

add PP for siesta



Known issues:

LDA+U is not implemented fully in Siesta (beta) will not be added to the program

wannier90 convergence, wobbly bands. No way to fix, increase num_iter

Abinit SOC case, .mmn not output correctly, error at w90 interface
