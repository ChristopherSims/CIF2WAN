# CIF2WAN v0.65
CIF2WAN

Edit input.Py then run main.py

See https://arxiv.org/abs/2006.12647 for details

VASP: not well tested

Quantum espresso: well tested

Siesta: testing

Abinit: well tested

Todo:

add GUI

add pure CIF interface

add PP for siesta

add LDA+U for siesta 

Known issues:

wannier90 convergence, wobbly bands. No way to fix, increase num_iter

Abinit SOC case, .mmn not output correctly, error at w90 interface
