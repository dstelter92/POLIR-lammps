# POLIR-lammps
polarizible water potential optimized for IR spectroscopy in LAMMPS

## Please cite:
Mankoo, P. K., & Keyes, T. (2008). POLIR: Polarizable, flexible, transferable
water potential optimized for IR spectroscopy. Journal of Chemical Physics,
129(3). https://doi.org/10.1063/1.2948966

# WARNING!
This code is an unfinished WIP and is not guarenteed to work. Updates will be
put here as they happen.

6/5/18 - Initial Commit, only Bond & Pair potential implemented.

## Install

    git clone https://github.com/dstelter92/POLIR-lammps.git
    cp POLIR-lammps/src/* /path/to/lammps/src/

then re-compile LAMMPS.


## src/

pair_polir_vdw (polir/vdw) is the pairwise O-O potential to attract waters, all other nonbonded
interactions are electrostatic

O-H bond potential uses a 'table' style, see bond/table for details. The
tabulated file is found in the examples directory as "PS_bond_potential.table"


## To-Do list:

- [ ] Implement electrostatics with Thole damping as fix
- [ ] Implement electrostatics with MPI
- [ ] Test against previous work



## Issues:

Currently only works with 1 MPI rank.
