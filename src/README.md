# Please cite:
Mankoo, P. K., & Keyes, T. (2008). POLIR: Polarizable, flexible, transferable
water potential optimized for IR spectroscopy. Journal of Chemical Physics,
129(3). https://doi.org/10.1063/1.2948966

# POLIR

pair_polir is the pairwise O-O potential to attract waters, all other nonbonded
interactions are electrostatic.

O-H bond potential uses a 'table' style, see bond/table for details. The
tabulated file is found in this directory as "PS_bond_potential.table"


## To-Do list:

[ ] Implement electrostatics with Thole damping as fix
[ ] Implement electrostatics with MPI
[ ] Test against previous work


## Known Issues:

None, but code is unfinished.
