/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: David Stelter (Boston University)
------------------------------------------------------------------------- */
#include <cstring>
#include <cmath>
#include "stdlib.h"
#include "compute_polir_charge_atom.h"
#include "fix_polir.h"
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define POLIR_DEBUG 0

/* ---------------------------------------------------------------------- */

ComputePolirChargeAtom::ComputePolirChargeAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 5)
    error->all(FLERR,"Illegal compute POLIR/CHARGE/ATOM command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  // get parameters from input
  int len = strlen(arg[3])+1;
  fix_polir = new char[len];
  strcpy(fix_polir,arg[3]);

  len = strlen(arg[4])+1;
  compute_lbond = new char[len];
  strcpy(compute_lbond,arg[4]);

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&np);
    
  nmax = 0;
  if (POLIR_DEBUG)
    fprintf(screen,"DEBUG-mode for compute POLIR/CHARGE/ATOM is ON\n");
}

/* ---------------------------------------------------------------------- */

ComputePolirChargeAtom::~ComputePolirChargeAtom()
{
  memory->destroy(qpolir);
  memory->destroy(qH);
  memory->destroy(qO);
}

/* ---------------------------------------------------------------------- */

void ComputePolirChargeAtom::init()
{
  // check for duplicate computes
  int count=0;
  for (int i=0; i<modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"POLIR/CHARGE/ATOM") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute POLIR/CHARGE/ATOM");

  // find polir fix
  ifix = modify->find_fix(fix_polir);
  if (ifix < 0)
    error->all(FLERR,"Fix polir ID for compute POLIR/CHARGE/ATOM does not exist");

  int dim;
  typeH = (int *)modify->fix[ifix]->extract("typeH",dim);
  typeO = (int *)modify->fix[ifix]->extract("typeO",dim);
  c1 = (double *)modify->fix[ifix]->extract("c1",dim);
  c3 = (double *)modify->fix[ifix]->extract("c3",dim);
  qeH = (double *)modify->fix[ifix]->extract("qeH",dim);
  reOH = (double *)modify->fix[ifix]->extract("reOH",dim);
    
  if (dim != 0)
    error->all(FLERR,"Cannot extract fix/polir inputs and constants");

  // find bond/local compute
  icompute = modify->find_compute(compute_lbond);
  if (icompute < 0)
    error->all(FLERR,"Compute lbond ID for compute bond/local does not exist");

  // memory init
  memory->create(qpolir,nmax,"POLIR/CHARGE/ATOM:qpolir");
  memory->create(qH,nmax,"POLIR/CHARGE/ATOM:qH");
  memory->create(qO,nmax,"POLIR/CHARGE/ATOM:qO");
}

/* ---------------------------------------------------------------------- */

void ComputePolirChargeAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow local qpolir array if necessary
  if (atom->nmax > nmax) {
    allocate();
    vector_atom = qpolir;
  }

  // system must of molecule IDs in LAMMPS.data file
  if (!atom->molecular)
    error->all(FLERR,"Atom style must include molecule IDs");

  /*
  int i,ii,j,n,nb,j1,j2;
  int partner,igx,igy,igz;
  double delx,dely,delz,rsq,roh_me,roh_partner;
  double qh1,qh2,q;
  */
  int i,j,k,m,nb,j1,j2;
  double qh1,qh2,q;

  double **x = atom->x;
  int *type = atom->type;
  int *num_bond = atom->num_bond;
  int *mask = atom->mask;
  int newton_bond = force->newton_bond;
  int nlocal = atom->nlocal;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  tagint *molid = atom->molecule;
  tagint *tag = atom->tag;

  // clear local array
  for (i=0; i<nmax; i++)
    qpolir[i] = 0.0;

  // get info from bond/local, let fix/polir run calculations
  // only access arrays here
  //modify->compute[icompute]->compute_local();
  lbond = modify->compute[icompute]->vector_local;
  nbonds = modify->compute[icompute]->size_local_rows;

  
  // loop over local atoms in this proc
  m = 0;
  for (i=0; i<nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    nb = num_bond[i];
    
    // Get ALL local atoms in specified group
    // only check for typeH or typeO, no check is made for 
    // atoms with those types or with different bond types.
    // This was done for flexibility with future polarizible
    // atoms and materials (eg H3O+)
    
    if (!((type[i] == (*typeO)) || (type[i] == (*typeH)))) {
      error->all(FLERR,"Atom types other than those specified are present in "
          "the current group\n Use a proper 'group' with POLIR");
    }

    for (k=0; k<nb; k++) {
      j = atom->map(bond_atom[i][k]); 

      // sanity check, should be caught during 'molecular check above
      if (molid[i] != molid[j]) {
        error->all(FLERR,"Atoms from same bond are not in same molecule,"
            " an accurate LAMMPS data file with moleculeIDs is needed");
      }

      if (POLIR_DEBUG) {
        double delx = x[i][0] - x[j][0];
        double dely = x[i][1] - x[j][1];
        double delz = x[i][2] - x[j][2];
        domain->minimum_image(delx,dely,delz);
        double rsq = delx*delx + dely*dely + delz*delz;
        double r = sqrt(rsq);

        fprintf(screen,"proc%d lbond[%d]=%g (%d-%d), r=%g\n",me,m,lbond[m],i,j,r);
      }
      m++;
    }

    if (type[i] == (*typeO)) {

      // 2nd local index in first 2 bonds
      j1 = atom->map(bond_atom[i][0]);
      j2 = atom->map(bond_atom[i][1]);

      // calculate charges for all atoms in molecule
      qh1 = (*qeH);
      qh1 += (*c1)*(lbond[m-1] + lbond[m-2] - (2*(*reOH)));
      qh1 += (*c3)*(lbond[m-1] - lbond[m-2]);

      qh2 = (*qeH);
      qh2 += (*c1)*(lbond[m-1] + lbond[m-2] - (2*(*reOH)));
      qh2 += (*c3)*(lbond[m-2] - lbond[m-1]);

      q = -(qh1 + qh2); // make negative later

      // Save to local charge array
      qO[i] = q;
      qH[j1] = qh1;
      qH[j2] = qh2;

      if (POLIR_DEBUG) {
        fprintf(screen,"molid:%d proc:%d\n  lbond1[%d][%d]=%f lbond2[%d][%d]=%f\n  "
            "q[%d]=%g qh1[%d]=%g qh2[%d]=%g\n",
            molid[i],me,i,j1,lbond[m],i,j2,lbond[m-1],i,q,j1,qh1,j2,qh2);
      }
    }
  }

  // Make sure all procs know charges on H
  //MPI_Allreduce(MPI_IN_PLACE,qO,nmax,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(MPI_IN_PLACE,qH,nmax,MPI_DOUBLE,MPI_MAX,world);

  // construct peratom charge array
  for (i=0; i<nmax; i++) {
    if (type[i] == (*typeO))
      qpolir[i] = qO[i];
    else if (type[i] == (*typeH))
      qpolir[i] = qH[i];
  }


}

/* ---------------------------------------------------------------------- */

void ComputePolirChargeAtom::allocate()
{
  nmax = atom->nmax;
  memory->destroy(qpolir);
  memory->destroy(qH);
  memory->destroy(qO);
  memory->create(qpolir,nmax,"POLIR/CHARGE/ATOM:qpolir");
  memory->create(qH,nmax,"POLIR/CHARGE/ATOM:qH");
  memory->create(qO,nmax,"POLIR/CHARGE/ATOM:qO");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePolirChargeAtom::memory_usage()
{
  double bytes = nmax*5 * sizeof(double);
  return bytes;
}

