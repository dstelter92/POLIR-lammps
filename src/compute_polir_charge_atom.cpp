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
  Compute(lmp, narg-1, arg)
{
  peratom_flag = 1;
  size_peratom_cols = 0;

  // last arg is id of fix/polir
  int len = strlen(arg[narg-1])+1;
  fix_polir = new char[len];
  strcpy(fix_polir,arg[narg-1]);

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&np);
    
  nmax = 0;
  if (POLIR_DEBUG)
    fprintf(screen,"DEBUG-mode for compute polir/charge/atom is ON\n");
}

/* ---------------------------------------------------------------------- */

ComputePolirChargeAtom::~ComputePolirChargeAtom()
{
  memory->destroy(qpolir);
  memory->destroy(qH);
  memory->destroy(qO);
  memory->destroy(roh);
}

/* ---------------------------------------------------------------------- */

void ComputePolirChargeAtom::init()
{
  int count=0;
  for (int i=0; i<modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"polir/charge/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute polir/charge/atom");

  // find polir fix
  int ifix = modify->find_fix(fix_polir);
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

  memory->create(qpolir,nmax,"polir/charge/atom:qpolir");
  memory->create(qH,nmax,"polir/charge/atom:qH");
  memory->create(qO,nmax,"polir/charge/atom:qO");
  memory->create(roh,nmax,2,"polir/charge/atom:roh");
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

  if (!atom->molecular)
    error->all(FLERR,"Atom style must include molecule IDs");

  int i,ii,j,n,nb,j1,j2;
  int partner,igx,igy,igz;
  double delx,dely,delz,rsq,roh_me,roh_partner;
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


  // Loop over all neighbors
  for (i=0; i<nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    nb = num_bond[i];

    // Get ALL local atoms in specified group
    // only check for typeH or typeO, no check is made for 
    // atoms with those types but with different bond types.
    // This was done for flexibility with future polarizible
    // atoms and materials (eg H3O+)

    if (!((type[i] == (*typeO)) || (type[i] == (*typeH)))) {
      error->all(FLERR,"Atom types other than those specified are present in "
          "the current group\n Use a proper 'group' with POLIR");
    }

    if (nb > 2) {
      error->all(FLERR,"More than 2 bonds present, must use this compute with"
          " water molecules");
    }
    
    // calculate bond lengths first
    for (n=0; n<nb; n++) {
      
      // get local index for other atom in bond
      j = atom->map(bond_atom[i][n]);

      // sanity check, should be caught during 'molecular' check above
      if (molid[i] != molid[j]) {
        error->all(FLERR,"Atoms from same bond are not in same molecule,"
            " an accurate LAMMPS data file with moleculeIDs is needed");
      }

      // Catches from compute/bond/local, not sure if needed
      if (j < 0 || !(mask[j] & groupbit)) continue;
      if (newton_bond == 0 && tag[i] > tag[j]) continue;
      //if (bond_type[i][j] == 0) continue;

      // get bond length, roh_me
      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      domain->minimum_image(delx,dely,delz);
      rsq = delx*delx + dely*dely + delz*delz;
      roh_me = sqrt(rsq);

      roh[i][j] = roh_me; // local array of O-H bond lengths

      /*
      partner = -1;
      if (np != 1) {
        if (j > nlocal) // then particle is ghost
          partner = comm->coord2proc(x[j],igx,igy,igz);
        else {
          j1 = atom->map(bond_atom[j][0]);
          partner = comm->coord2proc(x[j1],igx,igy,igx);
        }

        fprintf(screen,"me%d, partner%d",me,partner);

        // if multiproc, send/recv some data
        if (partner > 0) {
          if (me > partner) {
            MPI_Send(&roh_me,1,MPI_DOUBLE,partner,me,world);
            fprintf(screen,"proc%d sent to %d, roh=%g\n",me,partner,roh_me);
          }
          else {
            MPI_Recv(&roh_partner,1,MPI_DOUBLE,partner,me,world,MPI_STATUS_IGNORE);
            fprintf(screen,"proc%d recv from %d, roh=%g\n",me,partner,roh_partner);
          }
        }
      }
      */
    }
  
    if (type[i] == (*typeO)) {

      // use first two bonds to get local atom index
      j1 = atom->map(bond_atom[i][0]);
      j2 = atom->map(bond_atom[i][1]);

      // Calculate charges for all atoms in molecule
      qh1 = (*qeH);
      qh1 += (*c1)*(roh[i][j1] + roh[i][j2] - (2*(*reOH)));
      qh1 += (*c3)*(roh[i][j1] - roh[i][j2]);

      qh2 = (*qeH);
      qh2 += (*c1)*(roh[i][j1] + roh[i][j2] - (2*(*reOH)));
      qh2 += (*c3)*(roh[i][j2] - roh[i][j1]);

      q = -(qh1 + qh2); // make negative later

      // Save to local charge array
      qO[i] = q;
      qH[j1] = qh1;
      qH[j2] = qh2;

    
      if (POLIR_DEBUG) {
        fprintf(screen,"molid:%d proc:%d\n  roh1[%d][%d]=%f roh2[%d][%d]=%f\n  "
            "q[%d]=%g qh1[%d]=%g qh2[%d]=%g\n",
            molid[i],me,i,j1,roh[i][j1],i,j2,roh[i][j2],i,q,j1,qh1,j2,qh2);
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
  memory->destroy(roh);
  memory->create(qpolir,nmax,"polir/charge/atom:qpolir");
  memory->create(qH,nmax,"polir/charge/atom:qH");
  memory->create(qO,nmax,"polir/charge/atom:qO");
  memory->create(roh,nmax,2,"polir/charge/atom:roh");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePolirChargeAtom::memory_usage()
{
  double bytes = nmax*5 * sizeof(double);
  return bytes;
}

