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

#include <cstring>
#include <cmath>
#include "compute_polir_charge_atom.h"
#include "atom.h"
#include "neighbor.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePolirChargeAtom::ComputePolirChargeAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 9) error->all(FLERR,"Illegal compute polir/charge/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  typeH = force->inumeric(FLERR,arg[3]);
  typeO = force->inumeric(FLERR,arg[4]);
  qeH = force->numeric(FLERR,arg[5]);
  reOH = force->numeric(FLERR,arg[6]);
  c1 = force->numeric(FLERR,arg[7]);
  c3 = force->numeric(FLERR,arg[8]);

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputePolirChargeAtom::~ComputePolirChargeAtom()
{
  memory->destroy(qpolir);
}

/* ---------------------------------------------------------------------- */

void ComputePolirChargeAtom::init()
{
  int count = 0;
  for (int i=0; i<modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"polir/charge/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute polir/charge/atom");

  memory->create(qpolir,nmax,"polir/charge/atom:qpolir");
}

/* ---------------------------------------------------------------------- */

void ComputePolirChargeAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow local qpolir array if necessary
  // needs to be atom->nmax in length
  if (atom->nmax > nmax) {
    allocate();
    vector_atom = qpolir;
  }

  int i,n,i1,i2,j1,j2;
  int mol,last_molid;
  double delix,deliy,deliz,deljx,deljy,deljz;
  double rsq,roh,roh_other;
  double qh,qo;

  double **x = atom->x;
  int *type = atom->type;
  tagint *molid = atom->molecule;

  int nbondlist = neighbor->nbondlist;
  int **bondlist = neighbor->bondlist;

  // clear local qpolir array
  for (i=0; i<nmax; i++) qpolir[i] = 0.0;

  for (n=0; n<nbondlist; n++) {
    // indicies for this bond 
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    // indicies for partner bond (in same H2O molec)
    if (n == 0) {
      j1 = bondlist[n+1][0];
      j2 = bondlist[n+1][1];
    }
    else {
      j1 = bondlist[n-1][0];
      j2 = bondlist[n-1][1];
   }

    mol = molid[i1]; // set this molecule id

    // sanity check on bonds
    if ((molid[i1] != molid[i2]) || (molid[j1] != molid[j2])) {
      error->all(FLERR,"Atoms within bond do not belong to same molid,"
          "accurate LAMMPS data file is needed");
    }

    // Check other neighboring bonds if wrong on 1st guess
    if ((molid[i1] != molid[j1]) && (n>0)) {
      j1 = bondlist[n+1][0];
      j2 = bondlist[n+1][1];
    }

    // Kill if cannot find bond in same molid
    if (molid[i1] != molid[j1])
      error->all(FLERR,"Unable to identify bondlength for charge calculation");

    // get bond lengths, roh and roh_other
    delix = x[i1][0] - x[i2][0];
    deliy = x[i1][1] - x[i2][1];
    deliz = x[i1][2] - x[i2][2];
    rsq = delix*delix + deliy*deliy + deliz*deliz;
    roh = sqrt(rsq);

    deljx = x[j1][0] - x[j2][0];
    deljy = x[j1][1] - x[j2][1];
    deljz = x[j1][2] - x[j2][2];
    rsq = deljx*deljx + deljy*deljy + deljz*deljz;
    roh_other = sqrt(rsq);
    //fprintf(screen,"molid:%d MyOH:%g, otherOH:%g\n",molid[i1],roh,roh_other);

    // reset charges
    qh = 0.0;
    if (mol != last_molid)
      qo = 0.0;

    // Check atom types, calc charges based on bond lengths
    // Save to qpolir array
    if (type[i1] == typeH) {
      qh = qeH;
      qh += c1*(roh + roh_other - (2*reOH));
      qh += c3*(roh - roh_other);
      qo += -qh;

      qpolir[i1] = qh;
      qpolir[i2] = qo;
    }
    else if (type[i2] == typeH) {
      qh = qeH;
      qh += c1*(roh + roh_other - (2*reOH));
      qh += c3*(roh - roh_other);
      qo += -qh;

      qpolir[i2] = qh;
      qpolir[i1] = qo;
    }
    else
      error->all(FLERR,"No input atom types match existing types");

    last_molid = mol; // set current molid as last_molid
  }
}

/* ---------------------------------------------------------------------- */

void ComputePolirChargeAtom::allocate()
{
  nmax = atom->nmax;
  memory->destroy(qpolir);
  memory->create(qpolir,nmax,"polir/charge/atom:qpolir");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePolirChargeAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}

