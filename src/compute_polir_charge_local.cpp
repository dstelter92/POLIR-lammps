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
#include "compute_polir_charge_local.h"
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
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DEBUG 0

/* ---------------------------------------------------------------------- */

ComputePolirChargeLocal::ComputePolirChargeLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 9) error->all(FLERR,"Illegal compute polir/charge/local command");

  local_flag = 1;
  size_local_cols = 0;

  typeH = force->inumeric(FLERR,arg[3]);
  typeO = force->inumeric(FLERR,arg[4]);
  qeH = force->numeric(FLERR,arg[5]);
  reOH = force->numeric(FLERR,arg[6]);
  c1 = force->numeric(FLERR,arg[7]);
  c3 = force->numeric(FLERR,arg[8]);

  nlocal = 0;
  if (DEBUG)
    fprintf(screen,"DEBUG-mode for compute polir/charge/local is ON\n");
}

/* ---------------------------------------------------------------------- */

ComputePolirChargeLocal::~ComputePolirChargeLocal()
{
  memory->destroy(qpolir);
}

/* ---------------------------------------------------------------------- */

void ComputePolirChargeLocal::init()
{
  int count = 0;
  for (int i=0; i<modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"polir/charge/local") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute polir/charge/local");

  memory->create(qpolir,nlocal,"polir/charge/local:qpolir");
}

/* ---------------------------------------------------------------------- */

void ComputePolirChargeLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // grow local qpolir array if necessary
  if (atom->nlocal > nlocal) {
    allocate();
    vector_local = qpolir;
  }

  if (!atom->molecular)
    error->all(FLERR,"Atom style must include molecule IDs");

  int i,i1,i2;
  int nb,mol,last_molid;
  double delx,dely,delz,rsq,roh1,roh2;
  double qh1,qh2,qo;

  double **x = atom->x;
  int *type = atom->type;
  int *num_bond = atom->num_bond;
  int *mask = atom->mask;
  int newton_bond = force->newton_bond;
  tagint **bond_atom = atom->bond_atom;
  tagint *molid = atom->molecule;
  tagint *tag = atom->tag;


  // clear local qpolir array
  for (i=0; i<nlocal; i++) qpolir[i] = 0.0;

  for (i=0; i<nlocal; i++) {
    if (mask[i] & groupbit) {

      // Get local atoms that are O and have 2 bonds, skip others
      nb = num_bond[i];
      if ((nb != 2) && (type[i] != typeO))
        continue;
    
      // atom indices (H atoms) bonded to n
      i1 = atom->map(bond_atom[i][0]);
      i2 = atom->map(bond_atom[i][1]);
      if (DEBUG)
        fprintf(screen,"atoms: O:%d H1:%d H2:%d\n",i,i1,i2);

      if ((type[i1] != typeH) || (type[i2] != typeH))
        error->all(FLERR,"Atom types do not match those specified in compute");

      mol = molid[i];
      if ((molid[i1] != molid[i]) || (molid[i2] != molid[i])) {
        error->all(FLERR,"Atoms within same bond do not belong to the same molid,"
            "an accurate LAMMPS data file with moleculeIDs is needed");
      }

      // Catches from compute/bond/local, not sure if needed
      if (i1 < 0 || !(mask[i1] & groupbit)) continue;
      if (i2 < 0 || !(mask[i2] & groupbit)) continue;
      if (newton_bond == 0 && tag[i] > tag[i1]) continue;
      if (newton_bond == 0 && tag[i] > tag[i2]) continue;

      // get bond lengths, roh and roh_other
      delx = x[i][0] - x[i1][0];
      dely = x[i][1] - x[i1][1];
      delz = x[i][2] - x[i1][2];
      domain->minimum_image(delx,dely,delz);
      rsq = delx*delx + dely*dely + delz*delz;
      roh1 = sqrt(rsq);

      delx = x[i][0] - x[i2][0];
      dely = x[i][1] - x[i2][1];
      delz = x[i][2] - x[i2][2];
      domain->minimum_image(delx,dely,delz);
      rsq = delx*delx + dely*dely + delz*delz;
      roh2 = sqrt(rsq);

      if (DEBUG)
        fprintf(screen,"molid:%d rOH1=%g, rOH2=%g\n",molid[i],roh1,roh2);

      // calc charges based on bond lengths
      qh1 = qeH;
      qh1 += c1*(roh1 + roh2 - (2*reOH));
      qh1 += c3*(roh1 - roh2);

      qh2 = qeH;
      qh2 += c1*(roh1 + roh2 - (2*reOH));
      qh2 += c3*(roh2 - roh1);

      qo = -(qh1 + qh2);
    
      // Save to qpolir array
      qpolir[i] = qo;
      qpolir[i1] = qh1;
      qpolir[i2] = qh2;
      if (DEBUG)
        fprintf(screen,"qo=%g, qh1=%g qh2=%g\n",qo,qh1,qh2);

    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePolirChargeLocal::allocate()
{
  nlocal = atom->nlocal;
  memory->destroy(qpolir);
  memory->create(qpolir,nlocal,"polir/charge/local:qpolir");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePolirChargeLocal::memory_usage()
{
  double bytes = nlocal * sizeof(double);
  return bytes;
}

