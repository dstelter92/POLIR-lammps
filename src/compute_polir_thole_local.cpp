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

#include <cmath>
#include <cstring>
#include "compute_polir_thole_local.h"
#include "fix_polir.h"
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define POLIR_DEBUG 1

/* ---------------------------------------------------------------------- */

ComputePolirTholeLocal::ComputePolirTholeLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg-1, arg),
  list(NULL)
{
  // default calculate all values (7) of thole potential for all local pairs
  nvalues = 7; 

  local_flag = 1;
  size_local_cols = nvalues;
  
  // last arg is id of fix/polir
  int len = strlen(arg[narg-1])+1;
  fix_polir = new char[len];
  strcpy(fix_polir,arg[narg-1]);
  
  nmax = 0;
  if (POLIR_DEBUG)
    fprintf(screen,"DEBUG-mode for compute POLIR/THOLE/LOCAL is ON\n");
}

/* ---------------------------------------------------------------------- */

ComputePolirTholeLocal::~ComputePolirTholeLocal()
{
  memory->destroy(thole);
  memory->destroy(damping);
}

/* ---------------------------------------------------------------------- */

void ComputePolirTholeLocal::init()
{
  // make multiple computes don't exist
  int count = 0;
  for (int i=0; i<modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"POLIR/THOLE/LOCAL") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute POLIR/THOLE/LOCAL");

  // need occasional full neighbor list build
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
  
  memory->create(thole,nmax,nvalues,"POLIR/THOLE/LOCAL:thole");
  memory->create(damping,nvalues,"POLIR/THOLE/LOCAL:damping");

  // find polir fix
  int ifix = modify->find_fix(fix_polir);
  if (ifix < 0)
    error->all(FLERR,"Fix polir ID for compute POLIR/THOLE/LOCAL does not exist");

  int dim;
  // get values needed here in compute/polir/thole
  typeH = (int *)modify->fix[ifix]->extract("typeH",dim);
  typeO = (int *)modify->fix[ifix]->extract("typeO",dim);
  CD_intra_OH = (double *)modify->fix[ifix]->extract("CD_intra_OH",dim);
  CD_intra_HH = (double *)modify->fix[ifix]->extract("CD_intra_HH",dim);
  DD_intra_OH = (double *)modify->fix[ifix]->extract("DD_intra_OH",dim);
  DD_intra_HH = (double *)modify->fix[ifix]->extract("DD_intra_HH",dim);
  CC_inter = (double *)modify->fix[ifix]->extract("CC_inter",dim);
  CD_inter = (double *)modify->fix[ifix]->extract("CD_inter",dim);
  DD_inter = (double *)modify->fix[ifix]->extract("DD_inter",dim);

  // pack constants into array
  for (int i=0; i<nvalues; i++) {
    if (i=0)
      damping[i] = (*CD_intra_OH);
    else if (i=1)
      damping[i] = (*CD_intra_HH);
    else if (i=2)
      damping[i] = (*DD_intra_OH);
    else if (i=3)
      damping[i] = (*DD_intra_HH);
    else if (i=4)
      damping[i] = (*CC_inter);
    else if (i=5)
      damping[i] = (*CD_inter);
    else if (i=6)
      damping[i] = (*DD_inter);
    else
      error->all(FLERR,"Too many damping coeffs");
  }

  if (dim != 0)
    error->all(FLERR,"Cannot extract fix/polir inputs and constants");

}

/* ---------------------------------------------------------------------- */

void ComputePolirTholeLocal::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePolirTholeLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // check that nmax is same, otherwise update arrays
  if (atom->nmax > nmax) {
    allocate();
    //array_local = thole;
  }

  // invoke neighbor list build/copy as occasional
  neighbor->build_one(list);

  int ii,jj,m,i,j,k;
  double alphai,alphaj;
  double delx,dely,delz,rsq,r;
  double t1,t2,ra4,ra;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *mask = atom->mask;
  int *type = atom->type;
  double **x = atom->x;


  // loop over all neighbors, calc thole damping
  m=0;
  for (ii=0; ii<inum; ii++) {
    if (!(mask[i] & groupbit)) continue;
    i = ilist[ii]; // local index1

    if (type[i] == (*typeH))
      alphai = alphaH;
    else if (type[i] == (*typeO))
      alphai = alphaO;
    else 
      error->all(FLERR,"No polarizibility to assign to atom type");

    for (jj=0; jj<inum; jj++) {
      if (!(mask[j] & groupbit)) continue;
      j = ilist[jj]; // local index2

      if (type[j] == (*typeH))
        alphaj = alphaH;
      else if (type[i] == (*typeO))
        alphaj = alphaO;
      else
        error->all(FLERR,"No polarizibility to assign to atom type");

      // calc polarizbility constant
      A_const = pow(alphai*alphaj,1/6);

      // get distance between particles
      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);

      // other constants
      ra = (r / A_const);
      ra4 = pow(r,4);

      for (k=0; k<nvalues; k++) {
        t1 = exp(-damping[k]*ra4);
        t2 = pow(damping[k],1/4)*ra*0; // for now no gamma fctn implemented
        
        thole[m][k] = (1 - t1 + t2) / r;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePolirTholeLocal::allocate()
{
  nmax = atom->nmax;
  memory->destroy(thole);
  memory->create(thole,nmax,nvalues,"POLIR/THOLE/LOCAL:thole");
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputePolirTholeLocal::memory_usage()
{
  double bytes = nmax*2 * sizeof(double);
  return bytes;
}
