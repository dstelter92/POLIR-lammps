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
  Compute(lmp, narg, arg),
  list(NULL)
{
  if (narg != 14) error->all(FLERR,"Illegal compute polir/thole/local command");

  local_flag = 1;
  size_local_cols = 0;
  //comm_forward = 3;
  
  typeH = force->inumeric(FLERR,arg[3]);
  typeO = force->inumeric(FLERR,arg[4]);
  alphaH = force->numeric(FLERR,arg[5]);
  alphaO = force->numeric(FLERR,arg[6]);
  
  // Thole damping params
  CD_intra_OH = force->numeric(FLERR,arg[7]);
  CD_intra_HH = force->numeric(FLERR,arg[8]);
  DD_intra_OH = force->numeric(FLERR,arg[9]);
  DD_intra_HH = force->numeric(FLERR,arg[10]);
  CC_inter = force->numeric(FLERR,arg[11]);
  CD_inter = force->numeric(FLERR,arg[12]);
  DD_inter = force->numeric(FLERR,arg[13]);

  nmax = 0;
  if (POLIR_DEBUG)
    fprintf(screen,"DEBUG-mode for compute polir/thole/local is ON\n");
}

/* ---------------------------------------------------------------------- */

ComputePolirTholeLocal::~ComputePolirTholeLocal()
{
  memory->destroy(thole);
}

/* ---------------------------------------------------------------------- */

void ComputePolirTholeLocal::init()
{
  // make multiple computes don't exist
  int count = 0;
  for (int i=0; i<modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"polir/thole/local") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute polir/thole/local");

  // need occasional full neighbor list build
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  memory->create(thole,nmax,"polir/thole/local:thole");
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
    vector_local = thole;
  }

  // invoke neighbor list build/copy as occasional
  neighbor->build_one(list);

  int ii,jj,i,j;
  int inum = list->inum;
  int *ilist = list->ilist;

  // loop over all neighbors, calc thole damping
  for (ii=0; ii<inum; ii++) {
    i = ilist[ii]; // local index1
    for (jj=0; jj<inum; jj++) {
      j = ilist[jj]; // local index2

      // BLAH!

    }
  }


}

/* ---------------------------------------------------------------------- */

void ComputePolirTholeLocal::allocate()
{
  nmax = atom->nmax;
  memory->destroy(thole);
  memory->create(thole,nmax,"polir/thole/local:thole");
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputePolirTholeLocal::memory_usage()
{
  double bytes = nmax*1 * sizeof(double);
  return bytes;
}
