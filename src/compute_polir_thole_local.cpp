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
#include <boost/math/special_functions/gamma.hpp>
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

#define POLIR_DEBUG 0

/* ---------------------------------------------------------------------- */

ComputePolirTholeLocal::ComputePolirTholeLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg-1, arg),
  list(NULL)
{
  // default calculate all values (30) of thole potential for all local pairs
  nvalues = 2 + (7*4);
  ndamp = 7;

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
  // make sure multiple computes don't exist
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
  memory->create(damping,7,"POLIR/THOLE/LOCAL:damping");

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

  gam = boost::math::tgamma(0.75);

  // put constants into array
  damping[0] = *CD_intra_OH;
  damping[1] = *CD_intra_HH;
  damping[2] = *DD_intra_OH;
  damping[3] = *DD_intra_HH;
  damping[4] = *CC_inter;
  damping[5] = *CD_inter;
  damping[6] = *DD_inter;

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

  npairs = compute_pairs(0);
  if (npairs > nmax) {
    allocate();
    array_local = thole;
  }

  size_local_rows = npairs;
  compute_pairs(1);

}

/* ---------------------------------------------------------------------- */

int ComputePolirTholeLocal::compute_pairs(int flag)
{
  // invoke neighbor list build/copy as occasional
  if (flag == 0)
    neighbor->build_one(list);

  int ii,jj,m,i,j,k,n;
  double alphai,alphaj;
  double delx,dely,delz,rsq,r;
  double ra,ra4,rr,rr2,rr3,rr4,expfac,a14,igam;

  inum = list->inum;
  ilist = list->ilist;
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
    else {
      if (POLIR_DEBUG)
        fprintf(screen,"atom%d has type%d\n",i,type[i]);
      error->all(FLERR,"No polarizibility to assign to atom type");
    }

    for (jj=0; jj<inum; jj++) {
      j = ilist[jj]; // local index 2
      if (i == j) continue;
      if (!(mask[j] & groupbit)) continue;

      if (type[j] == (*typeH))
        alphaj = alphaH;
      else if (type[j] == (*typeO))
        alphaj = alphaO;
      else {
        if (POLIR_DEBUG)
          fprintf(screen,"atom%d has type%d\n",j,type[j]);
        error->all(FLERR,"No polarizibility to assign to atom type");
      }

      if (flag > 0) { // only count if flag == 0
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
        ra4 = pow(ra,4);
        rr = 1 / r;
        rr2 = rr * rr;
        rr3 = rr2 * rr;
        rr4 = rr2 * rr2;

        if (POLIR_DEBUG)
          fprintf(screen,"thole for pair%d (%d-%d) r=%g:\n  ",m,i,j,r);

        thole[m][0] = i;
        thole[m][1] = j;
        for (k=0; k<ndamp; k++) {
          n = 4*k + 2;
          // compute all derivs for all types, compute extra derivs to determine
          // forces later
          expfac = exp(-damping[k]*ra4);
          a14 = pow(damping[k],1/4);
          igam = boost::math::gamma_q(0.75,damping[k]*ra4);

          // 0th derive
          thole[m][n] = (1.0 - expfac);
          thole[m][n] += a14*ra*gam*igam;
          // 1st, ignore gamma fctn term, paper does this but IDK why as of
          // 6/20/18. It's also like this in the old fortran code...
          // keep it for now..
          thole[m][n+1] = rr2 * (-1*(1.0 - expfac));
          // 2nd
          thole[m][n+2] = rr3 * ((-4*ra4*damping[k]*expfac) + (2*(1.0 - expfac)));
          // 3rd
          thole[m][n+3] = rr4 * (-6*(1.0 - expfac) + (-4*ra4*damping[k]*expfac)
              + (16*damping[k]*damping[k]*ra4*ra4*expfac));


          if (POLIR_DEBUG)
            fprintf(screen,"t[%d]=%g  ",k,thole[m][k]);
        }

        if (POLIR_DEBUG)
          fprintf(screen,"\n");
      }
      m++;
    }
  }
  return m;
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
  double bytes = nmax*nvalues * sizeof(double);
  bytes += (nvalues-2) * sizeof(double);
  return bytes;
}
