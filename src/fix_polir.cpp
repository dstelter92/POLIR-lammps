/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Force scaling fix for gREM.
   Cite: http://dx.doi.org/10.1063/1.3432176
   Cite: http://dx.doi.org/10.1021/acs.jpcb.5b07614

------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: David Stelter (Boston University)
------------------------------------------------------------------------- */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include "comm.h"
#include "group.h"
#include "fix_polir.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "input.h"
#include "compute.h"
#include "memory.h"
#include "neighbor.h"
#include "error.h"
#include "universe.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

#define INVOKED_LOCAL 1
#define POLIR_DEBUG 1

/* ---------------------------------------------------------------------- */

FixPolir::FixPolir(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 11) error->all(FLERR,"Illegal fix polir command");

  local_flag = 1;

  typeH = force->inumeric(FLERR,arg[3]);  // atom type for H
  typeO = force->inumeric(FLERR,arg[4]);  // atom type for O
  qeH = force->numeric(FLERR,arg[5]);     // equilibrium charge on H
  reOH = force->numeric(FLERR,arg[6]);    // equilibrium O-H bond length
  uH = force->numeric(FLERR,arg[7]);      // H dipole along HOH bisector
  uO = force->numeric(FLERR,arg[8]);      // O dipole along OH bond
  alphaH = force->numeric(FLERR,arg[9]);  // polarizibility for typeH
  alphaO = force->numeric(FLERR,arg[10]); // polarizibility for type0


  // Set output defaults
  polir_output = 0;
  if ((comm->me == 0) && (logfile)) 
    polir_output = 1;

  if (POLIR_DEBUG) {
    if (comm->me == 0)
      fprintf(screen,"DEBUG-mode for fix polir is ON\n");
  }
}

/* ---------------------------------------------------------------------- */

FixPolir::~FixPolir()
{
  memory->destroy(charges);
}

/* ---------------------------------------------------------------------- */

int FixPolir::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPolir::init()
{
  // Init stuff here...
  int i;
  int nmax = atom->nmax;

  // Constants... 
  // see fix_polir.h for details

  // bond length coeff
  c1 = 0.120;
  c3 = -0.136;

  // Thole damping coeff
  CD_intra_OH = 0.650;
  CD_intra_HH = 0.050;
  DD_intra_OH = 0.690;
  DD_intra_HH = 0.050;
  CC_inter = 0.50;
  CD_inter = 0.15;
  DD_inter = 0.30;

  // MPI init
  me_universe = universe->me;
  MPI_Comm_rank(world,&me);
  nworlds = universe->nworlds;
  iworld = universe->iworld;
  
  // memory management
  memory->create(charges,nmax,"polir:charges");

  // Search for charge compute
  q_compute_id = -1;
  for (i=0; i<modify->ncompute; i++) {
    if (strcmp(modify->compute[i]->style,"polir/charge/atom") == 0) {
      q_compute_id = i;
      break;
    }
  }
  // If doesn't exist, create a new one
  if (q_compute_id < 0) {
    int n = strlen(id) + 4;
    id_q = new char[n];
    strcpy(id_q,id);
    strcat(id_q,"_polir_charge_atom");
    char param1[50];
    snprintf(param1,50,"%d",typeH);
    char param2[50];
    snprintf(param2,50,"%d",typeO);
    char param3[50];
    snprintf(param3,50,"%f",qeH);
    char param4[50];
    snprintf(param4,50,"%f",reOH);
    char param5[50];
    snprintf(param5,50,"%f",c1);
    char param6[50];
    snprintf(param6,50,"%f",c3);

    char **newarg = new char*[9];
    newarg[0] = id_q;
    newarg[1] = group->names[igroup];
    newarg[2] = (char *) "polir/charge/atom";
    newarg[3] = param1;
    newarg[4] = param2;
    newarg[5] = param3;
    newarg[6] = param4;
    newarg[7] = param5;
    newarg[8] = param6;

    modify->add_compute(9,newarg);
    delete [] newarg;

    q_compute_id = modify->ncompute -1;
  }
  compute_pca = modify->compute[q_compute_id];



}

/* ---------------------------------------------------------------------- */

void FixPolir::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) {
    post_force(vflag);

    // Invoke compute to run on each step
    compute_pca->invoked_flag |= INVOKED_LOCAL;
    modify->addstep_compute(update->ntimestep + 1);
  }
  else
    error->all(FLERR,"Currently expecting run_style verlet");
}

/* ---------------------------------------------------------------------- */

void FixPolir::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPolir::pre_force(int vflag)
{
  int i;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nmax = atom->nmax;

  // Re-calculate atomic charges, invoke per-atom fix
  compute_pca->compute_peratom();
  charges = compute_pca->vector_atom;
  /*
   MPI_Allgather(charges,nmax,MPI_DOUBLE,global_vector,nmax*nworlds,
      MPI_DOUBLE,world);
  */

  //fprintf(screen,"nlocal=%d\n",nlocal);
  for (i=0; i<nlocal; i++) {
    if (mask[i] & groupbit) {
      fprintf(screen,"proc:%d atom%d:%g\n",me,tag[i],charges[i]);
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixPolir::post_force(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

}

/* ---------------------------------------------------------------------- */

void FixPolir::end_of_step()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

}

/* ---------------------------------------------------------------------- */

void FixPolir::allocate()
{
  int nmax = atom->nmax;
  memory->destroy(charges);
  memory->create(charges,nmax,"polir:charges");
}

/* ---------------------------------------------------------------------- */

double FixPolir::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax*1 * sizeof(double);
  bytes += nmax*nworlds*1 * sizeof(double);
  return bytes;
}
