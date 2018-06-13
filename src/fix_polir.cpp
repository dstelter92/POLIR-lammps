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

#define INVOKED 1
#define POLIR_DEBUG 0

/* ---------------------------------------------------------------------- */

FixPolir::FixPolir(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 11) error->all(FLERR,"Illegal fix polir command");

  //local_flag = 1;

  // fix/polir input variables
  typeH = force->inumeric(FLERR,arg[3]);  // atom type for H
  typeO = force->inumeric(FLERR,arg[4]);  // atom type for O
  qeH = force->numeric(FLERR,arg[5]);     // equilibrium charge on H
  reOH = force->numeric(FLERR,arg[6]);    // equilibrium O-H bond length
  uH = force->numeric(FLERR,arg[7]);      // H dipole along HOH bisector
  uO = force->numeric(FLERR,arg[8]);      // O dipole along OH bond
  alphaH = force->numeric(FLERR,arg[9]);  // polarizibility for typeH
  alphaO = force->numeric(FLERR,arg[10]); // polarizibility for type0

  // Constants 
  // see fix_polir.h for details

  // Thole damping parameters
  // these are NOT inputs to fix/polir since they are constant
  // in the model. Only requires small changes if the model is improved
  CD_intra_OH = 0.650;
  CD_intra_HH = 0.050;
  DD_intra_OH = 0.690;
  DD_intra_HH = 0.050;
  CC_inter = 0.50;
  CD_inter = 0.15;
  DD_inter = 0.30;

  // bond length coeff
  c1 = 0.120;
  c3 = -0.136;

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
  memory->destroy(thole);
  delete [] id_q;
  delete [] id_thole;
  //delete [] id_lbond;
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
  nmax = atom->nmax;
  count = 0;


  // MPI init
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  
  // memory management
  memory->create(charges,nmax,"polir:charges");
  memory->create(thole,nmax,7,"polir:thole");

  // Search for charge compute
  compute_id = -1;
  for (i=0; i<modify->ncompute; i++) {
    if (strcmp(modify->compute[i]->style,"POLIR/CHARGE/ATOM") == 0) {
      compute_id = i;
      break;
    }
  }
  // If doesn't exist, create a new one
  if (compute_id < 0) {
    int n = strlen(id) + 4;
    id_q = new char[n];
    strcpy(id_q,id);
    strcat(id_q,"_polir_charge_atom");

    char **newarg = new char*[4];
    newarg[0] = id_q;
    newarg[1] = group->names[igroup];
    newarg[2] = (char *) "POLIR/CHARGE/ATOM";
    newarg[3] = id;

    modify->add_compute(4,newarg);
    delete [] newarg;

    compute_id = modify->ncompute -1;
  }
  compute_pca = modify->compute[compute_id];


  // Search for thole compute
  compute_id = -1;
  for (i=0; i<modify->ncompute; i++) {
    if (strcmp(modify->compute[i]->style,"POLIR/THOLE/LOCAL") == 0) {
      compute_id = i;
      break;
    }
  }
  // If doesn't exist, create a new one
  if (compute_id < 0) {
    int n = strlen(id) + 4;
    id_thole = new char[n];
    strcpy(id_thole,id);
    strcat(id_thole,"_polir_thole_local");

    char **newarg = new char*[4];
    newarg[0] = id_thole;
    newarg[1] = group->names[igroup];
    newarg[2] = (char *) "POLIR/THOLE/LOCAL";
    newarg[3] = id;

    modify->add_compute(4,newarg);
    delete [] newarg;

    compute_id = modify->ncompute - 1;
  }
  compute_thole = modify->compute[compute_id];
  

  /*
  // compute_bond_local
  compute_id = -1;
  for (i=0; i<modify->ncompute; i++) {
    if (strcmp(modify->compute[i]->style,"bond/local") == 0) {
      compute_id = i;
      break;
    }
  }
  if (compute_id < 0) {
    int n = strlen(id) + 4;
    id_lbond = new char[n];
    strcpy(id_lbond,id);
    strcat(id_lbond,"_bond_local");

    char **newarg = new char*[4];
    newarg[0] = id_lbond;
    newarg[1] = group->names[igroup];
    newarg[2] = (char *) "bond/local";
    newarg[3] = (char *) "dist";

    modify->add_compute(4,newarg);
    delete [] newarg;

    compute_id = modify->ncompute - 1;
  }
  compute_lbond = modify->compute[compute_id];
  */
}

/* ---------------------------------------------------------------------- */

void FixPolir::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) {
    post_force(vflag);

    // Invoke all computes to run on each step
    compute_pca->invoked_flag |= INVOKED;
    //compute_lbond->invoked_flag |= INVOKED;
    compute_thole->invoked_flag |= INVOKED;
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
  int *mask = atom->mask;
  int *tag = atom->tag;
  double *q = atom->q;
  int i;

  nlocal = atom->nlocal;
  nmax = atom->nmax;

  allocate();

  // Re-calculate atomic charges based on changes from last step
  // invoke per-atom fix
  /*
  compute_lbond->compute_local();
  lbond = compute_lbond->vector_atom;

  for (i=0; i<nlocal; i++) {
    if (mask[i] & groupbit) {
      fprintf(screen,"proc%d atom%d lbond[%d]=%g\n",me,tag[i],i,lbond[i]);
    }
  }
  */

  compute_pca->compute_peratom();
  charges = compute_pca->vector_atom;

  // set atomic charges
  for (i=0; i<nlocal; i++) {
    if (mask[i] & groupbit) {
      if (POLIR_DEBUG)
        fprintf(screen,"proc:%d atom%d:%g\n",me,tag[i],charges[i]);
      q[i] = charges[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPolir::post_force(int vflag)
{
  int *mask = atom->mask;
  double **f = atom->f;
  int i;
  
  nlocal = atom->nlocal;
  nmax = atom->nmax;
  
  compute_thole->compute_local();
  thole = compute_thole->array_local;
  
  allocate();

}

/* ---------------------------------------------------------------------- */

void FixPolir::end_of_step()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  count++;

}

/* ---------------------------------------------------------------------- */

void FixPolir::allocate()
{
  nmax = atom->nmax;
  memory->destroy(charges);
  memory->create(charges,nmax,"polir:charges");
  memory->create(thole,nmax,7,"polir:thole");
}

/* ---------------------------------------------------------------------- */

double FixPolir::memory_usage()
{
  nmax = atom->nmax;
  double bytes = nmax*8 * sizeof(double);
  return bytes;
}

/* -------------------------------------------------------------------------
   extract input parameters and defined constants from Fix
------------------------------------------------------------------------- */

void *FixPolir::extract(const char *str, int &dim)
{
  dim=0;
  // input parameters
  if (strcmp(str,"typeH") == 0)
    return &typeH;
  else if (strcmp(str,"typeO") == 0)
    return &typeO;
  else if (strcmp(str,"qeH") == 0)
    return &qeH;
  else if (strcmp(str,"reOH") == 0)
    return &reOH;
  else if (strcmp(str,"uH") == 0)
    return &uH;
  else if (strcmp(str,"uO") == 0)
    return &uO;
  else if (strcmp(str,"alphaH") == 0)
    return &alphaH;
  else if (strcmp(str,"alphaO") == 0)
    return &alphaO;
  // return POLIR specific constants
  else if (strcmp(str,"c1") == 0)
    return &c1;
  else if (strcmp(str,"c3") == 0)
    return &c3;
  else if (strcmp(str,"CD_intra_OH") == 0)
    return &CD_intra_OH;
  else if (strcmp(str,"CD_intra_HH") == 0)
    return &CD_intra_HH;
  else if (strcmp(str,"DD_intra_OH") == 0)
    return &DD_intra_OH;
  else if (strcmp(str,"DD_intra_HH") == 0)
    return &DD_intra_HH;
  else if (strcmp(str,"CC_inter") == 0)
    return &CC_inter;
  else if (strcmp(str,"CD_inter") == 0)
    return &CD_inter;
  else if (strcmp(str,"DD_inter") == 0)
    return &DD_inter;

  return NULL;
}
