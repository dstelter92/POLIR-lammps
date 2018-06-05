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
   Contributing author: David Stelter (BU)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "pair_polir.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairPolir::PairPolir(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairPolir::~PairPolir()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(c16);
    memory->destroy(c14);
    memory->destroy(c12);
    memory->destroy(c6);
    memory->destroy(polir1);
    memory->destroy(polir2);
    memory->destroy(polir3);
    memory->destroy(polir4);
  }
}

/* ---------------------------------------------------------------------- */

void PairPolir::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,r4,r10,r16inv,r17inv,forcepolir,factor_polir;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_polir = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) 
        r = sqrt(rsq);
        r4 = rsq*rsq;
        r10 = r4*r4*rsq;
        r16inv = 1.0/(r10*r4*rsq);
        r17inv = 1.0/(r10*r4*rsq*r);

        fpair = polir1[itype][jtype] + polir2[itype][jtype]*rsq + 
          polir3[itype][jtype]*r4 + polir4[itype][jtype]*r10;
        fpair *= factor_polir*r17inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = c16[itype][jtype] + c14[itype][jtype]*rsq + 
            c12[itype][jtype]*r4 + c6[itype][jtype]*r10;
          evdwl *= factor_polir*r16inv;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairPolir::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(c16,n+1,n+1,"pair:c16");
  memory->create(c14,n+1,n+1,"pair:c14");
  memory->create(c12,n+1,n+1,"pair:c12");
  memory->create(c6,n+1,n+1,"pair:c6");
  memory->create(polir1,n+1,n+1,"pair:polir1");
  memory->create(polir2,n+1,n+1,"pair:polir2");
  memory->create(polir3,n+1,n+1,"pair:polir3");
  memory->create(polir4,n+1,n+1,"pair:polir4");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPolir::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPolir::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double c16_one = force->numeric(FLERR,arg[2]);
  double c14_one = force->numeric(FLERR,arg[3]);
  double c12_one = force->numeric(FLERR,arg[4]);
  double c6_one = force->numeric(FLERR,arg[5]);

  double cut_one = cut_global;
  if (narg == 7) cut_one = force->numeric(FLERR,arg[6]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      c16[i][j] = c16_one;
      c14[i][j] = c14_one;
      c12[i][j] = c12_one;
      c6[i][j] = c6_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPolir::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  polir1[i][j] = 16.0 * c16[i][j];
  polir2[i][j] = 14.0 * c14[i][j]; 
  polir3[i][j] = 12.0 * c12[i][j];
  polir4[i][j] = 6.0 * c6[i][j];

  cut[j][i] = cut[i][j];
  polir1[j][i] = polir1[i][j];
  polir2[j][i] = polir2[i][j];
  polir3[j][i] = polir3[i][j];
  polir4[j][i] = polir4[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairPolir::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&c16[i][j],sizeof(double),1,fp);
        fwrite(&c14[i][j],sizeof(double),1,fp);
        fwrite(&c12[i][j],sizeof(double),1,fp);
        fwrite(&c6[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPolir::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&c16[i][j],sizeof(double),1,fp);
          fread(&c14[i][j],sizeof(double),1,fp);
          fread(&c12[i][j],sizeof(double),1,fp);
          fread(&c6[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&c16[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&c14[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&c12[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&c6[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairPolir::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPolir::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairPolir::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,c16[i][i],c14[i][i],c12[i][i],c6[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairPolir::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,
              c16[i][j],c14[i][j],c12[i][j],c6[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairPolir::single(int i, int j, int itype, int jtype, double rsq,
                            double factor_coul, double factor_polir,
                            double &fforce)
{
  double forcepolir,phi;
  double r,r4,r10,r16inv,r17inv;

  r = sqrt(rsq);
  r4 = rsq*rsq;
  r10 = r4*r4*rsq;
  r16inv = 1.0/(r10*r4*rsq);
  r17inv = 1.0/(r10*r4*rsq*r);

  forcepolir = polir1[itype][jtype] + polir2[itype][jtype]*rsq + 
      polir3[itype][jtype]*r4 + polir4[itype][jtype]*r10;
  forcepolir *= factor_polir*r17inv;

  phi = c16[itype][jtype] + c14[itype][jtype]*rsq + 
    c12[itype][jtype]*r4 + c6[itype][jtype]*r10;

  return factor_polir*phi*r16inv;
}
