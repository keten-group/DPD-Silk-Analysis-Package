/* ----------------------------------------------------------------------
   // Caution, a modification in force field made by Seunghwa Ryu
   on Aug 08, 2013, to describe hydrogen bonding
   version 1, single type case //
   
   Modified again on Jan 09, 2013, to describe hydrogen bonding better
   with overall r0 = 0.778, E0 = 2.85 
   (rin, E0 in this file should be 0.756, 3.53, respectively )

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_soft.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSoft::PairSoft(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairSoft::~PairSoft()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(prefactor);
    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */

void PairSoft::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double r,rsq,arg,factor_lj,rin,rout,rin2,rout2,droi,k,k0,EB;
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
  rin=0.756;rout=1.0;
  rin2=0.571536;rout2=1.0;
  droi=rout2-rin2;
  k0 = 24.0*rin2/droi/droi;
  //  k0 = 45.213379469434827; // don't know why k=k*EB makes problem in C++;

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
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        arg = MY_PI*r/cut[itype][jtype];
        EB = prefactor[itype][jtype];
        k = k0*EB;
        
        if (r < rin)
          fpair = factor_lj * k * (rin-r)/(0.000000000001+r);
        else
          fpair = factor_lj*12.0*EB*r*(rsq-rin2)*(rsq-rout2)/droi/droi/droi/r;
        
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag)
          if ( r < rin)
            evdwl = factor_lj * 0.5 * k * (r-rin) * (r-rin) - EB;
          else
            evdwl = (-1.0)* factor_lj * EB * (rout2-rsq)*(rout2-rsq)*
                    (rout2+2*rsq-3*rin2)/droi/droi/droi;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSoft::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(prefactor,n+1,n+1,"pair:prefactor");
  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSoft::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSoft::coeff(int narg, char **arg)
{
  if (narg < 3 || narg > 4)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double prefactor_one = force->numeric(FLERR,arg[2]);

  double cut_one = cut_global;
  if (narg == 4) cut_one = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      prefactor[i][j] = prefactor_one;
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

double PairSoft::init_one(int i, int j)
{
  // always mix prefactors geometrically

  if (setflag[i][j] == 0) {
    prefactor[i][j] = sqrt(prefactor[i][i]*prefactor[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  prefactor[j][i] = prefactor[i][j];
  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSoft::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&prefactor[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSoft::read_restart(FILE *fp)
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
          fread(&prefactor[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&prefactor[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSoft::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSoft::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairSoft::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g\n",i,prefactor[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairSoft::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g\n",i,j,prefactor[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */
double PairSoft::single(int i, int j, int itype, int jtype, double rsq,
                        double factor_coul, double factor_lj,
                        double &fforce)
{
  double r,arg,philj,rin,rout,rin2,rout2,droi,k,k0,EB;
    
  r = sqrt(rsq);
  arg = MY_PI*r/cut[itype][jtype];
    
  rin=0.756;rout=1.0;
  rin2=0.571536;rout2=1.0;
  droi=rout2-rin2;
  k0 = 24.0*rin2/droi/droi;
    
  EB=prefactor[itype][jtype];
  k=k0*EB;
  
  if (r < rin)
    fforce = factor_lj * k * (rin-r)/(0.000000000001+r);
  else
    fforce = factor_lj*12.0*EB*r*(rsq-rin2)*(rsq-rout2)/droi/droi/droi/r;

  if ( r < rin)
    philj = 0.5 * k * (r-rin) * (r-rin) - EB;
  else
    philj = (-1.0) * EB * (rout2-rsq)*(rout2-rsq)*
            (rout2+2*rsq-3*rin2)/droi/droi/droi;
  
  return factor_lj*philj;
}

/* ---------------------------------------------------------------------- */

void *PairSoft::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"a") == 0) return (void *) prefactor;
  return NULL;
}
