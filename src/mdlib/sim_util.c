/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */
static char *SRCID_sim_util_c = "$Id$";
#include <stdio.h>
#include <time.h>
#include "typedefs.h"
#include "string2.h"
#include "smalloc.h"
#include "names.h"
#include "confio.h"
#include "mvdata.h"
#include "txtdump.h"
#include "pbc.h"
#include "vec.h"
#include "time.h"
#include "nrnb.h"
#include "mshift.h"
#include "mdrun.h"
#include "update.h"
#include "physics.h"
#include "main.h"
#include "network.h"
#include "calcmu.h"
#include "constr.h"

#define difftime(end,start) ((double)(end)-(double)(start))

void print_time(FILE *out,time_t start,int step,t_inputrec *ir)
{
  static real time_per_step;
  static time_t end;
  time_t finish;
  double dt;
  char buf[48];

  fprintf(out,"\rstep %d",step);
  if ((step >= ir->nstlist)) {
    if ((ir->nstlist == 0) || ((step % ir->nstlist) == 0)) {
      /* We have done a full cycle let's update time_per_step */
      end=time(NULL);
      dt=difftime(end,start);
      time_per_step=dt/(step+1);
    }
    dt=(ir->nsteps-step)*time_per_step;

    if (dt >= 300) {    
      finish = end+(time_t)dt;
      sprintf(buf,"%s",ctime(&finish));
      buf[strlen(buf)-1]='\0';
      fprintf(out,", will finish at %s",buf);
    }
    else
      fprintf(out,", remaining runtime: %5d s               ",(int)dt);
  }
#ifdef USE_MPI
  if (gmx_parallel)
    fprintf(out,"\n");
#endif
  fflush(out);
}

time_t print_date_and_time(FILE *log,int nodeid,char *title)
{
  int i;
  char *ts,time_string[STRLEN];
  time_t now;

  now=time(NULL);
  ts=ctime(&now);
  for (i=0; ts[i]>=' '; i++) time_string[i]=ts[i];
  time_string[i]='\0';
  fprintf(log,"%s on node %d %s\n",title,nodeid,time_string);
  return now;
}

static void pr_commrec(FILE *log,t_commrec *cr)
{
  fprintf(log,"commrec: nodeid=%d, nnodes=%d, left=%d, right=%d, threadid=%d, nthreads=%d\n",
	  cr->nodeid,cr->nnodes,cr->left,cr->right,cr->threadid,cr->nthreads);
}

static void sum_forces(int start,int end,rvec f[],rvec flr[])
{
  int i;
  
  for(i=start; (i<end); i++)
    rvec_inc(f[i],flr[i]);
}

static void reset_energies(t_grpopts *opts,t_groups *grp,
			   t_forcerec *fr,bool bNS,real epot[])
{
  int   i,j;
  
  /* First reset all energy components but the Long Range, except in
   * some special cases.
   */
  for(i=0; (i<egNR); i++)
    if (((i != egLR) && (i != egLJLR)) ||
	(fr->bTwinRange && bNS) || (!fr->bTwinRange))
      for(j=0; (j<grp->estat.nn); j++)
	grp->estat.ee[i][j]=0.0;
  
  /* Normal potential energy components */
  for(i=0; (i<=F_EPOT); i++)
    epot[i] = 0.0;
  epot[F_DVDL]    = 0.0;
  epot[F_DVDLKIN] = 0.0;
}

/* 
 * calc_f_el calculates forces due to an electric field.
 *
 * force is kJ mol^-1 nm^-1 = e * kJ mol^-1 nm^-1 / e 
 *
 * Et[] contains the parameters for the time dependent 
 * part of the field (not yet used). 
 * Ex[] contains the parameters for
 * the spatial dependent part of the field. You can have cool periodic
 * fields in principle, but only a constant field is supported
 * now. 
 * The function should return the energy due to the electric field
 * (if any) but for now returns 0.
 */

static real calc_f_el(int start,int homenr,real charge[],rvec f[],t_cosines Ex[])
{
  real Emu,fmu,strength;
  int  i,m;
  
  Emu = 0;
  for(m=0; (m<DIM); m++)
    if (Ex[m].n) {
      /* Convert the field strength from V/nm to MD-units */
      strength = Ex[m].a[0]*FIELDFAC;
      for(i=start; (i<start+homenr); i++) {
	fmu      = charge[i]*strength;
	f[i][m] += fmu;
      } 
    }
  
  return Emu;
}

void do_force(FILE *log,t_commrec *cr,
	      t_parm *parm,t_nsborder *nsb,tensor vir_part,tensor pme_vir,
	      int step,t_nrnb *nrnb,t_topology *top,t_groups *grps,
	      rvec x[],rvec v[],rvec f[],rvec buf[],
	      t_mdatoms *mdatoms,real ener[],bool bVerbose,
	      real lambda,t_graph *graph,
	      bool bNS,bool bNBFonly,t_forcerec *fr, rvec mu_tot,
	      bool bGatherOnly)
{
  static rvec box_size;
  static real dvdl_lr = 0;
  int    cg0,cg1,i,j;
  int    start,homenr;
  static real mu_and_q[DIM+1]; 
  real   qsum;
  
  start  = START(nsb);
  homenr = HOMENR(nsb);
  cg0    = CG0(nsb);
  cg1    = CG1(nsb);
  
  update_forcerec(log,fr,parm->box);

  /* Calculate total (local) dipole moment in a temporary common array. 
   * This makes it possible to sum them over nodes faster.
   */
  calc_mu_and_q(nsb,x,mdatoms->chargeT,mu_and_q,mu_and_q+DIM);

  if (fr->ePBC != epbcNONE) { 
    /* Compute shift vectors every step, because of pressure coupling! */
    if (parm->ir.epc != epcNO)
      calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
    
    if (bNS) { 
      put_charge_groups_in_box(log,cg0,cg1,FALSE,
			       parm->box,box_size,&(top->blocks[ebCGS]),x,
			       fr->cg_cm);
      inc_nrnb(nrnb,eNR_RESETX,homenr);
    } else if (parm->ir.eI==eiSteep || parm->ir.eI==eiCG)
      unshift_self(graph,parm->box,x);

  }
  else if (bNS)
    calc_cgcm(log,cg0,cg1,&(top->blocks[ebCGS]),x,fr->cg_cm);
 
  if (bNS) {
    inc_nrnb(nrnb,eNR_CGCM,cg1-cg0);
    if (PAR(cr))
      move_cgcm(log,cr,fr->cg_cm,nsb->workload);
    if (debug)
      pr_rvecs(debug,0,"cgcm",fr->cg_cm,nsb->cgtotal);
  }
  
  /* Communicate coordinates and sum dipole and net charge if necessary */
  if (PAR(cr)) {
    move_x(log,cr->left,cr->right,x,nsb,nrnb);
    gmx_sum(DIM+1,mu_and_q,cr);
  }
  for(i=0;i<DIM;i++)
    mu_tot[i]=mu_and_q[i];
  qsum=mu_and_q[DIM];
  
  /* Reset energies */
  reset_energies(&(parm->ir.opts),grps,fr,bNS,ener);    
  if (bNS) {
    if (fr->ePBC != epbcNONE)
      /* Calculate intramolecular shift vectors to make molecules whole */
      mk_mshift(log,graph,parm->box,x);
	       
    /* Reset long range forces if necessary */
    if (fr->bTwinRange) {
      clear_rvecs(nsb->natoms,fr->f_twin);
      clear_rvecs(SHIFTS,fr->fshift_twin);
    }
    /* Do the actual neighbour searching and if twin range electrostatics
     * also do the calculation of long range forces and energies.
     */
    dvdl_lr = 0; 

    ns(log,fr,x,f,parm->box,grps,&(parm->ir.opts),top,mdatoms,
       cr,nrnb,nsb,step,lambda,&dvdl_lr);
  }
  /* Reset PME/Ewald forces if necessary */
  if (EEL_LR(fr->eeltype)) 
    clear_rvecs(homenr,fr->f_pme+start);
    
  /* Copy long range forces into normal buffers */
  if (fr->bTwinRange) {
    for(i=0; i<nsb->natoms; i++)
      copy_rvec(fr->f_twin[i],f[i]);
    for(i=0; i<SHIFTS; i++)
      copy_rvec(fr->fshift_twin[i],fr->fshift[i]);
  } 
  else {
    clear_rvecs(nsb->natoms,f);
    clear_rvecs(SHIFTS,fr->fshift);
  }
  
  /* Compute the forces */    
  force(log,step,fr,&(parm->ir),&(top->idef),nsb,cr,nrnb,grps,mdatoms,
	top->atoms.grps[egcENER].nr,&(parm->ir.opts),
	x,f,ener,bVerbose,parm->box,lambda,graph,&(top->atoms.excl),
	bNBFonly,pme_vir,mu_tot,qsum,bGatherOnly);
	
  /* Take long range contribution to free energy into account */
  ener[F_DVDL] += dvdl_lr;
  
  /* Compute forces due to electric field */
  calc_f_el(start,homenr,mdatoms->chargeT,f,parm->ir.ex);
  
#ifdef DEBUG
  if (bNS)
    print_nrnb(log,nrnb);
#endif

  /* The short-range virial from surrounding boxes */
  clear_mat(vir_part);
  calc_vir(log,SHIFTS,fr->shift_vec,fr->fshift,vir_part,cr);
  inc_nrnb(nrnb,eNR_VIRIAL,SHIFTS);

  if (debug) 
    pr_rvecs(debug,0,"vir_shifts",vir_part,DIM);

  /* When using PME/Ewald we compute the long range virial (pme_vir) there.
   * otherwise we do it based on long range forces from twin range
   * cut-off based calculation (or not at all).
   */
  
  /* Communicate the forces */
  if (PAR(cr))
    move_f(log,cr->left,cr->right,f,buf,nsb,nrnb);
}

void sum_lrforces(rvec f[],t_forcerec *fr,int start,int homenr)
{
  /* Now add the forces from the PME calculation. Since this only produces
   * forces on the local atoms, this can be safely done after the
   * communication step.
   */
  if (EEL_LR(fr->eeltype))
    sum_forces(start,start+homenr,f,fr->f_pme);
}

void calc_virial(FILE *log,int start,int homenr,rvec x[],rvec f[],
		 tensor vir_part,tensor pme_vir,
		 t_commrec *cr,t_graph *graph,matrix box,
		 t_nrnb *nrnb,t_forcerec *fr,bool bTweak)
{
  int i,j;
  tensor virtest;
  
  /* Now it is time for the short range virial. At this timepoint vir_part
   * already contains the virial from surrounding boxes.
   * Calculate partial virial, for local atoms only, based on short range. 
   * Total virial is computed in global_stat, called from do_md 
   */
  f_calc_vir(log,start,start+homenr,x,f,vir_part,cr,graph,box);
  inc_nrnb(nrnb,eNR_VIRIAL,homenr);

  /* Add up the long range forces if necessary */
  /* if (!bTweak) {
    sum_lrforces(f,fr,start,homenr);
    }*/
  
  /* Add up virial if necessary */  
  if (EEL_LR(fr->eeltype) && (fr->eeltype != eelPPPM)) {
    if (debug && bTweak) {
      clear_mat(virtest);
      f_calc_vir(log,start,start+homenr,x,fr->f_pme,virtest,cr,graph,box);
      pr_rvecs(debug,0,"virtest",virtest,DIM);
      pr_rvecs(debug,0,"pme_vir",pme_vir,DIM);
    }    
    /* PPPM virial sucks */
    if (!bTweak)
      for(i=0; (i<DIM); i++) 
	for(j=0; (j<DIM); j++) 
	  vir_part[i][j]+=pme_vir[i][j];
  }
  if (debug)
    pr_rvecs(debug,0,"vir_part",vir_part,DIM);
}

#ifdef NO_CLOCK 
#define clock() -1
#endif
static double runtime=0;
static clock_t cprev;

void start_time(void)
{
  cprev   = clock();
  runtime = 0.0;
}

void update_time(void)
{
  clock_t c;
  
  c        = clock();
  runtime += (c-cprev)/(double)CLOCKS_PER_SEC;
  cprev    = c;
}

double node_time(void)
{
  return runtime;
}

void do_shakefirst(FILE *log,bool bTYZ,real lambda,real ener[],
		   t_parm *parm,t_nsborder *nsb,t_mdatoms *md,
		   rvec x[],rvec vold[],rvec buf[],rvec f[],
		   rvec v[],t_graph *graph,t_commrec *cr,t_nrnb *nrnb,
		   t_groups *grps,t_forcerec *fr,t_topology *top,
		   t_edsamyn *edyn,t_pull *pulldata)
{
  int    i,m,start,homenr,end,step;
  tensor shake_vir;
  double mass,tmass,vcm[4];
  real   dt=parm->ir.delta_t;
  real   dt_1;

  if (count_constraints(top,cr)) {
    start  = START(nsb);
    homenr = HOMENR(nsb);
    end    = start+homenr;
    if (debug)
      fprintf(debug,"vcm: start=%d, homenr=%d, end=%d\n",start,homenr,end);
    /* Do a first SHAKE to reset particles... */
    step = -2;
    fprintf(log,"\nConstraining the starting coordinates (step %d)\n",step);
    clear_mat(shake_vir);
    update(nsb->natoms,start,homenr,step,lambda,&ener[F_DVDL],
	   parm,1.0,md,x,graph,
	   NULL,NULL,vold,NULL,x,top,grps,shake_vir,cr,nrnb,bTYZ,
	   FALSE,edyn,pulldata,FALSE);
    /* Compute coordinates at t=-dt, store them in buf */
    /* for(i=0; (i<nsb->natoms); i++) {*/
    for(i=start; (i<end); i++) {
      for(m=0; (m<DIM); m++) { 
	f[i][m]=x[i][m];
	buf[i][m]=x[i][m]-dt*v[i][m];
      }
    }
    
    /* Shake the positions at t=-dt with the positions at t=0
     * as reference coordinates.
     */
    step = -1;
    fprintf(log,"\nConstraining the coordinates at t0-dt (step %d)\n",step);
    clear_mat(shake_vir);
    update(nsb->natoms,start,homenr,
	   step,lambda,&ener[F_DVDL],parm,1.0,md,f,graph,
	   NULL,NULL,vold,NULL,buf,top,grps,shake_vir,cr,nrnb,bTYZ,FALSE,
	   edyn,pulldata,FALSE);
    
    /* Compute the velocities at t=-dt/2 using the coordinates at
     * t=-dt and t=0
     * Compute velocity of center of mass and total mass
     */
    for(m=0; (m<4); m++)
      vcm[m] = 0;
    dt_1=1.0/dt;
    for(i=start; (i<end); i++) {
      /*for(i=0; (i<nsb->natoms); i++) {*/
      mass = md->massA[i];
      for(m=0; (m<DIM); m++) {
	v[i][m]=(x[i][m]-f[i][m])*dt_1;
	vcm[m] += v[i][m]*mass;
      }
      vcm[3] += mass;
    }
    /* Compute the global sum of vcm */
    if (debug)
      fprintf(debug,"vcm: %8.3f  %8.3f  %8.3f,"
	      " total mass = %12.5e\n",vcm[XX],vcm[YY],vcm[ZZ],vcm[3]);
    if (PAR(cr))
      gmx_sumd(4,vcm,cr);
    tmass = vcm[3];
    for(m=0; (m<DIM); m++)
      vcm[m] /= tmass;
    if (debug) 
      fprintf(debug,"vcm: %8.3f  %8.3f  %8.3f,"
	      " total mass = %12.5e\n",vcm[XX],vcm[YY],vcm[ZZ],tmass);
    /* Now we have the velocity of center of mass, let's remove it */
    for(i=start; (i<end); i++) {
      for(m=0; (m<DIM); m++)
	v[i][m] -= vcm[m];
    }
  }
}

void calc_dispcorr(FILE *log,int eDispCorr,t_forcerec *fr,int natoms,
		   matrix box,tensor pres,tensor virial,real ener[])
{
  static bool bFirst=TRUE;
  real vol,rc3,spres,svir;
  int  m;
  
  ener[F_DISPCORR] = 0.0;
  ener[F_PRES]     = trace(pres)/3.0;
  
  if (eDispCorr != edispcNO) {
    vol           = det(box);
    /* Forget the (small) effect of the shift on the LJ energy *
     * when fr->bLJShift = TRUE                                */  
    rc3              = fr->rvdw*fr->rvdw*fr->rvdw;
    ener[F_DISPCORR] = -2.0*natoms*natoms*M_PI*fr->avcsix/(3.0*vol*rc3);
    if (eDispCorr == edispcEnerPres) {
      spres            = 2.0*ener[F_DISPCORR]*PRESFAC/vol;
      svir             = -6.0*ener[F_DISPCORR];
      ener[F_PRES]    += spres;
      for(m=0; (m<DIM); m++) {
	pres[m][m]    += spres;
	virial[m][m]  += svir;
      }
    }
    else {
      spres = 0;
      svir  = 0;
    }
    if (bFirst) {
      if (eDispCorr == edispcEner)
	fprintf(log,"Long Range LJ corr. to Epot: %10g\n",ener[F_DISPCORR]);
      else if (eDispCorr == edispcEnerPres)
	fprintf(log,"Long Range LJ corr. to Epot: %10g, Pres: %10g, Vir: %10g\n",
		ener[F_DISPCORR],spres,svir);
      bFirst = FALSE;
    }
  }
  ener[F_EPOT] += ener[F_DISPCORR];
  ener[F_ETOT] += ener[F_DISPCORR];
}

