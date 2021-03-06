//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code 
// contributors.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================
//! \file fluxrope_lf.cpp
//  \brief Problem generator for the fluxrope based on Lin-Forbes model.
//
// REFERENCE:
//  Lin & Forbes 2000
//==============================================================================

// C/C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

// Normalization parameters
Real Lchar, Bchar, Nechar, Rhochar, Pchar, Timechar, Tchar, Vchar;
Real Gsun_nondim;

// Constant parameters
Real Kb = 1.38e-23, Mp = 1.6726e-27, Mu0 = 1.25663706144e-6, LnA = 30.0;
Real gravitational_const = 6.672e-11; /* (N M^2 kg^-2)*/
Real solarmass = 1.99e+30;            /* (kg) */
Real solarradius = 6.96e+8;           /* (m) */

// functions
static Real func_pini(const Real x, const Real y);
static Real func_rhoini(const Real x, const Real y);
static Real func_teini(const Real x, const Real y);

static Real func_rho_pdepend(const Real y, const Real y0, const Real rho0, const Real T0);
static Real func_p_betadepend(Real y, const Real y0, const Real rho0, const Real T0);
static Real func_p_isoteconst(const Real x2, const Real y0,
                               const Real p0, const Real te0);
static Real presini_integral(const Real y, const Real y0, const Real p0);


Real func_azini(Real x, Real y);
Real func_bzini(Real x, Real y);

Real func_integ_bx(Real y, const Real xfix);
Real func_integ_by(Real x, const Real yfix);
Real func_integ_pphi(Real r, const Real phi);
Real func_integ_pphi_bg(Real r, const Real phi);
Real func_integ_pphi_rope(Real r, const Real phi);
Real func_integ_py(Real y, Real xfix);

static Real func_bphi(const Real r);
static Real func_bphi_bg(const Real r);
static Real func_bphi_rope(const Real r);

static Real func_rjphi(const Real r);
static Real func_rjphi_bg(const Real r);
static Real func_rjphi_rope(const Real r);

static Real func_uphi(const Real r);
static Real func_uphi_bg(const Real r);
static Real func_uphi_rope(const Real r);

static Real func_bmx(const Real x, const Real y);
static Real func_bmy(const Real x, const Real y);
static Real func_bmx_bg(const Real x, const Real y, const Real ls, const Real le, const int type);
static Real func_bmy_bg(const Real x, const Real y, const Real ls, const Real le, const int type);
static Real func_pbypxini(const Real x, const Real y);

double adaptiveSimpsons(double (*f)(double, double),   // ptr to function
                        double param, 
                        double a, double b,  // interval [a,b]
                        double epsilon,  // error tolerance
                        int maxRecursionDepth);
double adaptiveSimpsonsAux(double (*f)(double, double),
                          double param,
                          double a, double b,
                          double epsilon,
                          double S,
                          double fa, double fb, double fc, int bottom);

// Global array
AthenaArray<Real> az;

// Global parameters to define the initial fluxrope
static int sw_frbz, fr_mode;
static Real fr_d, fr_h, fr_ri, fr_del, fr_rja;
static Real gauss_fr, gauss_erope;

// Global parameters to define the source
static Real erope_rja, erope_ri;
static Real I_erope_ratio;
static Real I_ratio, I_source;
static Real lambda_sta, lambda_end; // X range of source distribution

// Global parameters to define the background gas
static Real scale_bgdens, scale_lowtcore;
static Real Te_corona, Te_bottom;
static Real y_base, p_base, rho_base, te_base;
static Real beta_floor;
static Real height_bottom, width_bottom;

// made global to share with BC functions
static Real gsun_y(const Real x2);         // Gsun(r)

// Global paramters to define the pre-eruption phase
static Real pre_etime;

// Global parameter to define refinements
static Real dplim_refine, jzlim_refine;

// Boundary conditions
void SymmInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, 
                FaceField &b, Real time, Real dt, 
                int is, int ie, int js, int je, int ks, int ke, int ngh);
void OpenOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                FaceField &b, Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh);
void LintInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void OpenOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, 
                FaceField &b, Real time, Real dt, 
                int is, int ie, int js, int je, int ks, int ke, int ngh);
void ReduceInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh);
void ReduceOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh);
// Refinement
int RefinementCondition(MeshBlock *pmb);

// Gravity source
void static_grav_source(MeshBlock *pmb, const Real time, const Real dt, 
                        const AthenaArray<Real> &prim,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

// Reduce Momenotum source
void initial_slowv_source(MeshBlock *pmb, const Real time, const Real dt,
                        const AthenaArray<Real> &prim,
                        const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

//==============================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also
//  be used to initialize variables which are global to (and therefore can be 
//  passed to) other functions in this file.  Called in Mesh constructor.
//==============================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Enroll boundary value function pointers
  EnrollUserBoundaryFunction(INNER_X1, SymmInnerX1);
  EnrollUserBoundaryFunction(OUTER_X1, OpenOuterX1);
  EnrollUserBoundaryFunction(INNER_X2, LintInnerX2);
  EnrollUserBoundaryFunction(OUTER_X2, OpenOuterX2);
  //EnrollUserBoundaryFunction(INNER_X3, ReduceInnerX3);
  //EnrollUserBoundaryFunction(OUTER_X3, ReduceOuterX3);
  
  // Diffusion
  //EnrollFieldDiffusivity(VariDiffusivity);

  if(adaptive==true)
      EnrollUserRefinementCondition(RefinementCondition);

  // Static gravity source in momentum equations
  EnrollUserExplicitSourceFunction(static_grav_source);
  
  // Pre-initialize phase
  // EnrollUserExplicitSourceFunction(initial_slowv_source);
  return;
}

//==============================================================================
//! \fn void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
//  \brief Function to clarify user outputs.
//==============================================================================
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(3);
  return;
}

//==============================================================================
//! \fn void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
//  \brief Function to define user outputs.
//==============================================================================
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        Real pmag = 0.5*(SQR(pfield->bcc(IB1,k,j,i))
                            +SQR(pfield->bcc(IB2,k,j,i))
                            +SQR(pfield->bcc(IB3,k,j,i)));
        user_out_var(0,k,j,i) = phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);
        user_out_var(1,k,j,i) = phydro->w(IPR,k,j,i)/pmag;
        user_out_var(2,k,j,i) = 0;
      }
    }
  }
  return;
}

//==============================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//==============================================================================
// This function is called at the end of every timestep 
// (note: it is not called at the half timestep). 
void MeshBlock::UserWorkInLoop(void) {
  /* (1) remove vz component 
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        //phydro->u(IEN,k,j,i) = phydro->u(IEN,k,j,i)
        //  -0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IM3,k,j,i) = 0.0;
      }
    }
  }
  */
  return;
}

//==============================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the flarecs.  The initial conditions are
//  constructed assuming the domain extends over [-1.0x1.0, 0.0x2.0].
//==============================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // Read input arguments
  Real gamma_const = 5./3.;
  
  // Read input parameters
  Lchar = pin->GetOrAddReal("problem", "Lchar", 1.0e8); // (200Mm)
  Bchar = pin->GetOrAddReal("problem", "Bchar", 0.01);  // (100G)
  Nechar = pin->GetOrAddReal("problem", "Nechar", 1.0e16); // (10^10 cm-3)

  Pchar = Bchar * Bchar / Mu0;
  Rhochar = Nechar * Mp;
  Tchar = Pchar / (Nechar * 2. * Kb); /* Total number is 2*ne */
  Vchar = Bchar / sqrt(Mu0 * Rhochar);
  Timechar = Lchar / Vchar;

  // Initialize the flux rop
  sw_frbz = pin->GetOrAddInteger("problem", "sw_frbz", 1);
  fr_mode = pin->GetOrAddInteger("problem", "fr_mode", 2);
  fr_d  = pin->GetOrAddReal("problem", "fr_d", 0.075);
  fr_h  = pin->GetOrAddReal("problem", "fr_h", 0.51);
  fr_ri = pin->GetOrAddReal("problem", "fr_ri", 0.05);
  fr_del = pin->GetOrAddReal("problem", "fr_del", 0.025);
  fr_rja = pin->GetOrAddReal("problem", "fr_rja", 200.0);
  gauss_fr = fr_ri*(2./4.29193);
  Real I_rope; 
  if (fr_mode == 1) {
    I_rope = PI * fr_rja*(SQR(fr_ri) + 0.25*SQR(fr_del) - 2.*SQR(fr_del / PI));
  } else if (fr_mode == 2) {
    I_rope = fr_rja*2.0*PI*SQR(gauss_fr)* (1.0 - exp(-0.5*SQR(10.0/gauss_fr)));
  } else {
    I_rope = 0;
  }

  // Initialize sources
  I_ratio = pin->GetOrAddReal("problem", "I_ratio", 0.90);
  I_source = I_rope/I_ratio;
  lambda_sta = pin->GetOrAddReal("problem", "lambda_sta", 0);
  lambda_end = pin->GetOrAddReal("problem", "lambda_end", 0.1);

  // Initialize the Extra-current around the fluxROPE
  erope_ri = pin->GetOrAddReal("problem", "erope_ri", 0.05);
  gauss_erope = erope_ri*(2./4.29193);
  erope_rja = pin->GetOrAddReal("problem", "erope_rja", 0);
  Real I_erope = erope_rja*2.0*PI*SQR(gauss_erope)* (1.0 - exp(-0.5*SQR(10.0/gauss_erope)));
  
  // Initialize background gas
  beta_floor = pin->GetOrAddReal("problem", "beta_floor", 0.005);
  scale_bgdens = pin->GetOrAddReal("problem", "scale_bgdens", 0.2);
  Te_corona = pin->GetOrAddReal("problem", "Te_corona", 2.0e6); // (unit: K)
  Te_bottom = pin->GetOrAddReal("problem", "Te_bottom", 1.0e4); // (unit: K)
  height_bottom = pin->GetOrAddReal("problem", "height_bottom", 0.1); // (non-dim)
  width_bottom = pin->GetOrAddReal("problem", "width_bottom", 0.05); // (non-dim)
  
  // Define the base density and temperature
  y_base = 0.0;
  rho_base = 1.0*scale_bgdens;
  te_base = Te_corona / Tchar;
  p_base = rho_base * te_base;
  
  // Define the pre-eruption phase
  pre_etime = pin->GetOrAddReal("problem","pre_etime", 0.1);

  // Define refinement parameters
  dplim_refine = pin->GetOrAddReal("problem","dplim_refine", 2.5);
  jzlim_refine = pin->GetOrAddReal("problem","jzlim_refine", 10.0);

  //if (my_rank == 0) {
    printf("---------------\n");
    printf("sw_frbz=%d\n", sw_frbz);
    printf("h = %f\n", fr_h);
    printf("d = %f\n", fr_d);
    printf("ri = %f\n", fr_ri);
    printf("rja = %f\n", fr_rja);
    printf("I_rope = %f\n", I_rope);
    printf("erope_ri =%f\n", erope_ri);
    printf("erope_rja =%f\n", erope_rja);
    printf("I_erope = %f\n", I_erope);
    printf("scale_bgdens =%f\n", scale_bgdens);
    printf("scale_lowtcore =%f\n", scale_lowtcore);
    printf("rho_base=%f, Te_base=%f, p_base=%f\n", rho_base, te_base, p_base);
    printf("Lchar=%.3e, Vchar=%.3e, Timechar=%.3e, Tchar=%.3e, Nechar=%.3e, Bchar=%.3e\n", Lchar, Vchar, Timechar, Tchar, Nechar, Bchar);
    printf("---------------\n");
  //}

  // Define the local temporary array: az & pgas
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  Real xc, yc;
  AthenaArray<Real> pbzf;
  AthenaArray<Real> rhof;
  AthenaArray<Real> pf;
  pbzf.NewAthenaArray(nx2+1,nx1+1);
  rhof.NewAthenaArray(nx2+1,nx1+1);
  pf.NewAthenaArray(nx2+1,nx1+1);
  az.NewAthenaArray(nx2+1,nx1+1);

  // Initialize vector potential
  for (int j=js; j<=je+1; ++j) {
  for (int i=is; i<=ie+1; ++i) {
    az(j,i) = func_azini(pcoord->x1f(i), pcoord->x2f(j));
  }}

  // Initialize the pressure_fluxrop, density and pressure at the cell interface
  for (int j=js; j<=je+1; ++j) {
  for (int i=is; i<=ie+1; ++i) {
    // Pressure for Bz
    Real rf = sqrt(pow(pcoord->x1f(i),2) + pow(pcoord->x2f(j)-fr_h, 2));
    if (sw_frbz == 1) {
      pbzf(j,i) = fabs(func_uphi_rope(rf));
    } else {
      pbzf(j,i) = 0;
    }
    // Density
    rhof(j,i) = func_rhoini(pcoord->x1f(i), pcoord->x2f(j));
    // Pressure
    pf(j,i) = func_pini(pcoord->x1f(i), pcoord->x2f(j));
  }}

  // Initialize density, momentum, face-centered fields
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real rhoavg = 0.25*(rhof(j,i) + rhof(j,i+1) + rhof(j+1,i) + rhof(j+1,i+1));
    phydro->u(IDN,k,j,i) = rhoavg;
    phydro->u(IM1,k,j,i) = 0.0;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;
  }}}

  // initialize interface B
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
    //pfield->b.x1f(k,j,i) = (az(j+1,i) - az(j,i))/pcoord->dx2f(j);
    pfield->b.x1f(k,j,i) = func_bmx(pcoord->x1f(i),
                                    0.5*(pcoord->x2f(j)+pcoord->x2f(j+1)));
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    //pfield->b.x2f(k,j,i) = (az(j,i) - az(j,i+1))/pcoord->dx1f(i);
    pfield->b.x2f(k,j,i) = func_bmy(0.5*(pcoord->x1f(i)+pcoord->x1f(i+1)),
                                    pcoord->x2f(j));
  }}}
  
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
      Real p_fluxrope_avg = 0.25*(pbzf(j,i) + pbzf(j,i+1)
                               + pbzf(j+1,i) + pbzf(j+1,i+1));
      pfield->b.x3f(k,j,i) = sqrt(2.0*p_fluxrope_avg);
  }}}

  // initialize total energy
  if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      Real pavg = 0.25*(pf(j,i) + pf(j,i+1) + pf(j+1,i) + pf(j+1,i+1));
      phydro->u(IEN,k,j,i) = pavg/(gamma_const - 1.0) +
          0.5*(SQR(0.5*(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))) +
               SQR(0.5*(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i))) +
               SQR(0.5*(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i)))) + (0.5)*
          (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i))
           + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
    }}}
  }

  // Release memory
  pbzf.DeleteAthenaArray();
  rhof.DeleteAthenaArray();
  pf.DeleteAthenaArray();
  az.DeleteAthenaArray();
  return;
}

//==============================================================================
// function: initial p, te, and rho(x, y)
//==============================================================================
static Real func_pini(const Real x, const Real y) {
  Real p;
  Real r = sqrt(x*x + (y-fr_h)*(y-fr_h));

  Real y0 = 0.5;
  Real p0 = func_p_isoteconst(y0, y_base, p_base, te_base);

  p = presini_integral(y, y0, p0);
  
  if (sw_frbz == 0) {
    p = p + fabs(func_uphi_rope(r));
  }
  if (erope_rja >= 0.0) {
    p = p + fabs(func_uphi_bg(r));
  }
  
  return p;
}

static Real func_rhoini(const Real x, const Real y) {
  Real rho;
  Real y0 = 0.5;
  Real p0 = func_p_isoteconst(y0, y_base, p_base, te_base);
  Real p_bg = presini_integral(y, y0, p0);
  Real rho0 = p_bg/func_teini(x, y);
  
  // case 1
  rho = rho0*pow(func_pini(x, y)/p_bg, 0.6);
  /*
  // case 2
  Real r = sqrt(x*x + (y-fr_h)*(y-fr_h));
  Real p_temp = p_bg + fabs(func_uphi_bg(r)); // do not include the realpfluxrope
  Real rho_temp = rho0*pow(p_temp/p_bg, 0.6);
  Real te_temp = p_temp/rho_temp;
  rho = func_pini(x, y)/te_temp;
  */
  // case 3
  // rho = func_pini(x, y)/func_teini(x, y);
  return rho;
}

static Real func_teini(const Real x, const Real y)
{
  Real Te, t1, t2;
  
  /* Temperature distribution in y-direction: in non-dimensional forms */
  Real tecor = Te_corona / Tchar, tebot = Te_bottom / Tchar;
  t1 = 0.5 * (tecor - tebot);
  t2 = 0.5 * (tecor + tebot);
  Te = t1 * tanh((y - height_bottom) / width_bottom) + t2;
  
  return Te;
}

/*----------------------------------------------------------------------------*/
/*  \fn gas pressure assuming constant temperature */
static Real func_p_isoteconst(const Real y, const Real y0,
                              const Real p0, const Real te0)
{
  Real pgas, gm, r, r0;
  
  gm = gravitational_const * solarmass * Rhochar / (Lchar * Pchar);
  r = solarradius / Lchar + y;
  r0 = solarradius / Lchar + y0;
  pgas = p0 * exp(gm / te0 * (1.0 / r - 1.0 / r0));
  
  return pgas;
}

/*----------------------------------------------------------------------------*/
/*  \fn gas pressure performing numerical integration */
static Real presini_integral(const Real y, const Real y0, const Real p0)
{
  Real pgas, Iout;
  Iout = adaptiveSimpsons(func_integ_py, 
                          0.0,
                          y0, y, 1.0e-9, 10000);
  pgas = p0 * exp(Iout);
  return pgas;
}

//==============================================================================
// function: initial az(x, y)
//==============================================================================
Real func_azini(Real x, Real y) {  
  Real Ixx, Iyy, az_out;
  Real az0=0.0, x0=0.5, y0=fr_h;
  // Integrate (y=y0 -> y) at x = x0 in the y-direction
  Iyy = adaptiveSimpsons(func_integ_bx,
                        x0, 
                        y0, y, 1.0e-9, 1000000);

  // Integrate (x=x0 -> x) at y = y in the x-direction
  Ixx = adaptiveSimpsons(func_integ_by,
                        y,
                        x0, x, 1.0e-9, 1000000);

  az_out = Ixx + Iyy;
  return az_out;
}

/*
//==============================================================================
// function: initial bz(x, y)
//==============================================================================
Real func_bzini(Real x, Real y) {
  Real r = sqrt(x*x + (y-fr_h)*(y-fr_h));
  Real bz, p_fluxrope;

  if (sw_frbz == 0) {
    bz = 0;
  } else {
    p_fluxrope = -func_uphi_rope(r);
    bz = sqrt(2.0*p_fluxrope);
  }
  
  return bz;
}
*/

//==============================================================================
//  Functions for fluxrope
//==============================================================================
static Real func_bmx(const Real x, const Real y)
{
  /*c model field x-component */
  Real rs, rm;
  Real bmx_rope, bmx_bg;
  rs = sqrt(pow(x, 2) + pow(y - fr_h, 2));
  rm = sqrt(pow(x, 2) + pow(y + fr_h, 2));
  
  // flux rope (and its mirros)
  if (rs > 0.0) {
    bmx_rope = +func_bphi(rs)*(y-fr_h)/rs - func_bphi(rm)*(y+fr_h)/rm;
  } else {
    bmx_rope = 0.0;
  }
  
  // background
  bmx_bg = func_bmx_bg(x, y, lambda_sta, lambda_end, 2);
  
  return bmx_rope + bmx_bg;
}

static Real func_bmx_bg(const Real x, const Real y, const Real ls, const Real le, const int type)
{
  int i;
  Real yc = y+fr_d;
  Real bmx_bg = 0;
  Real coeff;
  if (type == 0){
    // single source
    coeff = -0.25*I_source/PI;
    Real lmx = lambda_end-x;
    Real lpx = lambda_end+x;
    Real numera = lmx*(pow(yc,2)+pow(lpx,2)) + lpx*(pow(yc,2)+pow(lmx,2));
    Real denomi = (pow(yc,2)+pow(lmx,2))*(pow(yc,2)+pow(lpx,2));
    bmx_bg = coeff*(numera/denomi);
  } else if (type == 2){
    // line source
    coeff = (0.125*I_source/PI)/(le-ls);
    Real x4 = pow(x,4) + 2.0*pow(x*yc, 2) + pow(yc,4);
    Real log_0 = -log(pow(le,4) - 2.0*pow(le,2)*(x*x-yc*yc) + x4);
    Real log_1 = +log(pow(ls,4) - 2.0*pow(ls,2)*(x*x-yc*yc) + x4);
    bmx_bg = coeff*(log_0 + log_1);
  } else {
    bmx_bg = 0;
  }
  return bmx_bg;
}

static Real func_bmy(const Real x, const Real y)
{
  /*  model field y-component */
  Real rs, rm;
  Real bmy_rope, bmy_bg;
  rs = sqrt(pow(x, 2) + pow(y - fr_h, 2));
  rm = sqrt(pow(x, 2) + pow(y + fr_h, 2));
  
  // fluxrope
  if (rs > 0.0) {
    bmy_rope = -func_bphi(rs)*x/rs + func_bphi(rm)*x/rm;
  } else {
    bmy_rope = 0;
  }
  
  // background
  bmy_bg = func_bmy_bg(x, y, lambda_sta, lambda_end, 2);
  
  return bmy_rope + bmy_bg;
}

static Real func_bmy_bg(const Real x, const Real y, const Real ls, const Real le, const int type)
{
  Real yc = y+fr_d;
  Real bmy_bg = 0;
  Real coeff;
  if (type == 0){
    // single source
    coeff = -0.25*I_source/PI;
    Real lmx = lambda_end-x;
    Real lpx = lambda_end+x;
    Real numera = yc*(pow(lmx,2) - pow(lpx,2));
    Real denomi = (pow(yc,2) + pow(lmx,2))*(pow(yc,2) + pow(lpx,2));
    bmy_bg = coeff*numera/denomi;
  } else if (type == 2){
    // line source
    Real inner_0 = 0.5*(le*le-x*x+yc*yc)/(x*yc);
    Real inner_1 = 0.5*(ls*ls-x*x+yc*yc)/(x*yc);
    coeff = (0.25*I_source/PI)/(le-ls);
    bmy_bg = coeff*(atan(inner_0) - atan(inner_1));
  } else {
    bmy_bg = 0;
  }
  return bmy_bg;
}
  
static Real func_bphi(const Real r)
{
  return func_bphi_rope(r) + func_bphi_bg(r);
}

static Real func_bphi_rope(const Real r)
{
  Real bphi, bphi_rope;
  if (fr_mode == 1) {
    /* mode 1: Cylindrical field function */
    Real riq, delq, piq, t1, t2, t3, bphi;
    Real ro, r1, r2;
    r1 = fr_ri - 0.5 * fr_del;
    r2 = fr_ri + 0.5 * fr_del;
    riq = fr_ri * fr_ri;
    delq = fr_del * fr_del;
    piq = PI * PI;
    if (r <= r1) {
      bphi_rope = -0.5 * fr_rja * r;
    } else if (r <= r2) {
      t1 = 0.5 * r1 * r1 - delq / piq + 0.5 * r * r;
      t2 = (fr_del * r / PI) * sin((PI / fr_del) * (r - r1));
      t3 = (delq / piq) * cos((PI / fr_del) * (r - r1));
      bphi_rope = -0.5 * fr_rja * (t1 + t2 + t3) / r;
    } else {
      bphi_rope = -0.5 * fr_rja * (riq + 0.25 * delq - 2. * delq / piq) / r;
    }
  } else if (fr_mode == 2) {
    /* mode 2: Gaussian distribtion */
    Real I_rope = fr_rja * 2.0 * PI * SQR(gauss_fr) * (1.0 - exp(-0.5 * SQR(r / gauss_fr)));
    if (r > 0.0) {
      bphi_rope = -I_rope / (2.0 * PI * r);
    } else {
      bphi_rope = 0.0;
    }
  } else {
    bphi_rope = 0.0;
  }
  return bphi_rope;
}

static Real func_bphi_bg(const Real r)
{
  /* bphi_bg */
  Real bphi_bg;
  Real Ibg = erope_rja * 2.0 * PI * SQR(gauss_erope) * (1.0 - exp(-0.5 * SQR(r / gauss_erope)));
  if (r > 0.0) {
    bphi_bg = -Ibg / (2.0 * PI * r);
  } else {
    bphi_bg = 0.0;
  }
  return bphi_bg;
}

static Real func_rjphi(const Real r)
{
  return func_rjphi_rope(r)+func_rjphi_bg(r);
}

static Real func_rjphi_rope(const Real r)
{
  /*  current density */
  Real rjphi_rope;
  if (fr_mode == 1) {
    /* Case 1: Cylindrical rope */
    Real r1, r2;
    r1 = fr_ri - 0.5 * fr_del;
    r2 = fr_ri + 0.5 * fr_del;
    if (r <= r1) {
      rjphi_rope = fr_rja;
    } else if (r <= r2) {
      rjphi_rope = 0.5 * fr_rja * (cos((PI / fr_del) * (r - r1)) + 1.);
    } else {
      rjphi_rope = 0.0;
    }
  } else if (fr_mode == 2) {
    /* Case 2: Gaussion distribution */
    rjphi_rope = fr_rja * exp(-0.5 * pow((r / gauss_fr), 2));
  } else {
    rjphi_rope = 0;
  }
  return rjphi_rope;
}

static Real func_rjphi_bg(const Real r)
{
  /*  current density */
  Real rjphi_bg;
  rjphi_bg = erope_rja*exp(-0.5*pow((r/gauss_erope), 2));
  return rjphi_bg;
}

static Real func_uphi(const Real r)
{
  Real uphi;
  Real rend;
  if (fr_mode == 1) {
    rend = fr_ri + fr_del;
  } else if (fr_mode == 2) {
    rend = 2.0;
  } else {
    rend = 1.0;
  }
  if (r <= rend) {
    uphi = adaptiveSimpsons(func_integ_pphi,
                            0,
                            r, rend, 1.0e-9, 1000000);
  } else {
    uphi = 0;
  }
  return uphi;
}

static Real func_uphi_rope(const Real r)
{
  Real uphi;
  Real rend;
  if (fr_mode == 1) {
    rend = fr_ri + fr_del;
  } else if (fr_mode == 2) {
    rend = 2.0;
  } else {
    rend = 1.0;
  }
  if (r <= rend) {
    uphi = adaptiveSimpsons(func_integ_pphi_rope,
                            0,
                            r, rend, 1.0e-9, 1000000);
  } else {
    uphi = 0;
  }
  return uphi;
}

static Real func_uphi_bg(const Real r)
{
  Real uphi;
  Real rend = 10.0;
  if (r <= rend) {
    uphi = adaptiveSimpsons(func_integ_pphi_bg,
                            0,
                            r, rend, 1.0e-9, 1000000);
  } else {
    uphi = 0;
  }
  return uphi;
}

//==============================================================================
// function: integrations
//==============================================================================
Real func_integ_bx(Real y, const Real xfix) {
  return func_bmx(xfix, y);
}

Real func_integ_by(Real x, const Real yfix) {
  return -func_bmy(x, yfix);
}

Real func_integ_pphi(Real r, const Real phi) {
  return func_rjphi(r) * func_bphi(r);
}

Real func_integ_pphi_rope(Real r, const Real phi) {
  return func_rjphi_rope(r) * func_bphi(r);
}

Real func_integ_pphi_bg(Real r, const Real phi) {
  return func_rjphi_bg(r) * func_bphi(r);
}

Real func_integ_py(Real y, Real xfix)
{
  return -gsun_y(y) / func_teini(xfix, y);
}

//==============================================================================
// Adaptive Simpson's Rule
//==============================================================================
double adaptiveSimpsons(double (*f)(double, double),   // ptr to function
                        double param,
                        double a, double b, // interval [a,b]
                        double epsilon,  // error tolerance
                        int maxRecursionDepth) {   // recursion cap
  double c = (a + b)/2, h = b - a;
  double fa = f(a, param), fb = f(b, param), fc = f(c, param);
  double S = (h/6)*(fa + 4*fc + fb);
  return adaptiveSimpsonsAux(f,
                            param,
                            a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}

//==============================================================================
// Recursive auxiliary function for adaptiveSimpsons() function below
//==============================================================================
double adaptiveSimpsonsAux(double (*f)(double, double),
                          double param,
                          double a, double b,
                          double epsilon,
                          double S,
                          double fa, double fb, double fc, int bottom) {
  
  double c = (a + b)/2, h = b - a;
  double d = (a + c)/2, e = (c + b)/2;
  double fd = f(d, param), fe = f(e, param);
  double Sleft = (h/12)*(fa + 4*fd + fc);
  double Sright = (h/12)*(fc + 4*fe + fb);
  double S2 = Sleft + Sright;
  if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon) // magic 15 comes 
                                                 // from error analysis
    return S2 + (S2 - S)/15;
  return adaptiveSimpsonsAux(f,
                            param,
                            a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1)
        +adaptiveSimpsonsAux(f,
                            param,
                            c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);
}

//==============================================================================
// SymmInnerX1 boundary condution
//==============================================================================
void SymmInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, 
                FaceField &b, Real time, Real dt, 
                int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        prim(n,k,j,is-i) = prim(n,k,j,is+i-1);
      }
    }}
  }

  // Set velocity Vx
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        prim(IVX,k,j,is-i) = -prim(IVX,k,j,is+i-1);
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        b.x1f(k,j,(is-i)) = b.x1f(k,j,is+i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
      for (int i=1; i<=ngh; ++i) {
        b.x2f(k,j,(is-i)) = -b.x2f(k,j,is+i-1);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        b.x3f(k,j,(is-i)) = b.x3f(k,j,is+i-1);
      }
    }}
  }
}

//==============================================================================
// Open boudnary condition at the right edge
//==============================================================================
//! \fn void OpenOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                         FaceField &b, Real time, Real dt,
//                         int is, int ie, int js, int je, int ks, int ke, int ngh)
//  \brief Open boundary conditions, outer x1 boundary
void OpenOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,ie+i) = prim(n,k,j,ie-i+1);
        }
      }
    }
  }

  // inflow restriction
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        if (prim(IVX,k,j,ie+i) < 0.0) {
          prim(IVX,k,j,ie+i) = 0.0;
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
      for (int i=1; i<=ngh; ++i) {
        b.x2f(k,j,(ie+i)) = 2.0*b.x2f(k,j,(ie+i-1))-b.x2f(k,j,(ie+i-2));
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+i))
        -(pco->dx1f(ie+i)/pco->dx2f(j))
        *(b.x2f(k,(j+1),(ie+i)) - b.x2f(k,j,(ie+i)));
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        b.x3f(k,j,(ie+i)) = b.x3f(k,j,(ie-i+1));
      }
    }}
  }

  return;
}

//==============================================================================
// LintInnerX2 boundary condition
//==============================================================================
//! \fn void LintInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                       FaceField &b, Real time, Real dt,
//                       int is, int ie, int js, int je, int ks, int ke, int ngh)
//  \brief Linetied boundary conditions at the bottom.
void LintInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, 
                    FaceField &b, Real time, Real dt, 
                    int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // (a) First extroplate all primary values
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          prim(n,k,js-j,i) = prim(n,k,js,i);
        }
      }
    }
  }
  // (b) Set velocity and density
  Real xc, yc;
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=is; i<=ie; ++i) {
        prim(IVX,k,js-j,i) = 0;
        prim(IVY,k,js-j,i) = 0;
        prim(IVZ,k,js-j,i) = 0;
        // case 1: depends on p
        xc = 0.5*(pco->x1f(i) + pco->x1f(i+1));
        yc = 0.5*(pco->x2f(js-j) + pco->x2f(js-j+1));
        //Real grav_acc = -gsun_y(yc);
        //prim(IPR,k,js-j,i) = prim(IPR,k,js+j-1,i) - prim(IDN,k,js+j-1,i)*grav_acc*(2*j-1)*pco->dx2f(j);
        
        // case 2: fixed p & rho
        prim(IDN,k,js-j,i) = func_rhoini(xc, yc);
        prim(IPR,k,js-j,i) = func_pini(xc, yc);
      }
    }
  }
  // (c) Set face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    Real az_c, az_p, xc;
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          xc = 0.5*(pco->x1f(i) + pco->x1f(i+1));
          b.x2f(k,(js-j),i) = func_bmy(xc, pco->x2f(js-j));
        }
      }
    }
    Real pbypx;
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pbypx = func_pbypxini(pco->x1f(i), pco->x2f(js-j+1));
          b.x1f(k,(js-j),i) = b.x1f(k,(js-j+1),i) - pbypx*pco->dx2v(js-j+1);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          b.x3f(k,(js-j),i) = b.x3f(k,js,i);
        }
      } 
    }
  }
  return;
}

static Real func_pbypxini(const Real x, const Real y)
{
  // line source
  Real yc = y + fr_d;
  Real le = lambda_end, ls = lambda_sta;
  Real coeff = 0.125*I_source/(PI*(le-ls));
  Real xles = x*x*yc*yc + 0.25*pow(le*le-x*x+yc*yc, 2);
  Real xlss = x*x*yc*yc + 0.25*pow(ls*ls-x*x+yc*yc, 2);
  Real numera_1 = xles*(ls*ls+x*x+yc*yc);
  Real numera_2 = xlss*(le*le+x*x+yc*yc);
  Real numera = yc*(numera_1 - numera_2);
  Real denomi = xles*xlss;
  Real pbypx = coeff*numera/denomi;

  // Test
  //Real ddx = 1.0e-12;
  //pbypx = (func_bmy(x+ddx,y)-func_bmy(x-ddx,y))/(2.0*ddx);

  return pbypx;
}

//==============================================================================
// OpenOuterX2 boundary condition
//==============================================================================
void OpenOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                FaceField &b, Real time, Real dt,
                int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          prim(n,k,je+j,i) = prim(n,k,je-j+1,i);
          Real yc = 0.5*(pco->x2f(je+j) + pco->x2f(je+j+1));
          Real grav_acc = -gsun_y(yc);
          prim(IPR,k,je+j,i) = prim(IPR,k,je-j+1,i)
             + prim(IDN,k,je-j+1,i)*grav_acc*(2*j-1)*pco->dx2f(j);
        }
      }
    }
  }

  // Inflow restriction
  Real dn_ratio;
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (prim(IVY,k,je+j,i) < 0.0) {
          prim(IVY,k,je+j,i) = 0.0;
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          b.x1f(k,(je+j),i) = 2.0*b.x1f(k,(je+j-1),i) - b.x1f(k,(je+j-2),i);
        }
      }
    }

    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          b.x2f(k,(je+j+1),i) = b.x2f(k,(je+j),i)
          -pco->dx2f(je+j)/pco->dx1f(i)*(b.x1f(k,(je+j),i+1)-b.x1f(k,(je+j),i));
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          b.x3f(k,(je+j  ),i) = b.x3f(k,(je-j+1),i);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReduceInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke, int ngh)
//  \brief Reduce boundary conditions, inner x3 boundary

void ReduceInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,ks-k,j,i) = prim(n,ks,j,i);
        }
      }
    }
  }
  
  // Reduce v3
  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        prim(IVZ,ks-k,j,i) = 0;
      }
    }
  }
  
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f((ks-k),j,i) = b.x1f(ks,j,i);
        }
      }}
    
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f((ks-k),j,i) = b.x2f(ks,j,i);
        }
      }}
    
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f((ks-k),j,i) = b.x3f(ks,j,i);
        }
      }}
  }
  
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReduceOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke, int ngh)
//  \brief Reduce boundary conditions, outer x3 boundary

void ReduceOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,ke+k,j,i) = prim(n,ke,j,i);
        }
      }
    }
  }
  
  // Reduce v3 (or pressure)
  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        prim(IVZ,ke+k,j,i) = 0;
      }
    }
  }
  
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f((ke+k  ),j,i) = b.x1f((ke  ),j,i);
        }
      }}
    
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f((ke+k  ),j,i) = b.x2f((ke  ),j,i);
        }
      }}
    
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f((ke+k+1),j,i) = b.x3f((ke+1),j,i);
        }
      }}
  }
  
  return;
}

//------------------------------------------------------------------------------
// Magnetic diffusivity from variable coefficients
//------------------------------------------------------------------------------
/*void VariDiffusivity(FieldDiffusion *pfdif, MeshBlock *pmb, const AthenaArray<Real> &w,
                      const AthenaArray<Real> &bmag, const int is, const int ie, const int js,
                      const int je, const int ks, const int ke) {
  if (pfdif->eta_ohm > 0.0) { // Ohmic resistivity is turned on
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
#pragma omp simd
        for(int i=is; i<=ie; i++)
          pfdif->etaB(I_O, k,j,i) = pfdif->eta_ohm;
      }
    }
  }
  
  return;
}
*/

//------------------------------------------------------------------------------
// Static gravity source in y-direction
//------------------------------------------------------------------------------
void static_grav_source(MeshBlock *pmb, const Real time, const Real dt, 
  const AthenaArray<Real> &prim,
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real yc = 0.5*(pmb->pcoord->x2f(j) + pmb->pcoord->x2f(j+1));
        Real grav_acc = -gsun_y(yc); // negative direction in y-
        Real src = dt*prim(IDN,k,j,i)*grav_acc;
        cons(IM2,k,j,i) += src;
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVY,k,j,i);
        }
      }
    }
  return;
}

//------------------------------------------------------------------------------
// Pre-Initialize: slow-v
//------------------------------------------------------------------------------
void initial_slowv_source(MeshBlock *pmb, const Real time, const Real dt,
                          const AthenaArray<Real> &prim,
                          const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  /* reduce velocity components at the pre-eruption phase */
  if (time < pre_etime) {
    Real scale = (pre_etime - time)/pre_etime;
    Real ls = 1.0 - scale;
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          cons(IEN,k,j,i) = cons(IEN,k,j,i)
          -0.5*(SQR(scale*cons(IM1,k,j,i))
               +SQR(scale*cons(IM2,k,j,i))
               +SQR(scale*cons(IM3,k,j,i)))/cons(IDN,k,j,i);
          cons(IM1,k,j,i) = ls*cons(IM1,k,j,i);
          cons(IM2,k,j,i) = ls*cons(IM2,k,j,i);
          cons(IM3,k,j,i) = ls*cons(IM3,k,j,i);
        }
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real gsun_r(const Real x2)
 *  \brief Non-dimensional Gravity at height x2
 */
static Real gsun_y(const Real x2)
{
  Real g, gnondim;
  Real r;
  r = solarradius + x2 * Lchar;                    /* (m) */
  g = (gravitational_const * solarmass) / (r * r); /* (m s^-2)*/
  gnondim = g * Rhochar * Lchar / Pchar;
  //gnondim = g * Lchar / pow(Vchar, 2);
  return gnondim;
}

//------------------------------------------------------------------------------
// AMR refinement condition functions
//------------------------------------------------------------------------------
int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> &bx = pmb->pfield->b.x1f;
  AthenaArray<Real> &by = pmb->pfield->b.x2f;
  AthenaArray<Real> &bz = pmb->pfield->b.x3f;
  Real maxeps = 0;
  int k=pmb->ks;
  /* (1) By and Bx flag
  for(int j=pmb->js+1; j<=pmb->je-1; j++) {
    for(int i=pmb->is+1; i<=pmb->ie-1; i++) {
      Real bc = pow(bx(k,j,i)*bx(k,j,i)+by(k,j,i)*by(k,j,i), 0.5);
      Real epsbx= std::abs(bx(k,j+1,i)-2.0*bx(k,j,i)+bx(k,j-1,i))/bx(k,j,i);
      Real epsby= std::abs(by(k,j,i+1)-2.0*by(k,j,i)+by(k,j,i-1))/by(k,j,i);
      Real eps = std::max(epsbx*bc, epsby*bc);
      maxeps = std::max(maxeps, eps);
    }
  }
  if(maxeps > 1.0) return 1;
  if(maxeps < 0.5) return -1;
  return 0;
  */
  
  /* (2) Pressure flag */
  int pflag = 0;
  for(int j=pmb->js; j<=pmb->je; j++) {
    for(int i=pmb->is; i<=pmb->ie; i++) {
      Real eps= (std::abs(w(IEN,k,j,i+1)-2.0*w(IEN,k,j,i)+w(IEN,k,j,i-1))
                  +std::abs(w(IEN,k,j+1,i)-2.0*w(IEN,k,j,i)+w(IEN,k,j-1,i)))/w(IEN,k,j,i);
      maxeps = std::max(maxeps, eps);
    }
  }
  if(maxeps > dplim_refine) pflag = 1;
  if(maxeps < 0.2*dplim_refine) pflag = -1;
  
  /* (3) Jz_max */
  int jflag = 0;
  Real jz_max = 0;
  for(int j=pmb->js; j<=pmb->je; j++) {
    for(int i=pmb->is; i<=pmb->ie; i++) {
      Real jzc = (by(k,j,i) - by(k,j,i-1))/pmb->pcoord->dx1v(i)
               - (bx(k,j,i) - bx(k,j-1,i))/pmb->pcoord->dx2v(j);
      jz_max = std::max(jz_max, fabs(jzc));
    }
  }
  if(jz_max >= jzlim_refine) jflag = 1;
  if(jz_max <= 0.2*jzlim_refine) jflag = -1;

  /* (4) beta flag 
  Real minbeta=1.0;
  for(int j=pmb->js; j<=pmb->je; j++) {
    for(int i=pmb->is; i<=pmb->ie; i++) {
      Real bxc = bx(k,j,i);
      Real byc = by(k,j,i);
      Real bzc = bz(k,j,i);
      Real pmag = std::max(0.5*(bxc*bxc + byc*byc + bzc*bzc), 1.0e-7);
      Real beta = w(4,k,j,i)/pmag;
      minbeta = std::min(minbeta, beta);
    }
  }
  //int flag_beta = 0;
  //if (minbeta < 0.001) flag_beta = 1;
  //if (minbeta > 0.010) flag_beta = -1;
  */

  // return flag 
  if ((pflag == 1) || (jflag == 1)) return 1;
  if ((pflag == -1) && (jflag == -1)) return -1;
  return 0;
}
