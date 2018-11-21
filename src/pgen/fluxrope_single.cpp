//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code
// contributors.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================
//! \file fluxrope_single.c
//  \brief Problem generator for the single fluxrope test.
//
// REFERENCE: For example, see: Wang et al 2009, ApJ.
//==============================================================================

// C/C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
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

// functions
Real func_pini(Real x, Real y);
Real func_teini(Real x, Real y);
Real func_rhoini(Real x, Real y);
Real func_azini(Real x, Real y);
Real func_bzini(Real x, Real y);
Real func_pbypxini(Real x, Real y);
Real func_integ_bx(Real y, const Real xfix);
Real func_integ_by(Real x, const Real yfix);
Real func_integ_pphi(Real r, const Real phi);

static Real func_bphi(const Real r);
static Real func_rjphi(const Real r);
static Real func_uphi(const Real r);
static Real func_back(const Real r);
static Real func_bmx(const Real x1, const Real x2);
static Real func_bmy(const Real x1, const Real x2);

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

// Global parameters to define the initial fluxrope
static Real beta0, temperature0;
static int fr_case;
static Real fr_d, fr_h, fr_ri, fr_del, fr_rja;
static Real p_ambient=0.06, te_ambient=0.6, rho_ambient=0.1;

//==============================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also
//  be used to initialize variables which are global to (and therefore can be
//  passed to) other functions in this file.  Called in Mesh constructor.
//==============================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Enroll boundary value function pointers
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
//! \fn void MeshBlock::UserWorkInLoop(void)
//==============================================================================
// This function is called at the end of every timestep
// (note: it is not called at the half timestep).
void MeshBlock::UserWorkInLoop(void) {
  // (1) remove velocity vz components
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
      }
    }
  }
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
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the flarecs.  The initial conditions are
//  constructed assuming the domain extends over [-1.0x1.0, 0.0x2.0].
//==============================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  
  // Read input arguments
  Real gm1 = peos->GetGamma() - 1.0;
  Real gamma_const = 5./3.;
  
  // Initialize the flux rope
  fr_d = 0.0625;
  fr_h = 0.15;
  fr_ri = 0.01;
  fr_del = 0.0;
  fr_rja = 6000.0;
  
  // Define the local temporary array: az & pgas
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  AthenaArray<Real> az;
  az.NewAthenaArray(nx2,nx1);
  
  // Initialize vector potential
  for (int j=js; j<=je+1; ++j) {
    for (int i=is; i<=ie+1; ++i) {
      az(j,i) = func_azini(pcoord->x1f(i), pcoord->x2f(j));
    }}
  
  // Initialize density, momentum, face-centered fields
  Real xc, yc;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        xc = 0.5*(pcoord->x1f(i) + pcoord->x1f(i+1));
        yc = 0.5*(pcoord->x2f(j) + pcoord->x2f(j+1));
        phydro->u(IDN,k,j,i) = func_rhoini(xc, yc);
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
      }}}
  
  // initialize interface B
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie+1; i++) {
        //pfield->b.x1f(k,j,i) = (az(j+1,i) - az(j,i))/pcoord->dx2f(j);
        yc = 0.5*(pcoord->x2f(j) + pcoord->x2f(j+1));
        pfield->b.x1f(k,j,i) = func_bmx(pcoord->x1f(i), yc);
      }}}
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie; i++) {
        //pfield->b.x2f(k,j,i) = (az(j,i) - az(j,i+1))/pcoord->dx1f(i);
        xc = 0.5*(pcoord->x1f(i) + pcoord->x1f(i+1));
        pfield->b.x2f(k,j,i) = func_bmy(xc, pcoord->x2f(j));
      }}}
  for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        xc = 0.5*(pcoord->x1f(i) + pcoord->x1f(i+1));
        yc = 0.5*(pcoord->x2f(j) + pcoord->x2f(j+1));
        pfield->b.x3f(k,j,i) = func_bzini(xc, yc);
      }}}
  
  // initialize total energy
  if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          xc = 0.5*(pcoord->x1f(i) + pcoord->x1f(i+1));
          yc = 0.5*(pcoord->x2f(j) + pcoord->x2f(j+1));
          phydro->u(IEN,k,j,i) = func_pini(xc, yc)/gm1 +
          0.5*(SQR(0.5*(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))) +
               SQR(0.5*(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i))) +
               SQR(0.5*(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i)))) + (0.5)*
          (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i))
           + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }}}
  }
  
  az.DeleteAthenaArray();
  return;
}

//==============================================================================
// function: initial p, te, and rho(x, y)
//==============================================================================
Real func_pini(Real x, Real y) {
  Real r = sqrt(x*x + (y-fr_h)*(y-fr_h));
  Real p;
  //p = p_ambient - func_uphi(r);
  p = p_ambient;
  return p;
}

Real func_teini(Real x, Real y) {
  Real te;
  te = te_ambient;
  return te;
}

Real func_rhoini(Real x, Real y) {
  Real rho, p, te;
  // Get gas pressure and temperature
  p = func_pini(x, y);
  te = func_teini(x, y);
  
  /* Add a dense (low temperature) fluxrope core
   Real pi = 3.14159265358979;
   Real r, r1, r2, r_edge;
   Real t_inner, t_outer;
   t_outer = te;
   t_inner = t_outer*1.0e-2;
   r1 = fr_ri - 0.5*fr_del;
   r_edge = fr_del;
   r = sqrt(x*x + pow(y-fr_h,2));
   if (r <= r1) {
   te = t_inner;
   } else if (r <= r1 + r_edge) {
   te = t_inner + 0.5*(t_outer-t_inner)*(1.0-cos(pi*(r-r1)/r_edge));
   } else {
   te = te;
   }*/
  
  // Get non-dimensional density
  rho = p/te;
  return rho;
}

//==============================================================================
// function: initial az(x, y)
//==============================================================================
Real func_azini(Real x, Real y) {
  Real Ixx, Iyy, az_out;
  Real az0=0, x0=0.5, y0=0.5;
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

//==============================================================================
// function: initial bz(x, y)
//==============================================================================
Real func_bzini(Real x, Real y) {
  Real bz, p_fluxrope;
  Real r;
  r = sqrt(x*x + (y-fr_h)*(y-fr_h));
  p_fluxrope = -func_uphi(r);
  bz = sqrt(2.0*p_fluxrope);
  bz = 0;
  return bz;
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

//==============================================================================
//  Functions for fluxrope
//==============================================================================
static Real func_bmx(const Real x, const Real y) {
  /*c model field x-component */
  Real rs, rm;
  Real bmx;
  rs = sqrt(pow(x, 2) + pow(y - fr_h, 2));
  rm = sqrt(pow(x, 2) + pow(y + fr_h, 2));
  if (rs > 0.0) {
    bmx = +func_bphi(rs)*(y-fr_h)/rs - func_bphi(rm)*(y+fr_h)/rm;
  } else {
    bmx = 0.0;
  }
  return bmx;
}

static Real func_bmy(const Real x, const Real y) {
  /*  model field z-component */
  Real rs, rm;
  Real bmy;
  rs = sqrt(pow(x, 2) + pow(y - fr_h, 2));
  rm = sqrt(pow(x, 2) + pow(y + fr_h, 2));
  if (rs > 0.0) {
    bmy = -func_bphi(rs)*x/rs + func_bphi(rm)*x/rm;
  } else {
    bmy = 0.0;
  }
  return bmy;
}

static Real func_bphi(const Real r)
{
  /* cylindrical field function */
  Real riq, delq, piq, t1, t2, t3, bphi;
  Real pi = 3.14159265358979;
  Real ro, r1, r2;
  r1 = fr_ri - 0.5 * fr_del;
  r2 = fr_ri + 0.5 * fr_del;
  riq = fr_ri * fr_ri;
  delq = fr_del * fr_del;
  piq = pi * pi;
  if (r <= r1) {
    bphi = -0.5 * fr_rja * r;
  } else if (r <= r2) {
    t1 = 0.5 * r1 * r1 - delq / piq + 0.5 * r * r;
    t2 = (fr_del * r / pi) * sin((pi / fr_del) * (r - r1));
    t3 = (delq / piq) * cos((pi / fr_del) * (r - r1));
    bphi = -0.5 * fr_rja * (t1 + t2 + t3) / r;
  } else {
    bphi = -0.5 * fr_rja * (riq + 0.25 * delq - 2. * delq / piq) / r;
  }
  return bphi;
}

static Real func_uphi(const Real r)
{
  Real uphi;
  Real rend = fr_ri + fr_del;
  if (r < rend) {
    uphi = adaptiveSimpsons(func_integ_pphi,
                            0,
                            r, rend, 1.0e-9, 1000000);
  } else {
    uphi = 0;
  }
  return uphi;
}

static Real func_rjphi(const Real r)
{
  /*  current density */
  Real pi = 3.14159265358979;
  Real r1, r2;
  Real rjphi;
  r1 = fr_ri - 0.5 * fr_del;
  r2 = fr_ri + 0.5 * fr_del;
  if (r <= r1) {
    rjphi = fr_rja;
  } else if (r <= r2) {
    rjphi = 0.5 * fr_rja * (cos((pi / fr_del) * (r - r1)) + 1.);
  } else {
    rjphi = 0.;
  }
  return rjphi;
}

static Real func_back(const Real r)
{
  Real riq, delq, piq, rmm, back;
  Real pi = 3.14159265358979;
  riq = pow(fr_ri, 2);
  delq = pow(fr_del, 2);
  piq = pow(pi, 2);
  rmm = 0.5 * fr_rja * (riq + 0.25 * delq - 2. * delq / piq);
  back = rmm / pow(r, 3);
  return back;
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
