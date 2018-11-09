//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code 
// contributors.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================
//! \file fluxrope.c
//  \brief Problem generator for the fluxrope problem.
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
static Real gamma_const = 5./3.;
static int fr_case;
static Real fr_d, fr_h, fr_ri, fr_del, fr_rmom, fr_rja;

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
        phydro->u(IEN,k,j,i) = phydro->u(IEN,k,j,i)
          -0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IM3,k,j,i) = 0.0;
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
  
  // Initialize the flux rope
  fr_case = pin->GetInteger("problem","fr_case");
  if (fr_case == 1) {
    fr_d = 0.0625;
    fr_h = 0.25;
    fr_ri = 0.02;
    fr_del = 0.01;
    fr_rja = 2000.0;
    fr_rmom = 1.0;
  } else if (fr_case == 2) {
    fr_d  = pin->GetReal("problem","fr_d");
    fr_h = 1.25*fr_d;
    fr_ri = 0.15*fr_d;
    fr_del = 0.0;
    fr_rja = 3000.0;
    fr_rmom = fr_d*fr_d*125./32.* 0.80;
  } else {
    printf("Error: define fluxrope.");
  }

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
    pfield->b.x1f(k,j,i) = (az(j+1,i) - az(j,i))/pcoord->dx2f(j);
    //yc = 0.5*(pcoord->x2f(j) + pcoord->x2f(j+1));
    //pfield->b.x1f(k,j,i) = func_bmx(pcoord->x1f(i), yc);
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x2f(k,j,i) = (az(j,i) - az(j,i+1))/pcoord->dx1f(i);
    //xc = 0.5*(pcoord->x1f(i) + pcoord->x1f(i+1));
    //pfield->b.x2f(k,j,i) = func_bmy(xc, pcoord->x2f(j));
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
  p = 10.0/gamma_const;
  return p;
}

Real func_teini(Real x, Real y) {
  Real te;
  te = 0.5/gamma_const;
  return te;
}

Real func_rhoini(Real x, Real y) {
  Real rho, p, te;
  // Get gas pressure and temperature
  p = func_pini(x, y);
  te = func_teini(x, y);

  /* Add a dense (low temperature) fluxrope core */
  Real pi = 3.14159265358979;
  Real r, r1, r2, r_edge;
  Real t_inner, t_outer;
  Real bphi2, bphi2_max, ratio;
  t_outer = te;
  t_inner = t_outer*1.0e-2;
  r2 = fr_ri + 0.5*fr_del;
  r_edge = 0.1;
  r = sqrt(x*x + pow(y-fr_h,2)); 
  if (r <= r2) {
    te = t_inner;
  } else if (r <= r2 + r_edge) {
    te = t_inner + 0.5*(t_outer-t_inner)*(1.0-cos(pi*(r-r2)/r_edge));
  } else {
    te = te;
  }

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
static Real func_bmx(const Real x, const Real y)
{
  /*c model field x-component */
  Real rs, rm, rb, r2;
  Real bmx;
  rs = sqrt(pow(x, 2) + pow(y - fr_h, 2));
  rm = sqrt(pow(x, 2) + pow(y + fr_h, 2));
  rb = sqrt(pow(x, 2) + pow(y + fr_d, 2));
  r2 = fr_ri + 0.5*fr_del;
  if (rs > 0.0) {
    if (fr_case == 1) {
    /* dipole case */
    bmx = +func_bphi(rs)*(y-fr_h)/rs
        -func_bphi(rm)*(y+fr_h)/rm
        -func_bphi(r2)*fr_rmom*fr_d*r2*(pow(x,2)-pow(y+fr_d, 2))/pow(rb,4);
    } else {
    /* quadrupole */
    bmx = +func_bphi(rs)*(y-fr_h)/rs
          -func_bphi(rm)*(y+fr_h)/rm
          -fr_rmom*func_back(rb)*(pow(y+fr_d,3)-3.0*(y+fr_d)*pow(x,2))/pow(rb,3);
    
    }
  } else {
    bmx = 0.0;
  }
  return bmx;
}

static Real func_bmy(const Real x, const Real y)
{
  /*  model field z-component */
  Real rs, rm, rb, r2;
  Real bmy;
  rs = sqrt(pow(x, 2) + pow(y - fr_h, 2));
  rm = sqrt(pow(x, 2) + pow(y + fr_h, 2));
  rb = sqrt(pow(x, 2) + pow(y + fr_d, 2));
  r2 = fr_ri + 0.5*fr_del;
  if (rs > 0.0) {
    if (fr_case == 1) {
    /* dipole case */
    bmy = -func_bphi(rs)*x/rs
          +func_bphi(rm)*x/rm
          -func_bphi(r2)*fr_rmom*fr_d*r2*(2.0*x*(y+fr_d))/pow(rb,4);
    } else {
    /* quadrupole */
    bmy = -func_bphi(rs)*x/rs
          +func_bphi(rm)*x/rm
          -fr_rmom*(func_back(rb))*(pow(x,3) - 3.0*x*pow((y+fr_d),2))/pow(rb,3);
    }
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
        prim(n,k,j,ie+i) = prim(n,k,j,ie);
      }
    }}
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
        b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
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
        prim(IDN,k,js-j,i) = pow(prim(IPR,k,js-j,i), 3./5.);
      }
    }
  }
  // (c) Set face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    Real pbypx;
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pbypx = func_pbypxini(pco->x1f(i), pco->x2f(js-j+1));
          b.x1f(k,(js-j),i) = b.x1f(k,(js-j+1),i)
            - pbypx*pco->dx2v(js-j+1);
        }
      }
    }
    Real az_c, az_p, xc;
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          az_c = func_azini(pco->x1f(i), pco->x2f(js-j));
          az_p = func_azini(pco->x1f(i+1), pco->x2f(js-j));
          b.x2f(k,(js-j),i) = (az_c - az_p)/pco->dx1f(i);
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          b.x3f(k,(js-j),i) = b.x3f(k,js+j-1,i);
        }
      } 
    }
  }
  return;
}

Real func_pbypxini(Real x, Real y) {
  Real pbypx;
  pbypx = (func_bmy(x + 1.0e-9, y)-func_bmy(x - 1.0e-9, y))/2.0e-9;
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
          prim(n,k,je+j,i) = prim(n,k,je,i);
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          b.x1f(k,(je+j),i) = 2.0*b.x1f(k,(je+j-1),i)-b.x1f(k,(je+j-2),i);
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
          b.x3f(k,(je+j  ),i) = 2.0*b.x3f(k,(je+j-1),i)-b.x3f(k,(je+j-2),i);
        }
      }
    }
  }
  return;
}