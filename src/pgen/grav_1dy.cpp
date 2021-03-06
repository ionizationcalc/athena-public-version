//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code 
// contributors.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================
//! \file grav_1d.cpp
//  \brief Problem generator for the fluxrope based on Lin-Forbes model.
//
// REFERENCE:
//  
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
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

// Normalization parameters
Real Lchar, Bchar, Nechar, Rhochar, Pchar, Timechar, Tchar, Vchar;
Real Gsun_nondim;
Real kappa_nondim_coeff;

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


Real func_integ_py(Real y, Real xfix);

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
AthenaArray<Real> rhof;
AthenaArray<Real> pf;

// Global parameters to define the background gas
static Real scale_bgdens, scale_lowtcore;
static Real Te_corona, Te_bottom;
static Real y_base, p_base, rho_base, te_base;
static Real beta_floor;
static Real height_bottom, width_bottom;

// made global to share with BC functions
static Real gsun_y(const Real x2);         // Gsun(r)

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

// Temperature dependent kappa
void VariConductivity(HydroDiffusion *phdif, MeshBlock *pmb,
              const AthenaArray<Real> &w, const AthenaArray<Real> &bc,
              int il, int iu, int jl, int ju, int kl, int ku);

// Cooling and Heating terms
static Real Qt(const Real T);
void cooling_source(MeshBlock *pmb, const Real time, const Real dt, 
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
  //EnrollUserBoundaryFunction(INNER_X1, SymmInnerX1);
  //EnrollUserBoundaryFunction(OUTER_X1, OpenOuterX1);
  EnrollUserBoundaryFunction(INNER_X2, LintInnerX2);
  EnrollUserBoundaryFunction(OUTER_X2, OpenOuterX2);
  //EnrollUserBoundaryFunction(INNER_X3, ReduceInnerX3);
  //EnrollUserBoundaryFunction(OUTER_X3, ReduceOuterX3);
  
  // Diffusion
  //EnrollFieldDiffusivity(VariDiffusivity);
  EnrollConductionCoefficient(VariConductivity);

  if(adaptive==true)
      EnrollUserRefinementCondition(RefinementCondition);

  // Static gravity source in momentum equations
  EnrollUserExplicitSourceFunction(static_grav_source);
  EnrollUserExplicitSourceFunction(cooling_source);
  
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

  kappa_nondim_coeff = pow(Tchar, 3.5)/(Lchar*Pchar*Vchar);
  
  // Initialize background gas
  scale_bgdens = pin->GetOrAddReal("problem", "scale_bgdens", 1.0); //
  Te_corona = pin->GetOrAddReal("problem", "Te_corona", 1.0e6); // (unit: K)
  Te_bottom = pin->GetOrAddReal("problem", "Te_bottom", 1.0e6); // (unit: K)
  height_bottom = pin->GetOrAddReal("problem", "height_bottom", 0.1); // (non-dim)
  width_bottom = pin->GetOrAddReal("problem", "width_bottom", 0.05); // (non-dim)
  
  // Define the base density and temperature
  y_base = 0;
  rho_base = scale_bgdens;
  te_base = Te_corona / Tchar;
  p_base = rho_base * te_base;

  //if (my_rank == 0) {
    printf("---------------\n");
    printf("rho_base=%f, Te_base=%f, p_base=%f\n", rho_base, te_base, p_base);
    printf("Lchar=%.3e, Vchar=%.3e, Timechar=%.3e, Tchar=%.3e, Nechar=%.3e, Bchar=%.3e\n", Lchar, Vchar, Timechar, Tchar, Nechar, Bchar);
    printf("---------------\n");
  //}

  // Define the local temporary array: az & pgas
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  Real xc, yc;
  AthenaArray<Real> pbzf;
  //AthenaArray<Real> rhof;
  //AthenaArray<Real> pf;
  pbzf.NewAthenaArray(nx2+1,nx1+1);
  rhof.NewAthenaArray(nx2+1,nx1+1);
  pf.NewAthenaArray(nx2+1,nx1+1);
  az.NewAthenaArray(nx2+1,nx1+1);

  // Initialize the pressure_fluxrop, density and pressure at the cell interface
  for (int j=js; j<=je+1; ++j) {
  for (int i=is; i<=ie+1; ++i) {
    // Pressure for Bz
    pbzf(j,i) = 0;
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
    pfield->b.x1f(k,j,i) = 0;
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x2f(k,j,i) = 1.0;
  }}}
  
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
      pfield->b.x3f(k,j,i) = 0;
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
  //rhof.DeleteAthenaArray();
  //pf.DeleteAthenaArray();
  az.DeleteAthenaArray();
  return;
}

//==============================================================================
// function: initial p, te, and rho(x, y)
//==============================================================================
static Real func_pini(const Real x, const Real y) {
  Real p;
  Real y0 = 1.0;
  Real p0 = func_p_isoteconst(y0, y_base, p_base, te_base);
  p = presini_integral(y, y0, p0);
  return p;
}

static Real func_rhoini(const Real x, const Real y) {
  Real rho = func_pini(x, y)/func_teini(x, y);
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
                          y0, y, 1.0e-10, 10000);
  pgas = p0 * exp(Iout);
  return pgas;
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

        /*
        // case 3: isothermal 
        Real y0 = 0.5*(pco->x2f(js) + pco->x2f(js+1));
        Real p0 = prim(IPR,k,js,i);
        Real te0 = p0/prim(IDN,k,js,i);
        prim(IPR,k,js-j,i) = func_p_isoteconst(yc, y0, p0, te0);
        prim(IDN,k,js-j,i) = prim(IPR,k,js-j,i)/te0;
        */
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
          b.x2f(k,(js-j),i) = b.x2f(k,js,i);
        }
      }
    }
    Real pbypx;
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          b.x1f(k,(js-j),i) = b.x1f(k,js,i);
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
          //Real grav_acc = -gsun_y(yc);
          //prim(IPR,k,je+j,i) = prim(IPR,k,je-j+1,i)
          //   + prim(IDN,k,je-j+1,i)*grav_acc*(2*j-1)*pco->dx2f(j);
          
          Real y0 = 0.5*(pco->x2f(je) + pco->x2f(je+1));
          Real p0 = prim(IPR,k,je,i);
          Real te0 = p0/prim(IDN,k,je,i);
          prim(IPR,k,je+j,i) = func_p_isoteconst(yc, y0, p0, te0);
          prim(IDN,k,je+j,i) = prim(IPR,k,je+j,i)/te0;
          
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
  if(maxeps > 0.5) pflag = 1;
  if(maxeps < 0.05) pflag = -1;
  
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
  if(jz_max >= 10.0) jflag = 1;
  if(jz_max <= 1.0) jflag = -1;

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

//------------------------------------------------------------------------------
// Conductivity from variable coefficients
//------------------------------------------------------------------------------
void VariConductivity(HydroDiffusion *phdif, MeshBlock *pmb,
              const AthenaArray<Real> &w, const AthenaArray<Real> &bc,
              int il, int iu, int jl, int ju, int kl, int ku) {
  int i,j,k;

  // Compute temperauture at the cell center.
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);
  AthenaArray<Real> te;
  te.NewAthenaArray(ncells3,ncells2,ncells1);
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        te(k,j,i) = w(IPR,k,j,i)/w(IDN,k,j,i);
      }
    }
  }
  
  // kappa_iso
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        phdif->kappa(0,k,j,i) = 0;
      }
    }
  }

  // kappa_aniso
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        phdif->kappa(1,k,j,i) = 1.0e-11*pow(te(k,j,i), 2.5)*kappa_nondim_coeff;
        
        /* increase kappa if T<500,000K */
        Real T_lowlevel = 5.0e5;
        Real T = te(k,j,i)*Tchar;
        if (T <= T_lowlevel) {
          Real factor = pow(T_lowlevel/T, 3.0);
          phdif->kappa(1,k,j,i) *= factor;
        }
      }
    }
  }

  // Release memory
  te.DeleteAthenaArray();

  return;
}

void cooling_source(MeshBlock *pmb, const Real time, const Real dt, 
                    const AthenaArray<Real> &prim,
                    const AthenaArray<Real> &bcc, AthenaArray<Real> &cons) {
  Real g = pmb->peos->GetGamma();
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        // Cooling: rho*2 Q(T)
        Real temp = prim(IPR,k,j,i)/prim(IDN,k,j,i);
        Real ne = prim(IDN,k,j,i)*Nechar;
        Real coolrate = ne*ne*Qt(temp*Tchar);
        coolrate = coolrate*Timechar/Pchar;
        cons(IEN,k,j,i) -= dt*coolrate;
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*  \brief Calculate Qt
 *  piecewise linear approximation (Klimcuk et al. 2008)*/
/* T: dimensional variable, SI unit */
static Real Qt(const Real T) {
  Real q;
  Real factor;
  Real logt = log10(T);

  /* first in cgs: ergs sec^-1 cm^3 */
  if (logt <= 4.97) {
    q = 1.09e-31 * (pow(T, 2));
  } else if (logt <= 5.67) {
    q = 8.87e-17 * (pow(T, -1.0));
  } else if (logt <= 6.18) {
    q = 1.90e-22;
  } else if (logt <= 6.55) {
    q = 3.54e-13 * (pow(T, -3. / 2.));
  } else if (logt <= 6.90) {
    q = 3.46e-25 * (pow(T, 1. / 3.));
  } else if (logt <= 7.63) {
    q = 5.49e-16 * (pow(T, -1.0));
  } else {
    q = 1.96e-27 * (pow(T, 0.5));
  }

  /* Decrease Q(T) if T<500,000K */
  Real T_lowlevel = 5.0e5;
  if (T <= T_lowlevel) {
    factor = pow(T_lowlevel/T, 3.0);
    q /= factor;
  }

  /* to SI unit: W m^3 */
  q = q * 1.0e-13;

  return q;
}