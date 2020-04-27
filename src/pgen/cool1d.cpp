//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code 
// contributors.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================
//! \file cool1d.cpp
//  \brief Problem generator for cool 1d test.
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

// Global array
AthenaArray<Real> az;
AthenaArray<Real> rhof;
AthenaArray<Real> pf;

// Global parameters to define the background gas
static Real scale_bgdens, scale_lowtcore;
static Real Te_corona, Te_bottom;
static Real y_base, p_base, rho_base, te_base;
static Real height_bottom, width_bottom;

// Temperature dependent kappa
void VariConductivity(HydroDiffusion *phdif, MeshBlock *pmb,
              const AthenaArray<Real> &w, const AthenaArray<Real> &bc,
              int il, int iu, int jl, int ju, int kl, int ku);

// Cooling and Heating terms
static Real Qt(const Real T);
void coolingheating_source(MeshBlock *pmb, const Real time, const Real dt, 
                    const AthenaArray<Real> &prim,
                    const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

//==============================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also
//  be used to initialize variables which are global to (and therefore can be 
//  passed to) other functions in this file.  Called in Mesh constructor.
//==============================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  
  // Diffusion
  EnrollConductionCoefficient(VariConductivity);

  // Static gravity source in momentum equations
  EnrollUserExplicitSourceFunction(coolingheating_source);
  
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
  Te_bottom = pin->GetOrAddReal("problem", "Te_bottom", 1.0e4); // (unit: K)
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
  rhof.NewAthenaArray(nx2+1,nx1+1);
  pf.NewAthenaArray(nx2+1,nx1+1);
  az.NewAthenaArray(nx2+1,nx1+1);

  // Initialize the pressure_fluxrop, density and pressure at the cell interface
  for (int j=js; j<=je+1; ++j) {
  for (int i=is; i<=ie+1; ++i) {
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
  az.DeleteAthenaArray();
  return;
}

//==============================================================================
// function: initial p, te, and rho(x, y)
//==============================================================================
static Real func_pini(const Real x, const Real y) {
  return p_base;
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
        
        /* increase kappa if T<500,000K 
        Real T_lowlevel = 5.0e5;
        Real T = te(k,j,i)*Tchar;
        if (T <= T_lowlevel) {
          Real factor = pow(T_lowlevel/T, 3.0);
          phdif->kappa(1,k,j,i) *= factor;
        }
        */
      }
    }
  }

  // Release memory
  te.DeleteAthenaArray();

  return;
}

void coolingheating_source(MeshBlock *pmb, const Real time, const Real dt, 
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

        // Heating: H_const
        Real temp_const = func_teini(pmb->pcoord->x1f(i), pmb->pcoord->x2f(j));
        Real rho_const = func_rhoini(pmb->pcoord->x1f(i), pmb->pcoord->x2f(j));
        Real heatrate = ne*(rho_const*Nechar)*Qt(temp_const*Tchar);
        heatrate = heatrate*Timechar/Pchar;
        // At the bottom, H == cool
        if (temp < Te_corona/Tchar) {
          heatrate = coolrate;
        }

        // Return EN
        cons(IEN,k,j,i) -= dt*(coolrate-heatrate);
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

  /* Decrease Q(T) if T<500,000K 
  Real T_lowlevel = 5.0e5;
  if (T <= T_lowlevel) {
    factor = pow(T_lowlevel/T, 3.0);
    q /= factor;
  }
  */

  /* to SI unit: W m^3 */
  q = q * 1.0e-13;

  return q;
}