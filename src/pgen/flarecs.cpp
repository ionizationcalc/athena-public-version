//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file flarecs.c
//  \brief Problem generator for flare current sheet problem.
//
// REFERENCE: For example, see: Shen et al 2011, ApJ. 
//========================================================================================

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
Real func_azini(Real x, Real y);
Real func_pbypxini(Real x, Real y);

// Global parameters to define the initial CS
static Real cs_width, beta0;

// Boundary conditions
void LinetiedInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int is, int ie, int js, int je, int ks, int ke, int ngh);
void OpenOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, 
                FaceField &b, Real time, Real dt, 
                int is, int ie, int js, int je, int ks, int ke, int ngh);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Enroll boundary value function pointers
  EnrollUserBoundaryFunction(INNER_X2, LinetiedInnerX2);
  EnrollUserBoundaryFunction(OUTER_X2, OpenOuterX2);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the flarecs.  The initial conditions are
//  constructed assuming the domain extends over [-1.0x1.0, 0.0x2.0].
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // Read input arguments
  Real gm1 = peos->GetGamma() - 1.0;
  cs_width = pin->GetReal("problem","cs_width");
  beta0 = pin->GetReal("problem","beta0");

  AthenaArray<Real> az;
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  az.NewAthenaArray(nx2,nx1);

  AthenaArray<Real> pgas;
  pgas.NewAthenaArray(nx2,nx1);

  Real B0 = 1.0;
  Real d0 = 1.0;
  Real v0 = 1.0;
  Real p0 = 0.05;

  // Initialize vector potential
  for (int j=js; j<=je+1; ++j) {
  for (int i=is; i<=ie+1; ++i) {
    az(j,i) = func_azini(pcoord->x1f(i), pcoord->x2f(j));
  }}

  // Initialize gas pressure
  Real bxc, byc;
  for (int j=js; j<=je+1; ++j) {
  for (int i=is; i<=ie+1; ++i) {
    bxc = 0.5*((az(j+1,i) - az(j,i))/pcoord->dx2f(j) + (az(j+1,i+1) - az(j,i+1))/pcoord->dx2f(j));
    byc = 0.5*((az(j,i) - az(j,i+1))/pcoord->dx1f(i) + (az(j+1,i) - az(j+1,i+1))/pcoord->dx1f(i));
    // Assuming p_total = p0 + 0.5*1^2 and bz == 0
    pgas(j,i) = (p0 + 0.5) - 0.5*(bxc*bxc + byc*byc);
  }}

  // Initialize density, momentum, face-centered fields
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    phydro->u(IDN,k,j,i) = pgas(j,i)/(beta0/2.0);
    phydro->u(IM1,k,j,i) = 0.0;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;
  }}}

  // initialize interface B
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
    pfield->b.x1f(k,j,i) = (az(j+1,i) - az(j,i))/pcoord->dx2f(j);
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x2f(k,j,i) = (az(j,i) - az(j,i+1))/pcoord->dx1f(i);
  }}}
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x3f(k,j,i) = 0.0;
  }}}

  // initialize total energy
  if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      phydro->u(IEN,k,j,i) = pgas(j,i)/gm1 +
          0.5*(SQR(0.5*(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))) +
               SQR(0.5*(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i))) +
               SQR(0.5*(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i)))) + (0.5)*
          (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i))
           + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
    }}}
  }

  az.DeleteAthenaArray();
  pgas.DeleteAthenaArray();
  return;
}

//==============================================================================
// function of the initial az(x, y)
//==============================================================================
Real func_azini(Real x, Real y) {
  Real az0 = 0, az;
  az = -cs_width * log(cosh(x / cs_width)) + az0;
  return az;
}

//==============================================================================
// function: pbypx
//==============================================================================
Real func_pbypxini(Real x, Real y) {
  Real pbypx;
  pbypx = (1.0/cs_width)*(1.0 - pow(tanh(x/cs_width), 2));
  return pbypx;
}

//==============================================================================
// LinetiedInnerX2 boundary condition
//==============================================================================
void LinetiedInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
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
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=is; i<=ie; ++i) {
        prim(IVX,k,js-j,i) = 0;
        prim(IVY,k,js-j,i) = 0;
        prim(IVZ,k,js-j,i) = 0;
        //prim(IDN,k,js-j,i) = prim(IPR,k,js,i)/(beta0/2.0);
      }
    }
  }
  // (c) Set face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    Real pbypx;
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pbypx = func_pbypxini(pco->x1f(i), pco->x2f(j));
          b.x1f(k,(js-j),i) = b.x1f(k,(js-j+1),i) - pbypx*pco->dx2f(js-j+1);
        }
      }
    }
    Real az_c, az_p;
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
          b.x3f(k,(je+j  ),i) = b.x3f(k,(je  ),i);
        }
      }
    }
  }
  return;
}