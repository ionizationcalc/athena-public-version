//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code
// contributors.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================
//! \file conduction_2dtest.c
//  \brief Problem generator for 2d conduction tests.
//
// REFERENCE:
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

// Global parameters to define the initial fluxrope
static Real beta0, temperature0;
static Real p_ambient=0.6, te_ambient=0.6, rho_ambient=1.0;

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
  /* remove velocity and reset pressure to the unit.
  Real xc, yc, r;
  Real tec, pc, ek, eb;
  //Real gm1 = peos->GetGamma() - 1.0;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        //xc = 0.5*(pcoord->x1f(i) + pcoord->x1f(i+1));
        //yc = 0.5*(pcoord->x2f(j) + pcoord->x2f(j+1));
        //r = sqrt(xc*xc + yc*yc);
        ek = 0.5*(SQR(phydro->u(IM1,k,j,i))
                + SQR(phydro->u(IM2,k,j,i))
                + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = phydro->u(IEN,k,j,i) - ek;
      }}}
   */
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
  
  // Define the local temporary array: az & pgas
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  AthenaArray<Real> az;
  az.NewAthenaArray(nx2,nx1);
  
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
        pfield->b.x1f(k,j,i) = 0.1;
      }}}
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie; i++) {
        pfield->b.x2f(k,j,i) = 0.1;
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
  return p_ambient;
}

Real func_teini(Real x, Real y) {
  Real te, te_max;
  Real r, w;
  w = 0.05;
  r = sqrt(x*x + y*y);
  te_max = 20.0*te_ambient;
  
  if (r <= w) {
    te = te_max;
  } else {
    te = te_ambient;
  }
  
  //te = te_ambient + te_max*exp(-r*r/(2.0*w*w));
  return te;
}

Real func_rhoini(Real x, Real y) {
  return func_pini(x,y)/func_teini(x,y);
}
