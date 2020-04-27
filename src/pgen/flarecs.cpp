//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code 
// contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================
//! \file flarecs.c
//  \brief Problem generator for flare current sheet problem.
//
// REFERENCE: For example, see: Shen et al 2011, ApJ. 
// Update:
//  Symmbc: xc_symm = 0.5dx, and xmin = -0.5dx, xmax = 1.0-0.5dx;
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

// functions
static Real func_brini(const Real r, const Real phi);
static Real func_bxini(const Real x1, const Real x2);
static Real func_byini(const Real x1, const Real x2);
static Real func_pcs(const Real x1, const Real x2);
static Real func_azini(Real x, Real y);
static Real func_teini_rel(Real y);

/* Constant parameters */
Real Kb = 1.38e-23, Mp = 1.6726e-27, Mu0 = 1.25663706144e-6, LnA = 30.0;
Real gravitational_const = 6.672e-11; /* (N M^2 kg^-2)*/
Real solarmass = 1.99e+30;            /* (kg) */
Real solarradius = 6.96e+8;           /* (m) */

/* Normalization parameters */
Real Lchar, Bchar, Nechar, Rhochar, Pchar, Timechar, Tchar, Vchar;
Real Gsun_nondim;
Real kappa_nondim_coeff;

// Global parameters to define the initial CS
static int cs_mode, sw_bz;
static Real cs_width, depth;
static Real Te_corona = 2.0e6; // (K)
static Real yc_psipert = 0.25, ly_psipert = 3.5, lx_psipert=2.0, psi_pert=0.01; 

// Boundary conditions
void SymmInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, 
                FaceField &b, Real time, Real dt, 
                int is, int ie, int js, int je, int ks, int ke, int ngh);
void OpenInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, 
                FaceField &b, Real time, Real dt, 
                int is, int ie, int js, int je, int ks, int ke, int ngh);
void OpenOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                 FaceField &b, Real time, Real dt,
                 int is, int ie, int js, int je, int ks, int ke, int ngh);
void LinetiedInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
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
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//==============================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Enroll boundary value function pointers
  //EnrollUserBoundaryFunction(INNER_X1, SymmInnerX1);
  EnrollUserBoundaryFunction(INNER_X1, OpenInnerX1);
  EnrollUserBoundaryFunction(OUTER_X1, OpenOuterX1);
  EnrollUserBoundaryFunction(INNER_X2, LinetiedInnerX2);
  EnrollUserBoundaryFunction(OUTER_X2, OpenOuterX2);
  EnrollUserBoundaryFunction(INNER_X3, ReduceInnerX3);
  EnrollUserBoundaryFunction(OUTER_X3, ReduceOuterX3);
  
  // Enroll AMR
  if(adaptive==true)
      EnrollUserRefinementCondition(RefinementCondition);

  // Diffusion 
  EnrollConductionCoefficient(VariConductivity);

  // Other sources
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

        // Non-dimensional Temperature
        user_out_var(0,k,j,i) = phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);

        // Plasma Beta
        user_out_var(1,k,j,i) = phydro->w(IPR,k,j,i)/pmag;

        // Div(V)
        Real vxl = phydro->w(IM1,k,j,i-1)/phydro->w(IDN,k,j,i-1);
        Real vxr = phydro->w(IM1,k,j,i+1)/phydro->w(IDN,k,j,i+1);
        Real vyd = phydro->w(IM2,k,j-1,i)/phydro->w(IDN,k,j-1,i);
        Real vyu = phydro->w(IM2,k,j+1,i)/phydro->w(IDN,k,j+1,i);
        Real divvc = (vxr - vxl)/pcoord->dx1v(i)
                   + (vyu - vyd)/pcoord->dx2v(j);
        user_out_var(2,k,j,i) = divvc;
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
  /* Send/Recv divv_min
  int irank, nranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  printf("id=%d/%d\n", irank, nranks);
  double *divvmin_buf = NULL;
  divvmin_buf = (double *)malloc(sizeof(double) * nranks);
  
  // Compute local div(v) and find the local minimum
  Real divv_local = 0;
  for (int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        Real vxl = phydro->w(IM1,k,j,i-1)/phydro->w(IDN,k,j,i-1);
        Real vxr = phydro->w(IM1,k,j,i)/phydro->w(IDN,k,j,i);
        Real vyd = phydro->w(IM2,k,j-1,i)/phydro->w(IDN,k,j-1,i);
        Real vyu = phydro->w(IM2,k,j,i)/phydro->w(IDN,k,j,i);
      
        Real divvc = (vxr - vxl)/(pcoord->x1f(i+1) - pcoord->x1f(i))
                 + (vyu - vyd)/(pcoord->x2f(j+1) - pcoord->x2f(j));
        divv_local = std::min(divv_local, divvc);
      }
    }
  }
  // Send all divv to root node and return the minimum value 
  MPI_Gather(&divv_local, 1, MPI_DOUBLE, divvmin_buf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (irank == 0) {
    divv_minimum = -50.0;
    for (int i=0; i<=nranks-1; i++){
      divv_minimum = std::min(divvmin_buf[i], divv_minimum);
    }
    // Set divv_local to divv_minimum
    divv_local = divv_minimum;
  }
  
  //MPI_Bcast(&divv_local, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //printf("id=%d/%d, divv_local=%.3f, divv_min=%.3f\n",
  //  irank, nranks, divv_local, divv_minimum);
  free(divvmin_buf);
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
  Real gm1 = peos->GetGamma() - 1.0;
  
  sw_bz = pin->GetInteger("problem","sw_bz");
  cs_mode = pin->GetInteger("problem","cs_mode");
  
  cs_width = pin->GetReal("problem","cs_width");
  depth = pin->GetReal("problem","depth");

  Lchar = pin->GetOrAddReal("problem", "Lchar", 1.0e8); // (100Mm)
  Bchar = pin->GetOrAddReal("problem", "Bchar", 0.01);  // (100G)
  Nechar = pin->GetOrAddReal("problem", "Nechar", 1.0e16); // (10^10 cm-3)

  Pchar = Bchar * Bchar / Mu0;
  Rhochar = Nechar * Mp;
  Tchar = Pchar / (Nechar * 2. * Kb); /* Total number is 2*ne */
  Vchar = Bchar / sqrt(Mu0 * Rhochar);
  Timechar = Lchar / Vchar;
  kappa_nondim_coeff = pow(Tchar, 3.5)/(Lchar*Pchar*Vchar);

  if (Globals::my_rank == 0) {
    printf("Normalization parameters\n");
    printf("Lchar=%.3e (m), %.3e(Mm)\n", Lchar, Lchar*1.0e-6);
    printf("Vchar=%.3e (m/s), %.3e (km/s)\n", Vchar, Vchar*1.0e-3);
    printf("Timechar=%.3e (s)\n", Timechar);
    printf("Tchar=%.3e (K)\n", Tchar);
    printf("Nechar=%.3e (m^-3). %.3e(kg/m3)\n", Nechar, Rhochar);
    printf("Bchar=%.3e (Tesla)\n", Bchar);
  }
    
  AthenaArray<Real> az;
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  az.NewAthenaArray(nx2,nx1);

  AthenaArray<Real> pgas;
  pgas.NewAthenaArray(nx2,nx1);
  AthenaArray<Real> pcs;
  pcs.NewAthenaArray(nx2,nx1);
  Real p0 = (Te_corona/Tchar);

  // Initialize vector potential
  for (int j=js; j<=je+1; ++j) {
  for (int i=is; i<=ie+1; ++i) {
    az(j,i) = func_azini(pcoord->x1f(i), pcoord->x2f(j));
  }}

  // Initialize gas pressure
  for (int j=js; j<=je+1; ++j) {
  for (int i=is; i<=ie+1; ++i) {
    pcs(j,i) = 0.25*(func_pcs(pcoord->x1f(i), pcoord->x2f(j))
                    +func_pcs(pcoord->x1f(i+1), pcoord->x2f(j))
                    +func_pcs(pcoord->x1f(i), pcoord->x2f(j+1))
                    +func_pcs(pcoord->x1f(i+1), pcoord->x2f(j+1)));
    if (sw_bz == 1) {
      pgas(j,i) = p0;
    } else {
      pgas(j,i) = p0 + pcs(j,i);
    }
  }}

  // Initialize density, momentum, face-centered fields
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real temp = func_teini_rel(0.5*(pcoord->x2f(j)+pcoord->x2f(j+1)));
    phydro->u(IDN,k,j,i) = pgas(j,i)/temp;
    phydro->u(IM1,k,j,i) = 0.0;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;    
  }}}

  // Initialize interface B
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
    pfield->b.x1f(k,j,i) = func_bxini(pcoord->x1f(i), 0.5*(pcoord->x2f(j)+pcoord->x2f(j+1)));
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x2f(k,j,i) = func_byini(0.5*(pcoord->x1f(i)+pcoord->x1f(i+1)), pcoord->x2f(j));
  }}}
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    if (sw_bz == 1) {
      pfield->b.x3f(k,j,i) = sqrt(2.0*pcs(j,i));
    } else {
      pfield->b.x3f(k,j,i) = 0;
    }
  }}}

  // Initialize total energy
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
  pcs.DeleteAthenaArray();
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
// function of the initial br, bx, by, and bz(x, y)
//==============================================================================
static Real func_brini(const Real r, const Real phi)
{
  Real br, b0 = 1.0;
  Real w_phi, rs;
  w_phi = cs_width;
  rs = depth;
  
  /* cs_mode 1: sin type */
  if (cs_mode == 1) {
    if (phi >= w_phi) {
      br = b0 * (rs / r);
    } else if (phi <= -w_phi) {
      br = -b0 * (rs / r); 
    } else {
      br = b0 * (rs / r) * sin(phi * PI * 0.5 / w_phi);
    }
  } else if (cs_mode == 2) {
    /* cs_mode 2: harris sheet */
    br = b0 * (rs / r)*tanh(phi * PI * 0.5 / w_phi);
  } else {
    /* ideal */
    br = -b0 * (rs / r);
  }
  return br;
}

static Real func_bxini(const Real x1, const Real x2)
{
  Real bx0, br0;
  Real r, rs, phi;
  rs = depth;
  r = sqrt(x1 * x1 + (x2 + rs) * (x2 + rs));
  phi = atan(x1 / (x2 + rs));
  br0 = func_brini(r, phi);
  bx0 = br0 * sin(phi);
  
  /* 
  Real w_phi = cs_width;
  Real sx1 = x1*x2;
  Real srsx2 = pow(rs+x2,2);
  bx0 = rs*x1*tanh(0.5*PI*phi/w_phi)/((rs + x2)*sqrt(sx1 + srsx2)*sqrt(sx1/srsx2 + 1));
  */

  /*
  bx0 = 0;
  // Add perturbation
  Real pix, piy, bx1;
  pix = PI * x1 / lx_psipert;
  piy = 2.0 * PI * (x2 - yc_psipert) / ly_psipert;
  bx1 = (2.0 * PI / ly_psipert) * psi_pert * cos(pix) * sin(piy);
  bx0 = bx0 + bx1;
  */

  return bx0;
}

static Real func_byini(const Real x1, const Real x2)
{
  Real by0, br0;
  Real r, rs, phi;
  rs = depth;
  r = sqrt(x1 * x1 + (x2 + rs) * (x2 + rs));
  phi = atan(x1 / (x2 + rs));
  br0 = func_brini(r, phi);
  by0 = br0 * cos(phi);
  
  /*
  Real w_phi = cs_width;
  Real sx1 = x1*x2;
  Real srsx2 = pow(rs+x2,2);
  by0 = rs*tanh(0.5*PI*phi/w_phi)/(sqrt(sx1 + srsx2)*sqrt(sx1/srsx2 + 1));
  */

  /*
  by0 = tanh(x1/cs_width);
  // Add perturbation 
  Real pix, piy, by1;
  pix = PI * x1 / lx_psipert;
  piy = 2.0 * PI * (x2 - yc_psipert) / ly_psipert;
  by1 = (-PI / lx_psipert) * psi_pert * sin(pix) * cos(piy);
  by0 = by0 + by1;
  */

  return by0;
}

//==============================================================================
// Initial Pressure inside the CS
//==============================================================================
static Real func_pcs(const Real x1, const Real x2)
{
  Real pcs, r, phi, br;

  Real y = x2 + depth;
  r = sqrt(x1*x1 + y*y);
  phi = atan(x1/y);
  br = func_brini(r, phi);
  Real r_amb = sqrt(1.0 + y*y);
  Real phi_amb = atan(1.0/y);
  Real br_ambient = func_brini(r_amb, phi_amb);
  pcs = 0.5 * (SQR(br_ambient) - SQR(br));

  /*
  Real bxc = func_bxini(x1, x2);
  Real byc = func_byini(x1, x2);
  Real pmag_c = 0.5*(bxc*bxc + byc*byc);
  Real pmag_amb = 0.5;
  pcs = pmag_amb - pmag_c;
  */

  return pcs;
}

//==============================================================================
// function of the initial t(x, y)
//==============================================================================
Real func_teini_rel(Real y) {
  Real t;
  Real h = 0.045;
  Real w = 0.01;
  Real t1, t2, t_bott, t_coro, t_chro;
  t_coro = 1.0; // This is a scaled T, not the non-dimensional T.
  t_chro = 0.005;
  t1 = 0.5 * (t_coro - t_chro);
  t2 = 0.5 * (t_coro + t_chro);
  t = t1 * tanh((y - h) / w) + t2;
  // Use T0 = Te_corona K.
  t = t*(Te_corona/Tchar);
  return t;
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
        prim(n,k,j,is-i) = prim(n,k,j,is+i);
      }
    }}
  }

  // Set velocity Vx
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        prim(IVX,k,j,is-i) = -prim(IVX,k,j,is+i);
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        b.x1f(k,j,(is-i)) = b.x1f(k,j,is+i+1);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
      for (int i=1; i<=ngh; ++i) {
        b.x2f(k,j,(is-i)) = -b.x2f(k,j,is+i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        b.x3f(k,j,(is-i)) = b.x3f(k,j,is+i);
      }
    }}
  }
  
  return;
}

//==============================================================================
// OpenInnerX1 boundary condution
//==============================================================================
void OpenInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, 
                FaceField &b, Real time, Real dt, 
                int is, int ie, int js, int je, int ks, int ke, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        prim(n,k,j,is-i) = prim(n,k,j,is);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(is-i)) = 2.0*b.x2f(k,j,is-i+1) - b.x2f(k,j,is-i+2);
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=1; i<=ngh; ++i) {
          Real xl = pco->x1f(is-i);
          Real xr = pco->x1f(is-i+1);
          Real yd = pco->x2f(j);
          Real yu = pco->x2f(j+1);
          b.x1f(k,j,(is-i)) = b.x1f(k,j,(is-i+1))
            +(xr-xl)/(yu-yd)*(b.x2f(k,(j+1),(is-i)) - b.x2f(k,j,(is-i)));
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
        }
      }
    }
  }
  return;
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
          Real xl = pco->x1f(ie+i);
          Real xr = pco->x1f(ie+i+1);
          Real yd = pco->x2f(j);
          Real yu = pco->x2f(j+1);
          
          b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+i))
            -(xr-xl)/(yu-yd)*(b.x2f(k,(j+1),(ie+i)) - b.x2f(k,j,(ie+i)));
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
        Real temp = func_teini_rel(0.5*(pco->x2f(js-j)+pco->x2f(js-j+1)));
        prim(IDN,k,js-j,i) = prim(IPR,k,js,i)/temp;
        //Real temp = prim(IPR,k,js+j-1,i)/prim(IDN,k,js+j-1,i);
        //prim(IPR,k,js,i) = prim(IDN,k,js-j,i)*temp;
      }
    }
  }
  // (c) Set face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) { 
          b.x2f(k,(js-j),i) = func_byini(0.5*(pco->x1f(i)+pco->x1f(i+1)),
                                         pco->x2f(js-j));
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          Real xl = 0.5*(pco->x1f(i)+pco->x1f(i-1));
          Real xr = 0.5*(pco->x1f(i)+pco->x1f(i+1));
          Real xc = pco->x1f(i);
          Real yd = 0.5*(pco->x2f(js-j)+pco->x2f(js-j+1));
          Real yu = 0.5*(pco->x2f(js-j+1)+pco->x2f(js-j+2));
          Real yc = pco->x2f(js-j+1);
          Real ddx = 0.5*(pco->x1f(i) - pco->x1f(i-1));
           
          //Real pbypx = (b.x2f(k,(js-j+1),i)-b.x2f(k,(js-j+1),i-1))/(xr - xl);
          //Real pbypx = (func_byini(xc+ddx,yc)-func_byini(xc-ddx,yc))/(2.0*ddx);
          //b.x1f(k,(js-j),i) = b.x1f(k,(js-j+1),i) - pbypx*(yu - yd);
          b.x1f(k,(js-j),i) = func_bxini(pco->x1f(i),
                                         0.5*(pco->x2f(js-j)+pco->x2f(js-j+1)));
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

  // Inflows
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=is; i<=ie; ++i) {
        prim(IVY,k,je+j,i) = std::max(0.0, prim(IVY,k,je,i));
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          b.x1f(k,je+j,i) = b.x1f(k,je,i);
        }
      }
    }

    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          //Real xl = pco->x1f(i);
          //Real xr = pco->x1f(i+1);
          //Real yd = pco->x2f(je+j+1);
          //Real yu = pco->x2f(je+j);
          //b.x2f(k,(je+j+1),i) = b.x2f(k,(je+j),i)
          //-(yu-yd)/(xr-xl)*(b.x1f(k,(je+j),i+1)-b.x1f(k,(je+j),i));
          b.x2f(k,(je+j+1),i) = b.x2f(k,(je+1),i);
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          b.x3f(k,je+j,i) = b.x3f(k,je,i);
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
        prim(IVZ,ks-k,j,i) = std::min(0.0, prim(IVZ,ks-k,j,i));
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
        prim(IVZ,ke+k,j,i) = std::max(0.0,prim(IVZ,ke+k,j,i));
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
// AMR refinement condition functions
//------------------------------------------------------------------------------
int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  
  /* Div(v) xy-components and Velocity y-*/
  Real divv_min = 0;
  Real vy_min = 0;
  for(int k=pmb->ks; k<=pmb->ke; k++) {
    for(int j=pmb->js; j<=pmb->je; j++) {
      for(int i=pmb->is; i<=pmb->ie; i++) {
        
        Real vxl = w(IM1,k,j,i-1)/w(IDN,k,j,i-1);
        Real vxc = w(IM1,k,j,i)/w(IDN,k,j,i);
        //Real vxr = w(IM1,k,j,i+1)/w(IDN,k,j,i+1);

        Real xl = 0.5*(pmb->pcoord->x1f(i-1) + pmb->pcoord->x1f(i));
        Real xc = 0.5*(pmb->pcoord->x1f(i) + pmb->pcoord->x1f(i+1));
        //Real xr = 0.5*(pmb->pcoord->x1f(i+1) + pmb->pcoord->x1f(i+2));

        Real vyd = w(IM2,k,j-1,i)/w(IDN,k,j-1,i);
        Real vyc = w(IM2,k,j,i)/w(IDN,k,j,i);
        //Real vyu = w(IM2,k,j+1,i)/w(IDN,k,j+1,i);

        Real yd = 0.5*(pmb->pcoord->x2f(j-1) + pmb->pcoord->x2f(j));
        Real yc = 0.5*(pmb->pcoord->x2f(j) + pmb->pcoord->x2f(j+1));
        //Real yu = 0.5*(pmb->pcoord->x2f(j+1) + pmb->pcoord->x2f(j+2));
        
        Real divvc = (vxc - vxl)/(xc - xl) + (vyc - vyd)/(yc - yd);
        divv_min = std::min(divv_min, divvc);
        vy_min = std::min(vyc, vy_min);
      }
    }
  }
  if ((divv_min <= -10.0) && (vy_min <= -0.2)) {
    return 1;
  } else if ((divv_min >= -1.0) || (vy_min > -0.05)) {
    return -1;
  } else {
    return 0;
  }
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
        
        // Limit Te at NUN points
        if (te(k,j,i) >= 5.0e7/Tchar) {
          te(k,j,i) = 5.0e7/Tchar;
        }

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
        
      }
    }
  }

  // Release memory
  te.DeleteAthenaArray();

  return;
}

//-----------------------------------------------------------------------------
// Cooling and Heating sources
//-----------------------------------------------------------------------------
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
        Real temp_const = func_teini_rel(0.5*(pmb->pcoord->x2f(j)+pmb->pcoord->x2f(j+1)));
        Real pgas;
        Real p0 = (Te_corona/Tchar);
        if (sw_bz == 1) {
          pgas = p0;
        } else {
          Real pcs  = 0.25*(func_pcs(pmb->pcoord->x1f(i), pmb->pcoord->x2f(j))
                    +func_pcs(pmb->pcoord->x1f(i+1), pmb->pcoord->x2f(j))
                    +func_pcs(pmb->pcoord->x1f(i), pmb->pcoord->x2f(j+1))
                    +func_pcs(pmb->pcoord->x1f(i+1), pmb->pcoord->x2f(j+1)));
          pgas = p0 + pcs;
        }
        Real rho_const = pgas/temp_const;
        Real heatrate = (rho_const*Nechar)*(rho_const*Nechar)*Qt(temp_const*Tchar);
        heatrate = heatrate*Timechar/Pchar;

        /* Heating terms should not directly affect CS temperature
        if (temp > Te_corona/Tchar) {
          Real Te_range = 0.2e6;
          Real Te_frac = (temp*Tchar - Te_corona)/Te_range;
          Real factor = exp(-Te_frac*Te_frac);
          heatrate = factor*heatrate;
        }
        */

        // At the bottom, H always == cool
        if (temp < Te_corona/Tchar) {
          heatrate = coolrate;
        }

        // Total heating rate should be <= cooling rate
        if (heatrate > coolrate) {
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

  /* to SI unit: W m^3 */
  q = q * 1.0e-13;

  return q;
}

