//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code
// contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================

// Athena++ headers
#include "hydro_diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../eos/eos.hpp"

#include<algorithm>
using namespace std;
static Real limiter2(const Real A, const Real B);
static Real limiter4(const Real A, const Real B, const Real C, const Real D);
static Real vanleer (const Real A, const Real B);
static Real minmod  (const Real A, const Real B);

//------------------------------------------------------------------------------
// limiter2 and limiter4: call slope limiters to preserve monotonicity
static Real limiter2(const Real A, const Real B)
{
  /* van Leer slope limiter */
  return vanleer(A,B);
  
  /* monotonized central (MC) limiter */
  /* return minmod(2.0*minmod(A,B),0.5*(A+B)); */
}

static Real limiter4(const Real A, const Real B, const Real C, const Real D)
{
  return limiter2(limiter2(A,B),limiter2(C,D));
}

//------------------------------------------------------------------------------
// vanleer: van Leer slope limiter
static Real vanleer(const Real A, const Real B)
{
  if (A*B > 0) {
    return 2.0*A*B/(A+B);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// minmod: minmod slope limiter

static Real minmod(const Real A, const Real B)
{
  if (A*B > 0) {
    if (A > 0) {
      return min(A,B);
    } else {
      return max(A,B);
    }
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Calculate isotropic thermal conduction

void HydroDiffusion::ThermalFlux_iso(const AthenaArray<Real> &prim,
              const AthenaArray<Real> &cons, AthenaArray<Real> *cndflx) {
  AthenaArray<Real> &x1flux=cndflx[X1DIR];
  AthenaArray<Real> &x2flux=cndflx[X2DIR];
  AthenaArray<Real> &x3flux=cndflx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real kappaf, denf, dTdx, dTdy, dTdz;

  // i-direction
  jl=js, ju=je, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) {
      if(pmb_->block_size.nx3 == 1) // 2D
        jl=js-1, ju=je+1, kl=ks, ku=ke;
      else // 3D
        jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        kappaf = 0.5*(kappa(ISO,k,j,i)+kappa(ISO,k,j,i-1));
        denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j,i-1));
        dTdx = (prim(IPR,k,j,i)/prim(IDN,k,j,i) - prim(IPR,k,j,i-1)/
                prim(IDN,k,j,i-1))/pco_->dx1v(i-1);
        x1flux(k,j,i) -= kappaf*denf*dTdx;
      }
    }
  }

  // j-direction
  il=is, iu=ie, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx3 == 1) // 2D
      il=is-1, iu=ie+1, kl=ks, ku=ke;
    else // 3D
      il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
  }
  if(pmb_->block_size.nx2 > 1) { //2D or 3D
    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          kappaf = 0.5*(kappa(ISO,k,j,i)+kappa(ISO,k,j-1,i));
          denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j-1,i));
          dTdy = (prim(IPR,k,j,i)/prim(IDN,k,j,i)-prim(IPR,k,j-1,i)/
                    prim(IDN,k,j-1,i))/pco_->h2v(i)/pco_->dx2v(j-1);
          x2flux(k,j,i) -= kappaf*denf*dTdy;
        }
      }
    }
  } // zero flux for 1D

  // k-direction
  il=is, iu=ie, jl=js, ju=je;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) // 2D or 3D
      il=is-1, iu=ie+1, jl=js-1, ju=je+1;
    else // 1D
      il=is-1, iu=ie+1;
  }
  if(pmb_->block_size.nx3 > 1) { //3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          kappaf = 0.5*(kappa(ISO,k,j,i)+kappa(ISO,k-1,j,i));
          denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k-1,j,i));
          dTdz = (prim(IPR,k,j,i)/prim(IDN,k,j,i)-prim(IPR,k-1,j,i)/
                   prim(IDN,k-1,j,i))/pco_->dx3v(k-1)/pco_->h31v(i)/pco_->h32v(j);
          x3flux(k,j,i) -= kappaf*denf*dTdz;
        }
      }
    }
  } // zero flux for 1D/2D

  return;
}


//------------------------------------------------------------------------------
// Calculate anisotropic thermal conduction

void HydroDiffusion::ThermalFlux_aniso(const AthenaArray<Real> &prim,
              const AthenaArray<Real> &b, AthenaArray<Real> *cndflx) {
  AthenaArray<Real> &x1flux=cndflx[X1DIR];
  AthenaArray<Real> &x2flux=cndflx[X2DIR];
  AthenaArray<Real> &x3flux=cndflx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real kappaf, denf, kd, dTdx, dTdy, dTdz;
  Real Bx,By,Bz,B02, bDotGradT;

  // Problem must be at least 2D
  if (pmb_->block_size.nx2 == 1) return;

  // Compute temperature at cell centers.
  int ncells1 = pmb_->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb_->block_size.nx2 > 1) ncells2 = pmb_->block_size.nx2 + 2*(NGHOST);
  if (pmb_->block_size.nx3 > 1) ncells3 = pmb_->block_size.nx3 + 2*(NGHOST);

  AthenaArray<Real> te;
  te.NewAthenaArray(ncells3,ncells2,ncells1);
  
  if(pmb_->block_size.nx3 == 1) // 2D
    jl=js-1, ju=je+1, kl=ks, ku=ke;
  else // 3D
    jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
  
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
  #pragma omp simd
      for (int i=is-1; i<=ie+1; ++i) {
        te(k,j,i) = prim(IPR,k,j,i)/prim(IDN,k,j,i);
      }
    }
  }

  // i-direction
  jl=js, ju=je, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) {
      if(pmb_->block_size.nx3 == 1) // 2D
        jl=js-1, ju=je+1, kl=ks, ku=ke;
      else // 3D
        jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
  #pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        kappaf = 0.5*(kappa(ANI,k,j,i)+kappa(ANI,k,j,i-1));
        denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j,i-1));
        dTdx = (te(k,j,i) - te(k,j,i-1))/pco_->dx1v(i-1);
        
        // Monotonized temperature difference dT/dy
        dTdy = limiter4(te(k,j+1,i  ) - te(k,j  ,i  ),
                        te(k,j  ,i  ) - te(k,j-1,i  ),
                        te(k,j+1,i-1) - te(k,j  ,i-1),
                        te(k,j  ,i-1) - te(k,j-1,i-1));
        dTdy /= (pco_->h2v(i)*pco_->dx2v(j-1));
        
        // Monotonized temperature difference dT/dz, 3D problem ONLY
        if(pmb_->block_size.nx3 > 1) {
          dTdz = limiter4(te(k+1,j,i  ) - te(k,  j,i  ),
                          te(k  ,j,i  ) - te(k-1,j,i  ),
                          te(k+1,j,i-1) - te(k  ,j,i-1),
                          te(k,  j,i-1) - te(k-1,j,i-1));
          dTdz /= (pco_->dx3v(k-1)*pco_->h31v(i)*pco_->h32v(j));
        }
        
        // Add flux at x1-interface, 2D PROBLEM
        if(pmb_->block_size.nx3 == 1) {
          Bx = b(0,k,j,i);
          By = 0.5*(b(1,k,j,i-1) + b(1,k,j,i));
          B02 = SQR(Bx) + SQR(By);
          B02 = max(B02,TINY_NUMBER); /* limit in case B=0 */
          bDotGradT = Bx*dTdx + By*dTdy;
          kd = kappaf*denf;
          x1flux(k,j,i) -= kd*(Bx*bDotGradT)/B02;

        // Add flux at x1-interface, 3D PROBLEM
        } else {
          Bx = b(0,k,j,i);
          By = 0.5*(b(1,k,j,i-1) + b(1,k,j,i));
          Bz = 0.5*(b(2,k,j,i-1) + b(2,k,j,i));
          B02 = SQR(Bx) + SQR(By) + SQR(Bz);
          B02 = max(B02,TINY_NUMBER); // limit in case B=0
          bDotGradT = Bx*dTdx + By*dTdy + Bz*dTdz;
          kd = kappaf*denf;
          x1flux(k,j,i) -= kd*(Bx*bDotGradT)/B02;
        }
      }
    }
  }
  
  // j-direction
  il=is, iu=ie, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx3 == 1) // 2D
      il=is-1, iu=ie+1, kl=ks, ku=ke;
    else // 3D
      il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
  }
  if(pmb_->block_size.nx2 > 1) { //2D or 3D
    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
  #pragma omp simd
        for (int i=il; i<=iu; ++i) {
          kappaf = 0.5*(kappa(ANI,k,j,i)+kappa(ANI,k,j-1,i));
          denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j-1,i));
          dTdy = (te(k,j,i)-te(k,j-1,i))/pco_->h2v(i)/pco_->dx2v(j-1);
          
          // Monotonized temperature difference dT/dx
          dTdx = limiter4(te(k,j  ,i+1) - te(k,j  ,i  ),
                          te(k,j  ,i  ) - te(k,j  ,i-1),
                          te(k,j-1,i+1) - te(k,j-1,i  ),
                          te(k,j-1,i  ) - te(k,j-1,i-1));
          dTdx /= pco_->dx1v(i-1);
          
          // Monotonized temperature difference dT/dz, 3D problem ONLY
          if (pmb_->block_size.nx3 > 1) {
            dTdz = limiter4(te(k+1,j  ,i) - te(k  ,j  ,i),
                            te(k  ,j  ,i) - te(k-1,j  ,i),
                            te(k+1,j-1,i) - te(k  ,j-1,i),
                            te(k  ,j-1,i) - te(k-1,j-1,i));
            dTdz /= (pco_->dx3v(k-1)*pco_->h31v(i)*pco_->h32v(j));
          }
          
          // Add flux at x2-interface, 2D PROBLEM
          if (pmb_->block_size.nx3 == 1) {
            Bx = 0.5*(b(0,k,j-1,i) + b(0,k,j,i));
            By = b(1,k,j,i);
            B02 = SQR(Bx) + SQR(By);
            B02 = max(B02,TINY_NUMBER); // limit in case B=0
            bDotGradT = Bx*dTdx + By*dTdy;
            kd = kappaf*denf;
            x2flux(k,j,i) -= kd*(By*bDotGradT)/B02;

          // Add flux at x2-interface, 3D PROBLEM
          } else {
            Bx = 0.5*(b(0,k,j-1,i) + b(0,k,j,i));
            By = b(1,k,j,i);
            Bz = 0.5*(b(2,k,j-1,i) + b(2,k,j,i));
            B02 = SQR(Bx) + SQR(By) + SQR(Bz);
            B02 = max(B02,TINY_NUMBER); // limit in case B=0
            bDotGradT = Bx*dTdx + By*dTdy + Bz*dTdz;
            kd = kappaf*denf;
            x2flux(k,j,i) -= kd*(By*bDotGradT)/B02;
          }
        }
      }
    }
  }
  
  // k-direction
  il=is, iu=ie, jl=js, ju=je;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) // 2D or 3D
      il=is-1, iu=ie+1, jl=js-1, ju=je+1;
    else // 1D
      il=is-1, iu=ie+1;
  }
  if(pmb_->block_size.nx3 > 1) { //3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
  #pragma omp simd
        for (int i=il; i<=iu; ++i) {
          kappaf = 0.5*(kappa(ANI,k,j,i)+kappa(ANI,k-1,j,i));
          denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k-1,j,i));
          dTdz = (te(k,j,i)-te(k-1,j,i))/pco_->dx3v(k-1)/pco_->h31v(i)/pco_->h32v(j);
          
          // Monotonized temperature difference dT/dx
          dTdx = limiter4(te(k  ,j,i+1) - te(k  ,j,i  ),
                          te(k  ,j,i  ) - te(k  ,j,i-1),
                          te(k-1,j,i+1) - te(k-1,j,i  ),
                          te(k-1,j,i  ) - te(k-1,j,i-1));
          dTdx /= pco_->dx1v(i-1);
          
          // Monotonized temperature difference dT/dy
          dTdy = limiter4(te(k  ,j+1,i) - te(k  ,j  ,i),
                          te(k  ,j  ,i) - te(k  ,j-1,i),
                          te(k-1,j+1,i) - te(k-1,j  ,i),
                          te(k-1,j  ,i) - te(k-1,j-1,i));
          dTdy /= (pco_->h2v(i)*pco_->dx2v(j-1));
          
          // Add flux at x3-interface, 3D PROBLEM
          Bx = 0.5*(b(0,k-1,j,i) + b(0,k,j,i));
          By = 0.5*(b(1,k-1,j,i) + b(1,k,j,i));
          Bz = b(2,k,j,i);
          B02 = SQR(Bx) + SQR(By) + SQR(Bz);
          B02 = max(B02,TINY_NUMBER); // limit in case B=0
          bDotGradT = Bx*dTdx + By*dTdy + Bz*dTdz;
          kd = kappaf*denf;
          x3flux(k,j,i) -= kd*(Bz*bDotGradT)/B02;
        }
      }
    }
  }
  
  // Release memory
  te.DeleteAthenaArray();

  return;
}




//------------------------------------------------------------------------------
// constant conduction

void ConstConduction(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bcc, int is, int ie, int js, int je, int ks, int ke) {
  if (phdif->kappa_iso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->kappa(ISO,k,j,i) = phdif->kappa_iso;
      }
    }
  }
  if (phdif->kappa_aniso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->kappa(ANI,k,j,i) = phdif->kappa_aniso;
      }
    }
  }
  return;
}
