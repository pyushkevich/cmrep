#include "MedialAtom.h"

MedialAtom
::MedialAtom()
{
  this->SetAllDerivativeTermsToZero();
}

void
MedialAtom
::SetAllDerivativeTermsToZero()
{
  X.fill(0.0); Xu.fill(0.0); Xv.fill(0.0); Xuu.fill(0.0); Xuv.fill(0.0); Xvv.fill(0.0);
  F = 0.0; Fu = 0.0; Fv = 0.0; R = 0.0; Ru = 0.0; Rv = 0.0;
  Fuu = 0.0; Fuv = 0.0; Fvv = 0.0;

  G.g = G.gInv = 0.0;
  for(size_t a1 = 0; a1 < 2; a1++)
    for(size_t a2 = 0; a2 < 2; a2++)
      {
      G.xCovariantTensor[a1][a2] = 0.0;
      G.xContravariantTensor[a1][a2] = 0.0;
      for(size_t a3 = 0; a3 < 2; a3++)
        {
        G.xChristoffelFirst[a1][a2][a3] = 0.0;
        G.xChristoffelSecond[a1][a2][a3] = 0.0;
        }
      }

  aelt = 0.0;

  N.fill(0.0);
  xGradR.fill(0.0);
  xLapR = 0.0;
  xGradRMagSqr = 0.0; xNormalFactor = 0.0; xGradRMagSqrOrig = 0.0;
  Rs2 = 0;

  xBnd[0].X.fill(0.0); xBnd[0].N.fill(0.0);
  xBnd[1].X.fill(0.0); xBnd[1].N.fill(0.0);

  xBnd[0].curv_mean = 0.0; xBnd[0].curv_gauss = 0.0;
  xBnd[1].curv_mean = 0.0; xBnd[1].curv_gauss = 0.0;


  xMeanCurv = xGaussCurv = 0.0;

  // By default, atoms are assumed to depend on variation
  order = 0;

  // Reset the crest flag
  flagCrest = false;
  flagValid = false;
}

void
MedialAtom
::ComputeDifferentialGeometry()
{
  G.SetTwoJet( X.data_block(), Xu.data_block(), Xv.data_block(),
    Xuu.data_block(), Xuv.data_block(), Xvv.data_block());
}

void MedialAtom ::ComputeNormalVector()
{
  aelt = sqrt(G.g);
  N = vnl_cross_3d(Xu, Xv) / aelt;
}

bool MedialAtom::ComputeBoundaryAtoms(bool flagEdgeAtom)
{
  // Compute the Riemannian gradient of Phi
  double (*G2)[2] = G.xContravariantTensor;
  double Cu = (Fu * G2[0][0] + Fv * G2[0][1]);
  double Cv = (Fu * G2[0][1] + Fv * G2[1][1]);
  SMLVec3d xGradPhi = Xu * Cu + Xv * Cv;

  // Compute the radius of the atom
  R = sqrt(F);

  // Compute the Riemannian gradient of R
  xGradR = 0.5 * xGradPhi / R;

  // Compute the gradient magnitude of R
  xGradRMagSqr = 0.25 * (Fu * Cu + Fv * Cv) / F;

  // Split depending on whether this is an end atom
  if(flagEdgeAtom)
    {
    // We are at a crest, valid atom
    flagCrest = true;
    flagValid = true;

    // There is a zero badness value
    xGradRMagSqr = 1.0;
    xNormalFactor = 0.0;

    // Compute (dR / ds)^2
    Rs2 = Fv * Fv / (4 * F * G.xCovariantTensor[1][1]);

    // Basically, we expect each medial model to give us boundary atoms
    // with gradR ~= 1. We simply set xBnd.N = -gradR. If the magnitude
    // of xBnd.N is not exactly 1, that can be a problem - but it's the
    // problem that should be handled in the medial model, not here.
    xBnd[0].N = xBnd[1].N = - xGradR;
    xBnd[0].X = xBnd[1].X = X + R * xBnd[0].N;
    }
  else
    {
    // Otherwise, we have a normal atom
    flagCrest = false;

    // Compute the Riemannian gradient of R
    xGradR = 0.5 * xGradPhi / R;

    // Compute the gradient magnitude of R
    xGradRMagSqr = 0.25 * (Fu * Cu + Fv * Cv) / F;
    // xGradRMagSqr = xGradR.squared_magnitude();

    // Check if this is greater than one - should never be very close either
    if(xGradRMagSqr > 1.0)
      {
      // This is an invalid atom
      flagValid = false;

      // Set the boundary atoms to zero
      xNormalFactor = 0;
      xBnd[0].X = xBnd[1].X = X;
      xBnd[0].N = -N;
      xBnd[1].N = N;
      }
    else
      {
      // Finally, we have a clean internal atom
      flagValid = true;

      // Compute the normal component of the sail vectors
      xNormalFactor = sqrt(1.0 - xGradRMagSqr);
      SMLVec3d CN = N * xNormalFactor;

      // Compute the position and normals of boundary points
      xBnd[0].N = - xGradR - CN;
      xBnd[0].X = X + xBnd[0].N * R;
      xBnd[1].N = - xGradR + CN;
      xBnd[1].X = X + xBnd[1].N * R;
      }
    }

  return flagValid;
}

bool MedialAtom::ComputeBoundaryAtomsUsingR(bool flagEdgeAtom)
{
  // Compute the Riemannian gradient of R
  double (*G2)[2] = G.xContravariantTensor;
  double Cu = (Ru * G2[0][0] + Rv * G2[0][1]);
  double Cv = (Ru * G2[0][1] + Rv * G2[1][1]);

  // Compute gradR and its gradient magnitude
  xGradR       = Xu * Cu + Xv * Cv;
  xGradRMagSqr = Ru * Cu + Rv * Cv;

  // Compute F = r^2
  F = R * R;

  // Split depending on whether this is an end atom
  if(flagEdgeAtom)
    {
    // We are at a crest, valid atom
    flagCrest = true;

    // There is a zero badness value
    xNormalFactor = 0.0;

    // Compute (dR/ds) ^ 2
    Rs2 = Rv * Rv / G.xCovariantTensor[1][1];

    // If the magnitude of gradR is zero, we are in real trouble
    if(xGradRMagSqr == 0.0)
      {
      flagValid = false;
      xBnd[0].X = xBnd[1].X = X;
      xBnd[0].N = -N;
      xBnd[1].N = N;
      }
    else
      {
      flagValid = true;

      // Use a different computation for the boundary atoms. This
      // should be a correct computation if |gradR|==1 but it will
      // different for brute force models where gradR may not be 1
      // on the boundary
      xBnd[0].N = xBnd[1].N = - xGradR / sqrt(xGradRMagSqr);
      xBnd[0].X = xBnd[1].X = X + R * xBnd[0].N;
      }
    }
  else
    {
    // Otherwise, we have a normal atom
    flagCrest = false;

    // Check if this is greater than one - should never be very close either
    if(xGradRMagSqr > 1.0)
      {
      // This is an invalid atom
      flagValid = false;

      // Set the boundary atoms to zero
      xNormalFactor = 0;
      xBnd[0].X = xBnd[1].X = X;
      xBnd[0].N = -N;
      xBnd[1].N = N;
      }
    else
      {
      // Finally, we have a clean internal atom
      flagValid = true;

      // Compute the normal component of the sail vectors
      xNormalFactor = sqrt(1.0 - xGradRMagSqr);
      SMLVec3d CN = N * xNormalFactor;

      // Compute the position and normals of boundary points
      xBnd[0].N = - xGradR - CN;
      xBnd[0].X = X + xBnd[0].N * R;
      xBnd[1].N = - xGradR + CN;
      xBnd[1].X = X + xBnd[1].N * R;
      }
    }

  return flagValid;
}



// Prerequisites:
//  * The first jet for the atom and datom
//  * The contravariant tensor for the atom and datom
//
void 
MedialAtom::ComputeBoundaryAtomDerivatives(
  MedialAtom &dAtom, const DerivativeTerms &dt) const
{
  // Get the relevant elements of the atoms
  SMLVec3d &Y = dAtom.X, &Yu = dAtom.Xu, &Yv = dAtom.Xv;
  
  // Get the elements of the first fundamental form and its derivative
  const double &g11 = G.xContravariantTensor[0][0];
  const double &g12 = G.xContravariantTensor[0][1];
  const double &g22 = G.xContravariantTensor[1][1];
  const double &z11 = dAtom.G.xContravariantTensor[0][0];
  const double &z12 = dAtom.G.xContravariantTensor[0][1];
  const double &z22 = dAtom.G.xContravariantTensor[1][1];

  // Get the g's
  const double &g = G.g; const double &z = dAtom.G.g;

  // Get the partials of Phi and its variational derivative
  double &H = dAtom.F, &Hu = dAtom.Fu, &Hv = dAtom.Fv;

  // Get the derivatives of R
  double P = H * dt.x1_2R;
  double Pu = Hu * dt.x1_2R - H * dt.Ru_2F;
  double Pv = Hv * dt.x1_2R - H * dt.Rv_2F;

  dAtom.R = P; 
  
  // This is the derivative of the normal vector
  dAtom.N =  vnl_cross_3d(dt.Xu_aelt, Yv);
  dAtom.N += vnl_cross_3d(Yu, dt.Xv_aelt);
  vmuladd(dAtom.N, dt.N_2g, -z);

  // We will compute several intermediate terms
  double z1iRi = z11 * dt.Ru + z12 * dt.Rv;
  double z2iRi = z12 * dt.Ru + z22 * dt.Rv;

  // Compute the derivative of grad R
  vmulset(dAtom.xGradR, Xu, z1iRi + g11 * Pu + g12 * Pv);
  vmuladd(dAtom.xGradR, Xv, z2iRi + g12 * Pu + g22 * Pv);
  vmuladd(dAtom.xGradR, Yu, dt.g1iRi);
  vmuladd(dAtom.xGradR, Yv, dt.g2iRi);

  // Compute the derivative of grad R . grad R
  dAtom.xGradRMagSqr = z1iRi * dt.Ru + z2iRi * dt.Rv 
    + 2.0 * (dt.g1iRi * Pu + dt.g2iRi * Pv);

  /*
  double test = 
    z11 * dt.Ru * dt.Ru + 2 * g11 * Pu * dt.Ru +
    z22 * dt.Rv * dt.Rv + 2 * g22 * Pv * dt.Rv +
    2.0 * (z12 * dt.Ru * dt.Rv + g12 * Pu * dt.Rv + g12 * dt.Ru * Pv);
    */

  // Compute the 

  // Address the edge case first
  if(flagCrest) 
    {
    dAtom.xNormalFactor = 0.0;

    // Compute the derivative of (dR/ds)^2
    dAtom.Rs2 = (Fv == 0) ? 0 : Rs2 * (
      - H / F 
      - dAtom.G.xCovariantTensor[1][1] / G.xCovariantTensor[1][1]
      + 2 * Hv / Fv); 

    if(flagValid)
      {
      // Compute (gradR / |gradR|)'
      vmulset(dAtom.xBnd[0].N, dAtom.xGradR, -dt.inv_mag_gradR);
      vmuladd(dAtom.xBnd[0].N, xBnd[0].N, 
        - dot_product(dAtom.xBnd[0].N, xBnd[0].N));
      dAtom.xBnd[1].N = dAtom.xBnd[0].N;

      // Compute the boundary atoms
      dAtom.xBnd[1].X = dAtom.xBnd[0].X = 
        Y + R * dAtom.xBnd[0].N + P * xBnd[0].N;


      // The normal term vanishes
      // dAtom.xBnd[0].N = dAtom.xBnd[1].N = - dAtom.xGradR;
      // dAtom.xGradRMagSqr = 0.0;
      }
    else
      {
      dAtom.xBnd[0].X = dAtom.xBnd[1].X = dAtom.X;
      dAtom.xBnd[0].N = -dAtom.N;
      dAtom.xBnd[1].N = dAtom.N;
      }
    }
  else
    {
    if(flagValid)
      {
      // Set the normal factor
      dAtom.xNormalFactor = -0.5 * dAtom.xGradRMagSqr / xNormalFactor;

      // Compute the plus-minus term
      SMLVec3d dNormalTerm;
      vmulset(dNormalTerm, dAtom.N, xNormalFactor);
      vmuladd(dNormalTerm, dt.N_2nt, -dAtom.xGradRMagSqr);

      // Compute the boundary atom normals
      dAtom.xBnd[0].N = dAtom.xBnd[1].N = - dAtom.xGradR;
      dAtom.xBnd[0].N -= dNormalTerm;
      dAtom.xBnd[1].N += dNormalTerm;

      // Compute the boundary atoms
      dAtom.xBnd[0].X = dAtom.xBnd[1].X = Y;
      vmuladd(dAtom.xBnd[0].X, dAtom.xBnd[0].N, R);
      vmuladd(dAtom.xBnd[0].X, xBnd[0].N, P);
      vmuladd(dAtom.xBnd[1].X, dAtom.xBnd[1].N, R);
      vmuladd(dAtom.xBnd[1].X, xBnd[1].N, P);
      }
    else
      {
      dAtom.xBnd[0].X = dAtom.xBnd[1].X = dAtom.X;
      dAtom.xBnd[0].N = -dAtom.N;
      dAtom.xBnd[1].N = dAtom.N;
      }
    }
}


// Prerequisites:
//  * The first jet for the atom and datom
//  * The contravariant tensor for the atom and datom
//
void 
MedialAtom::ComputeBoundaryAtomDerivativesUsingR(
  MedialAtom &dAtom, const DerivativeTerms &dt) const
{
  // Get the relevant elements of the atoms
  SMLVec3d &Y = dAtom.X, &Yu = dAtom.Xu, &Yv = dAtom.Xv;
  
  // Get the elements of the first fundamental form and its derivative
  const double &g11 = G.xContravariantTensor[0][0];
  const double &g12 = G.xContravariantTensor[0][1];
  const double &g22 = G.xContravariantTensor[1][1];
  const double &z11 = dAtom.G.xContravariantTensor[0][0];
  const double &z12 = dAtom.G.xContravariantTensor[0][1];
  const double &z22 = dAtom.G.xContravariantTensor[1][1];

  // Get the g's
  const double &z = dAtom.G.g;

  // Get the partials of Phi and its variational derivative
  double P = dAtom.R;
  double Pu = dAtom.Ru;
  double Pv = dAtom.Rv;

  double Cu = (Ru * g11 + Rv * g12);
  double Cv = (Ru * g12 + Rv * g22);
  double dCu = (Pu * g11 + Ru * z11 + Pv * g12 + Rv * z12);
  double dCv = (Pu * g12 + Ru * z12 + Pv * g22 + Rv * z22);

  // Compute gradR and its gradient magnitude
  dAtom.xGradR       = Xu * dCu + Yu * Cu + Xv * dCv + Yv * Cv;
  dAtom.xGradRMagSqr = Ru * dCu + Pu * Cu + Rv * dCv + Pv * Cv;

  // Set the partials of F
  dAtom.F = 2.0 * R * P; 
  
  // This is the derivative of the normal vector
  dAtom.N =  vnl_cross_3d(dt.Xu_aelt, Yv);
  dAtom.N += vnl_cross_3d(Yu, dt.Xv_aelt);
  vmuladd(dAtom.N, dt.N_2g, -z);

  // Address the edge case first
  if(flagCrest) 
    {
    dAtom.xNormalFactor = 0.0;

    dAtom.Rs2 = (Rv == 0) ? 0 : Rs2 * (
      - dAtom.G.xCovariantTensor[1][1] / G.xCovariantTensor[1][1] 
      + 2 * Pv / Rv);

    if(flagValid)
      {
      // Compute (gradR / |gradR|)'
      vmulset(dAtom.xBnd[0].N, dAtom.xGradR, -dt.inv_mag_gradR);
      vmuladd(dAtom.xBnd[0].N, xBnd[0].N, 
        - dot_product(dAtom.xBnd[0].N, xBnd[0].N));
      dAtom.xBnd[1].N = dAtom.xBnd[0].N;

      // Compute the boundary atoms
      dAtom.xBnd[1].X = dAtom.xBnd[0].X = 
        Y + R * dAtom.xBnd[0].N + P * xBnd[0].N;
      }
    else
      {
      dAtom.xBnd[0].X = dAtom.xBnd[1].X = dAtom.X;
      dAtom.xBnd[0].N = -dAtom.N;
      dAtom.xBnd[1].N = dAtom.N;
      }
    }
  else
    {
    if(flagValid)
      {
      // Set the normal factor
      dAtom.xNormalFactor = -0.5 * dAtom.xGradRMagSqr / xNormalFactor;

      // Compute the plus-minus term
      SMLVec3d dNormalTerm;
      vmulset(dNormalTerm, dAtom.N, xNormalFactor);
      vmuladd(dNormalTerm, dt.N_2nt, -dAtom.xGradRMagSqr);

      // Compute the boundary atom normals
      dAtom.xBnd[0].N = dAtom.xBnd[1].N = - dAtom.xGradR;
      dAtom.xBnd[0].N -= dNormalTerm;
      dAtom.xBnd[1].N += dNormalTerm;

      // Compute the boundary atoms
      dAtom.xBnd[0].X = dAtom.xBnd[1].X = Y;
      vmuladd(dAtom.xBnd[0].X, dAtom.xBnd[0].N, R);
      vmuladd(dAtom.xBnd[0].X, xBnd[0].N, P);
      vmuladd(dAtom.xBnd[1].X, dAtom.xBnd[1].N, R);
      vmuladd(dAtom.xBnd[1].X, xBnd[1].N, P);
      }
    else
      {
      dAtom.xBnd[0].X = dAtom.xBnd[1].X = dAtom.X;
      dAtom.xBnd[0].N = -dAtom.N;
      dAtom.xBnd[1].N = dAtom.N;
      }
    }
}


double temp01(
  SMLVec3d &Au, SMLVec3d &Bu, SMLVec3d &Cu,
  SMLVec3d &Av, SMLVec3d &Bv, SMLVec3d &Cv)
{

  double rtn = 0.0;
  rtn += (dot_product(Au, Au) + 2*dot_product(Bu, Cu)) * dot_product(Cu, Cu);
  rtn += (dot_product(Av, Av) + 2*dot_product(Bv, Cv)) * dot_product(Cv, Cv);
  rtn -= (dot_product(Au, Av) + dot_product(Bu, Cv) + dot_product(Bv, Cu))
    * 0.5 * dot_product(Cu, Cv);
  rtn += 4 * dot_product(Au, Cu) * dot_product(Av, Cv);

  double t1 = dot_product(Au, Cv) + dot_product(Av, Cu);
  rtn -= t1 * t1;
  return rtn;
}

void MedialAtom::ComputeBoundaryCurvature()
{
  // Just set up an edge walk
  if(flagValid && !flagCrest)
    {
    // Create arrays to allow easier looping
    SMLVec3d Xi[2] = {Xu, Xv};
    SMLVec3d Xij[2][2] = {{Xuu, Xuv},{Xuv, Xvv}};
    double Fi[2] = {Fu, Fv};
    double Fij[2][2] = {{Fuu, Fuv},{Fuv, Fvv}};

    // Compute the partial derivatives of GradF
    size_t i,j,k,l;
    for(i=0;i<2;i++) 
      {
      // Derivative of gradF wrt U^i
      gradF_i[i].fill(0.0);

      // Sum over j, k 
      for(j=0;j<2;j++) for(k=0;k<2;k++)
        {
        // Compute dg^jk/du^i
        double dg_jk_d_ui = 0.0;
        for(l=0;l<2;l++)
          {
          dg_jk_d_ui += -(
            G.xContravariantTensor[l][k] * G.xChristoffelSecond[l][i][j] +
            G.xContravariantTensor[j][l] * G.xChristoffelSecond[l][i][k]);
          }

        // Compute contribution to gradFi
        gradF_i[i] += 
          (dg_jk_d_ui * Fi[k]) * Xi[j]
          + (G.xContravariantTensor[j][k] *  Fi[k]) * Xij[j][i] 
          + (G.xContravariantTensor[j][k] * Fij[k][i]) * Xi[j];
        }

      // Compute the partial derivative of gradR
      gradR_i[i] = (gradF_i[i] * R - xGradR * Fi[i]) / (2.0 * F);

      // Compute the partial derivative of the unit normal to the medial axis
      double dg_di = 2 * G.g * (
        G.xChristoffelSecond[0][i][0] + G.xChristoffelSecond[1][i][1]);
      SMLVec3d dN0_di = 
        vnl_cross_3d(Xij[0][i], Xi[1]) + vnl_cross_3d(Xi[0], Xij[1][i]);
      N_i[i] = (dN0_di * aelt - 0.5 * N * dg_di) / (G.g);
      }

    // Repeat for both sides of the atom
    for(size_t side = 0; side < 2; side++)
      {
      // Get the boundary atom
      BoundaryAtom &ba = xBnd[side];

      for(i=0;i<2;i++) 
        {
        // Compute the boundary node's normal
        ba.N_i[i] = 
          - gradR_i[i] + (side ? 1.0 : -1.0) * (
            - dot_product(xGradR, gradR_i[i]) * N / xNormalFactor +
            xNormalFactor * N_i[i]);

        // Compute the boundary node's position
        ba.X_i[i] = Xi[i] + (Fi[i] / (2 * R)) * ba.N + R * ba.N_i[i];
        }

      // Now that we have the first derivatives of the normal and X on the boundary, curvature
      // computation is easy
      SMLVec3d &Xu = ba.X_i[0], &Xv = ba.X_i[1], &Nu = ba.N_i[0], &Nv = ba.N_i[1];
      ba.ee = - dot_product(Nv,Xv);
      ba.ff = -0.5 * (dot_product(Nu, Xv) + dot_product(Nv, Xu));
      ba.gg = - dot_product(Nu, Xu);
      ba.E = dot_product(Xu, Xu);
      ba.F = dot_product(Xu, Xv);
      ba.G = dot_product(Xv, Xv);
      ba.det_g = (ba.E*ba.G - ba.F*ba.F);

      ba.curv_mean= (ba.ee * ba.G - 2 * ba.ff * ba.F + ba.gg * ba.E) / (2 * ba.det_g);
      ba.curv_gauss = (ba.ee * ba.gg - ba.ff * ba.ff) / ba.det_g;
      }
    }
  else
    {
    for(size_t side = 0; side < 2; side++)
      {
      BoundaryAtom &ba = xBnd[side];
      ba.curv_gauss = ba.curv_mean = 0.0;
      }
    }
}

void MedialAtom::ComputeBoundaryCurvatureDerivative(MedialAtom &da)
{
  if(flagValid && !flagCrest && da.order < 3)
    {
    // Create arrays to allow easier looping
    SMLVec3d Xi[2] = {Xu, Xv};
    SMLVec3d Xij[2][2] = {{Xuu, Xuv},{Xuv, Xvv}};
    SMLVec3d DXi[2] = {da.Xu, da.Xv};
    SMLVec3d DXij[2][2] = {{da.Xuu, da.Xuv},{da.Xuv, da.Xvv}};

    double Fi[2] = {Fu, Fv};
    double Fij[2][2] = {{Fuu, Fuv},{Fuv, Fvv}};
    double DFi[2] = {da.Fu, da.Fv};
    double DFij[2][2] = {{da.Fuu, da.Fuv},{da.Fuv, da.Fvv}};

    // Compute the partial derivatives of GradF
    size_t i,j,k,l;
    for(i=0;i<2;i++) 
      {
      // Derivative of gradF wrt U^i
      da.gradF_i[i].fill(0.0);

      // Sum over j, k 
      for(j=0;j<2;j++) for(k=0;k<2;k++)
        {
        // Compute dg^jk/du^i
        double dg_jk_d_ui = 0.0;
        double d_dg_jk_d_ui = 0.0;
        for(l=0;l<2;l++)
          {
          dg_jk_d_ui += -(
            G.xContravariantTensor[l][k] * G.xChristoffelSecond[l][i][j] +
            G.xContravariantTensor[j][l] * G.xChristoffelSecond[l][i][k]);

          d_dg_jk_d_ui += -(
            da.G.xContravariantTensor[l][k] * G.xChristoffelSecond[l][i][j] +
            G.xContravariantTensor[l][k] * da.G.xChristoffelSecond[l][i][j] +
            da.G.xContravariantTensor[j][l] * G.xChristoffelSecond[l][i][k] +
            G.xContravariantTensor[j][l] * da.G.xChristoffelSecond[l][i][k]);
          }

        // Holy cow, what a mess!
        da.gradF_i[i] += 
          (d_dg_jk_d_ui * Fi[k]) * Xi[j]
          + (dg_jk_d_ui * DFi[k]) * Xi[j]
          + (dg_jk_d_ui * Fi[k]) * DXi[j]
          + (da.G.xContravariantTensor[j][k] *  Fi[k]) * Xij[j][i] 
          + (G.xContravariantTensor[j][k] *  DFi[k]) * Xij[j][i] 
          + (G.xContravariantTensor[j][k] *  Fi[k]) * DXij[j][i] 
          + (da.G.xContravariantTensor[j][k] * Fij[k][i]) * Xi[j]
          + (G.xContravariantTensor[j][k] * DFij[k][i]) * Xi[j]
          + (G.xContravariantTensor[j][k] * Fij[k][i]) * DXi[j];
        }

      // cout << "d_gradF_i[" << i << " = " << d_gradF_i[i] << endl;

      // Compute the partial derivative of gradR
      da.gradR_i[i] = 
        ((da.gradF_i[i] * R + gradF_i[i] * da.R 
          - (da.xGradR * Fi[i] + xGradR * DFi[i])) * (F)
         - (gradF_i[i] * R - xGradR * Fi[i]) * (da.F))
        / (2.0 * F * F);

      // cout << "d_gradR_i[" << i << " = " << d_gradR_i[i] << endl;

      // Compute the partial derivative of the unit normal to the medial axis
      double dg_di = 2 * G.g * (
        G.xChristoffelSecond[0][i][0] + G.xChristoffelSecond[1][i][1]);
      double d_dg_di = 
        2 * da.G.g * (
          G.xChristoffelSecond[0][i][0] + G.xChristoffelSecond[1][i][1])
        + 2 * G.g * (
          da.G.xChristoffelSecond[0][i][0] + da.G.xChristoffelSecond[1][i][1]);

      SMLVec3d dN0_di = 
        vnl_cross_3d(Xij[0][i], Xi[1]) + vnl_cross_3d(Xi[0], Xij[1][i]);
      SMLVec3d d_dN0_di = 
        vnl_cross_3d(DXij[0][i], Xi[1]) 
        + vnl_cross_3d(Xij[0][i], DXi[1]) 
        + vnl_cross_3d(DXi[0], Xij[1][i])
        + vnl_cross_3d(Xi[0], DXij[1][i]);

      da.N_i[i] = (
        (d_dN0_di * aelt + dN0_di * da.aelt - 0.5 * 
         (da.N * dg_di + N * d_dg_di)) * G.g
        - (dN0_di * aelt - 0.5 * N * dg_di) * da.G.g) / (G.g * G.g);
    
      // Compute some intermediate terms first
      double tmp1 = - dot_product(xGradR, gradR_i[i]);
      double d_tmp1 = - dot_product(da.xGradR, gradR_i[i]) 
        - dot_product(xGradR, da.gradR_i[i]);

      SMLVec3d tmp2 = N / xNormalFactor;
      SMLVec3d d_tmp2 = (da.N * xNormalFactor - N * da.xNormalFactor) /
        (xNormalFactor * xNormalFactor);

      // Compute the boundary node's position
      double Ri = 0.5 * Fi[i] / R;
      double d_Ri = 0.5 * (DFi[i] * R - Fi[i] * da.R) / F;

      // Repeat for each side
      for(size_t side = 0; side < 2; side++)
        {
        BoundaryAtom &ba = xBnd[side];
        BoundaryAtom &dba = da.xBnd[side];

        // Compute the boundary node's normal
        dba.N_i[i] = 
          - da.gradR_i[i] + (side ? 1.0 : -1.0) * (
            d_tmp1 * tmp2 + tmp1 * d_tmp2 
            + da.xNormalFactor * N_i[i] + xNormalFactor * da.N_i[i]);

        dba.X_i[i] = DXi[i] + d_Ri * ba.N + Ri * dba.N + da.R * ba.N_i[i] + R * dba.N_i[i];
        }
      }

    // Compute the curvatures
    for(size_t side = 0; side < 2; side++)
      {
      BoundaryAtom &ba = xBnd[side];
      BoundaryAtom &dba = da.xBnd[side];

      // Now that we have the first derivatives of the normal and X on the boundary, curvature
      // computation is easy
      SMLVec3d &Xu = ba.X_i[0], &Xv = ba.X_i[1], &Nu = ba.N_i[0], &Nv = ba.N_i[1];
      SMLVec3d &dXu = dba.X_i[0], &dXv = dba.X_i[1], &dNu = dba.N_i[0], &dNv = dba.N_i[1];

      dba.ee = - dot_product(dNv,Xv) - dot_product(Nv, dXv);
      dba.ff = - 0.5 * (
        dot_product(dNu, Xv) + dot_product(Nu, dXv) + 
        dot_product(dNv, Xu) + dot_product(Nv, dXu));
      dba.gg = - dot_product(dNu, Xu) - dot_product(Nu, dXu);

      dba.E = 2.0 * dot_product(Xu, dXu);
      dba.F = dot_product(dXu, Xv) + dot_product(Xu, dXv);
      dba.G = 2.0 * dot_product(Xv, dXv);

      dba.det_g = (dba.E*ba.G + ba.E*dba.G - 2*ba.F*dba.F);

      // cout << "Mean Curv " << (ee * G - 2 * ff * F + gg * E) / (2 * det);
      // cout << "Gauss Crv " << (ee * gg - ff * ff) / det;
      dba.curv_mean = 0.5 *
        ((dba.ee * ba.G + ba.ee * dba.G - 
          2 * (dba.ff * ba.F + ba.ff * dba.F) + 
          dba.gg * ba.E + ba.gg * dba.E) * ba.det_g 
         - (ba.ee * ba.G - 2 * ba.ff * ba.F + ba.gg * ba.E) * dba.det_g) / (ba.det_g * ba.det_g);

      dba.curv_gauss = (
        (dba.ee * ba.gg + ba.ee * dba.gg - 2 * ba.ff * dba.ff) * ba.det_g
        - (ba.ee * ba.gg - ba.ff * ba.ff) * dba.det_g) / (ba.det_g * ba.det_g);
      }
    }
}


void
MedialAtom
::ComputeCommonDerivativeTerms(MedialAtom::DerivativeTerms &dt) const
{
  // Get the elements of the first fundamental form and its derivative
  const double &g11 = G.xContravariantTensor[0][0];
  const double &g12 = G.xContravariantTensor[0][1];
  const double &g22 = G.xContravariantTensor[1][1];

  // Get the derivatives of R
  dt.x1_2R = 0.5 / R;
  dt.Ru = Fu * dt.x1_2R, dt.Rv = Fv * dt.x1_2R;
  dt.Ru_R = dt.Ru / R;
  dt.Rv_R = dt.Rv / R;
  dt.Ru_2F = 0.5 * dt.Ru / F;
  dt.Rv_2F = 0.5 * dt.Rv / F;
  
  // Terms used to compute the derivative of the normal vector
  dt.Xu_aelt = Xu / aelt; dt.Xv_aelt = Xv / aelt;
  dt.N_2g = 0.5 * N / G.g;

  // We will compute several intermediate terms
  dt.g1iRi = g11 * dt.Ru + g12 * dt.Rv;
  dt.g2iRi = g12 * dt.Ru + g22 * dt.Rv;

  // Compute the plus-minus term
  dt.N_2nt = 0.5 * N / xNormalFactor;

  // Compute the inverse magnitude of grad R
  dt.inv_mag_gradR = (flagValid) ?  1.0 / xGradR.magnitude() : 0.0;
}

// Compute directional derivative of the contravariant tensor given the
// directional derivatives of X, Xu and Xv contained in dAtom
void
MedialAtom
::ComputeMetricTensorDerivatives(MedialAtom &dAtom) const
{
  // Compute the derivatives of the covariant tensor
  dAtom.G.xCovariantTensor[0][0] = 2 * dot(dAtom.Xu, Xu);
  dAtom.G.xCovariantTensor[1][1] = 2 * dot(dAtom.Xv, Xv);
  dAtom.G.xCovariantTensor[0][1] = 
    dAtom.G.xCovariantTensor[1][0] = dot(dAtom.Xv, Xu) + dot(dAtom.Xu, Xv);

  // Compute the derivative of g
  dAtom.G.g = 
    dAtom.G.xCovariantTensor[0][0] * G.xCovariantTensor[1][1] +
    dAtom.G.xCovariantTensor[1][1] * G.xCovariantTensor[0][0] -
    2 * dAtom.G.xCovariantTensor[0][1] * G.xCovariantTensor[0][1];
  dAtom.G.gInv = - G.gInv * (G.gInv * dAtom.G.g);

  // Compute the derivatives of the contravariant tensor
  dAtom.G.xContravariantTensor[0][0] = G.gInv *
    (dAtom.G.xCovariantTensor[1][1] - G.xContravariantTensor[0][0] * dAtom.G.g);

  dAtom.G.xContravariantTensor[1][1] = G.gInv *
    (dAtom.G.xCovariantTensor[0][0] - G.xContravariantTensor[1][1] * dAtom.G.g);

  dAtom.G.xContravariantTensor[0][1] = dAtom.G.xContravariantTensor[1][0] = - G.gInv *
    (dAtom.G.xCovariantTensor[0][1] + G.xContravariantTensor[0][1] * dAtom.G.g);

  // Compute the area element
  dAtom.aelt = 0.5 * dAtom.G.g / aelt;
}

void
MedialAtom
::ComputeChristoffelDerivatives(MedialAtom &dAtom) const
{
  // Compute the derivative of the first kind
  dAtom.G.xChristoffelFirst[0][0][0] = 
    dot_product(dAtom.Xuu, Xu) + dot_product(Xuu, dAtom.Xu);
  dAtom.G.xChristoffelFirst[0][0][1] = 
    dot_product(dAtom.Xuu, Xv) + dot_product(Xuu, dAtom.Xv);
  dAtom.G.xChristoffelFirst[0][1][0] = 
    dot_product(dAtom.Xuv, Xu) + dot_product(Xuv, dAtom.Xu);
  dAtom.G.xChristoffelFirst[0][1][1] = 
    dot_product(dAtom.Xuv, Xv) + dot_product(Xuv, dAtom.Xv);
  dAtom.G.xChristoffelFirst[1][1][0] = 
    dot_product(dAtom.Xvv, Xu) + dot_product(Xvv, dAtom.Xu);
  dAtom.G.xChristoffelFirst[1][1][1] = 
    dot_product(dAtom.Xvv, Xv) + dot_product(Xvv, dAtom.Xv);
  dAtom.G.xChristoffelFirst[1][0][0] = dAtom.G.xChristoffelFirst[0][1][0];
  dAtom.G.xChristoffelFirst[1][0][1] = dAtom.G.xChristoffelFirst[0][1][1];

  // Compute the derivative of the second kind
  size_t i, j, k;
  for(i = 0; i < 2; i++) for(j = 0; j < 2; j++) for(k = 0; k < 2; k++)
  	dAtom.G.xChristoffelSecond[i][j][k] = 
	  	dAtom.G.xContravariantTensor[k][0] * G.xChristoffelFirst[i][j][0] + 
	  	dAtom.G.xContravariantTensor[k][1] * G.xChristoffelFirst[i][j][1] +	  	 
	  	G.xContravariantTensor[k][0] * dAtom.G.xChristoffelFirst[i][j][0] + 
	  	G.xContravariantTensor[k][1] * dAtom.G.xChristoffelFirst[i][j][1]; 	  	 
}






void AddScaleMedialAtoms(
  const MedialAtom &A, const MedialAtom &B, double p, MedialAtom &C)
{
  // The u and v coordinates stay the same
  C.u = A.u;
  C.v = A.v;
  C.uIndex = A.uIndex;
  C.vIndex = A.vIndex;
  C.flagCrest = A.flagCrest;
  C.flagValid = A.flagValid;

  // The other objects are simply added and scaled
  C.X = A.X + p * B.X;
  C.Xu = A.Xu + p * B.Xu;
  C.Xv = A.Xv + p * B.Xv;
  C.Xuu = A.Xuu + p * B.Xuu;
  C.Xuv = A.Xuv + p * B.Xuv;
  C.Xvv = A.Xvv + p * B.Xvv;

  C.R = A.R + p * B.R;
  C.Ru = A.Ru + p * B.Ru;
  C.Rv = A.Rv + p * B.Rv;
  // C.Ruu = A.Ruu + p * B.Ruu;
  // C.Ruv = A.Ruv + p * B.Ruv;
  // C.Rvv = A.Rvv + p * B.Rvv;

  C.F = A.F + p * B.F;
  C.Fu = A.Fu + p * B.Fu;
  C.Fv = A.Fv + p * B.Fv;

  C.N = A.N + p * B.N;
  C.xGradR = A.xGradR + p * B.xGradR;
  C.Rs2 = A.Rs2 + p * B.Rs2;
  // C.xGradPhi = A.xGradPhi + p * B.xGradPhi;
  C.xGradRMagSqr = A.xGradRMagSqr + p * B.xGradRMagSqr;
  C.xGradRMagSqrOrig = A.xGradRMagSqrOrig + p * B.xGradRMagSqrOrig;
  C.xNormalFactor = A.xNormalFactor + p * B.xNormalFactor;

  // The differential geometry is also added and scaled 
  C.G.g = A.G.g + p * B.G.g;
  C.G.gInv = A.G.gInv + p * B.G.gInv;

  // Curvatures
  C.xMeanCurv = A.xMeanCurv + p * B.xMeanCurv;
  C.xGaussCurv = A.xGaussCurv + p * B.xGaussCurv;

  for(size_t i=0; i<2; i++) for(size_t j=0; j<2; j++)
    {
    C.G.xContravariantTensor[i][j] = 
      A.G.xContravariantTensor[i][j] + p * B.G.xContravariantTensor[i][j];
    C.G.xCovariantTensor[i][j] = 
      A.G.xCovariantTensor[i][j] + p * B.G.xCovariantTensor[i][j];
    C.G.xChristoffelSecond[i][j][0] = 
      A.G.xChristoffelSecond[i][j][0] + p * B.G.xChristoffelSecond[i][j][0];
    C.G.xChristoffelFirst[i][j][0] = 
      A.G.xChristoffelFirst[i][j][0] + p * B.G.xChristoffelFirst[i][j][0];
    C.G.xChristoffelSecond[i][j][1] = 
      A.G.xChristoffelSecond[i][j][1] + p * B.G.xChristoffelSecond[i][j][1];
    C.G.xChristoffelFirst[i][j][1] = 
      A.G.xChristoffelFirst[i][j][1] + p * B.G.xChristoffelFirst[i][j][1];
    }

  // The boundary atoms are scaled as well
  for(size_t k=0; k<2; k++)
    {
    C.xBnd[k].X = A.xBnd[k].X + p * B.xBnd[k].X;
    C.xBnd[k].N = A.xBnd[k].N + p * B.xBnd[k].N;
    }
}

void MedialAtomCentralDifference(
  const MedialAtom &A, const MedialAtom &B, double eps, MedialAtom &C)
{
  double p = 0.5 / eps;
  
  // The u and v coordinates stay the same
  C.u = A.u;
  C.v = A.v;
  C.uIndex = A.uIndex;
  C.vIndex = A.vIndex;
  C.flagCrest = A.flagCrest;
  C.flagValid = A.flagValid;

  // The other objects are simply added and scaled
  C.X = p * (A.X - B.X);
  C.Xu = p * (A.Xu - B.Xu);
  C.Xv = p * (A.Xv - B.Xv);
  C.Xuu = p * (A.Xuu - B.Xuu);
  C.Xuv = p * (A.Xuv - B.Xuv);
  C.Xvv = p * (A.Xvv - B.Xvv);

  C.R = p * (A.R - B.R);
  C.Ru = p * (A.Ru - B.Ru);
  C.Rv = p * (A.Rv - B.Rv);
  // C.Ruu = A.Ruu + p * B.Ruu;
  // C.Ruv = A.Ruv + p * B.Ruv;
  // C.Rvv = A.Rvv + p * B.Rvv;

  C.F = p * (A.F - B.F);
  C.Fu = p * (A.Fu - B.Fu);
  C.Fv = p * (A.Fv - B.Fv);

  C.N = p * (A.N - B.N);
  C.xGradR = p * (A.xGradR - B.xGradR);
  // C.xGradPhi = p * (A.xGradPhi - B.xGradPhi);
  C.xGradRMagSqr = p * (A.xGradRMagSqr - B.xGradRMagSqr);
  C.xGradRMagSqrOrig = p * (A.xGradRMagSqrOrig - B.xGradRMagSqrOrig);
  C.xNormalFactor = p * (A.xNormalFactor - B.xNormalFactor);
  C.Rs2 = p * (A.Rs2 - B.Rs2);

  // The differential geometry is also added and scaled 
  C.G.g = p * (A.G.g - B.G.g);
  C.G.gInv = p * (A.G.gInv - B.G.gInv);

  // Gaussian
  C.xMeanCurv = p * (A.xMeanCurv - B.xMeanCurv);
  C.xGaussCurv = p * (A.xGaussCurv - B.xGaussCurv);

  for(size_t i=0; i<2; i++) for(size_t j=0; j<2; j++)
    {
    C.G.xContravariantTensor[i][j] = 
      p * (A.G.xContravariantTensor[i][j] - B.G.xContravariantTensor[i][j]);
    C.G.xCovariantTensor[i][j] = 
      p * (A.G.xCovariantTensor[i][j] - B.G.xCovariantTensor[i][j]);
    C.G.xChristoffelSecond[i][j][0] = 
      p * (A.G.xChristoffelSecond[i][j][0] - B.G.xChristoffelSecond[i][j][0]);
    C.G.xChristoffelFirst[i][j][0] = 
      p * (A.G.xChristoffelFirst[i][j][0] - B.G.xChristoffelFirst[i][j][0]);
    C.G.xChristoffelSecond[i][j][1] = 
      p * (A.G.xChristoffelSecond[i][j][1] - B.G.xChristoffelSecond[i][j][1]);
    C.G.xChristoffelFirst[i][j][1] = 
      p * (A.G.xChristoffelFirst[i][j][1] - B.G.xChristoffelFirst[i][j][1]);
    }

  // The boundary atoms are scaled as well
  for(size_t k=0; k<2; k++)
    {
    C.xBnd[k].X = p * (A.xBnd[k].X - B.xBnd[k].X);
    C.xBnd[k].N = p * (A.xBnd[k].N - B.xBnd[k].N);
    }
}
