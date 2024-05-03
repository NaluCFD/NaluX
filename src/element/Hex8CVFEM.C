/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <element/Hex8CVFEM.h>
#include <element/Hex8GeometryFunctions.h>
#include <element/FORTRAN_Proto.h>
#include <element/Element.h>
#include <element/ElementFunctions.h>

#include <NaluEnv.h>

#include <cmath>
#include <iostream>

namespace sierra{
namespace nalu{

//-------- hex8_derivative -------------------------------------------------
template <typename DerivType>
void hex8_derivative(
  const int npts,
  const double *intgLoc,
  DerivType &deriv)
{
  const DoubleType half = 0.50;
  const DoubleType one4th = 0.25;
  for (int  ip = 0; ip < npts; ++ip) {
    const DoubleType s1 = intgLoc[ip*3];
    const DoubleType s2 = intgLoc[ip*3+1];
    const DoubleType s3 = intgLoc[ip*3+2];
    const DoubleType s1s2 = s1*s2;
    const DoubleType s2s3 = s2*s3;
    const DoubleType s1s3 = s1*s3;

    // shape function derivative in the s1 direction -
    deriv(ip,0,0) = half*( s3 + s2 ) - s2s3 - one4th;
    deriv(ip,1,0) = half*(-s3 - s2 ) + s2s3 + one4th;
    deriv(ip,2,0) = half*(-s3 + s2 ) - s2s3 + one4th;
    deriv(ip,3,0) = half*( s3 - s2 ) + s2s3 - one4th;
    deriv(ip,4,0) = half*(-s3 + s2 ) + s2s3 - one4th;
    deriv(ip,5,0) = half*( s3 - s2 ) - s2s3 + one4th;
    deriv(ip,6,0) = half*( s3 + s2 ) + s2s3 + one4th;
    deriv(ip,7,0) = half*(-s3 - s2 ) - s2s3 - one4th;

    // shape function derivative in the s2 direction -
    deriv(ip,0,1) = half*( s3 + s1 ) - s1s3 - one4th;
    deriv(ip,1,1) = half*( s3 - s1 ) + s1s3 - one4th;
    deriv(ip,2,1) = half*(-s3 + s1 ) - s1s3 + one4th;
    deriv(ip,3,1) = half*(-s3 - s1 ) + s1s3 + one4th;
    deriv(ip,4,1) = half*(-s3 + s1 ) + s1s3 - one4th;
    deriv(ip,5,1) = half*(-s3 - s1 ) - s1s3 - one4th;
    deriv(ip,6,1) = half*( s3 + s1 ) + s1s3 + one4th;
    deriv(ip,7,1) = half*( s3 - s1 ) - s1s3 + one4th;

    // shape function derivative in the s3 direction -
    deriv(ip,0,2) = half*( s2 + s1 ) - s1s2 - one4th;
    deriv(ip,1,2) = half*( s2 - s1 ) + s1s2 - one4th;
    deriv(ip,2,2) = half*(-s2 - s1 ) - s1s2 - one4th;
    deriv(ip,3,2) = half*(-s2 + s1 ) + s1s2 - one4th;
    deriv(ip,4,2) = half*(-s2 - s1 ) + s1s2 + one4th;
    deriv(ip,5,2) = half*(-s2 + s1 ) - s1s2 + one4th;
    deriv(ip,6,2) = half*( s2 + s1 ) + s1s2 + one4th;
    deriv(ip,7,2) = half*( s2 - s1 ) - s1s2 + one4th;
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HexSCV::HexSCV()
  : Element()
{
  nDim_ = 3;
  nodesPerElement_ = 8;
  numIntPoints_ = 8;

  // define ip node mappings
  ipNodeMap_.resize(8);
  ipNodeMap_[0] = 0; ipNodeMap_[1] = 1; ipNodeMap_[2] = 2; ipNodeMap_[3] = 3;
  ipNodeMap_[4] = 4; ipNodeMap_[5] = 5; ipNodeMap_[6] = 6; ipNodeMap_[7] = 7;

  // standard integration location
  intgLoc_.resize(24);
  intgLoc_[0]  = -0.25; intgLoc_[1]  = -0.25; intgLoc_[2]  = -0.25;
  intgLoc_[3]  = +0.25; intgLoc_[4]  = -0.25; intgLoc_[5]  = -0.25;
  intgLoc_[6]  = +0.25; intgLoc_[7]  = +0.25; intgLoc_[8]  = -0.25;
  intgLoc_[9]  = -0.25; intgLoc_[10] = +0.25; intgLoc_[11] = -0.25;
  intgLoc_[12] = -0.25; intgLoc_[13] = -0.25; intgLoc_[14] = +0.25;
  intgLoc_[15] = +0.25; intgLoc_[16] = -0.25; intgLoc_[17] = +0.25;
  intgLoc_[18] = +0.25; intgLoc_[19] = +0.25; intgLoc_[20] = +0.25;
  intgLoc_[21] = -0.25; intgLoc_[22] = +0.25; intgLoc_[23] = +0.25;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
HexSCV::~HexSCV()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
HexSCV::ipNodeMap(
  int /*ordinal*/)
{
  // define scv->node mappings
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void HexSCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(hex_scv_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords,
      volume, error, &lerr );
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void HexSCV::determinant(
  SharedMemView<DoubleType**>& coords,
  SharedMemView<DoubleType*>& volume)
{
  constexpr int subDivisionTable[8][8] = {
      {  0,  8, 12, 11, 19, 20, 26, 25},
      {  8,  1,  9, 12, 20, 18, 24, 26},
      { 12,  9,  2, 10, 26, 24, 22, 23},
      { 11, 12, 10,  3, 25, 26, 23, 21},
      { 19, 20, 26, 25,  4, 13, 17, 16},
      { 20, 18, 24, 26, 13,  5, 14, 17},
      { 26, 24, 22, 23, 17, 14,  6, 15},
      { 25, 26, 23, 21, 16, 17, 15,  7}
  };

  DoubleType coordv[27][3];
  subdivide_hex_8(coords, coordv);

  constexpr int numSCV = 8;
  for (int ip = 0; ip < numSCV; ++ip) {
    DoubleType scvHex[8][3];
    for (int n = 0; n < 8; ++n) {
      const int subIndex = subDivisionTable[ip][n];
      for (int d = 0; d < 3; ++d) {
        scvHex[n][d] = coordv[subIndex][d];
      }
    }
    volume(ip) = hex_volume_grandy(scvHex);
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HexSCS::HexSCS()
  : Element()
{
  nDim_ = 3;
  nodesPerElement_ = 8;
  numIntPoints_ = 12;
  scaleToStandardIsoFac_ = 2.0;

  // define L/R mappings
  lrscv_.resize(24);
  lrscv_[0]  = 0; lrscv_[1]  = 1;
  lrscv_[2]  = 1; lrscv_[3]  = 2;
  lrscv_[4]  = 2; lrscv_[5]  = 3;
  lrscv_[6]  = 0; lrscv_[7]  = 3;
  lrscv_[8]  = 4; lrscv_[9]  = 5;
  lrscv_[10] = 5; lrscv_[11] = 6;
  lrscv_[12] = 6; lrscv_[13] = 7;
  lrscv_[14] = 4; lrscv_[15] = 7;
  lrscv_[16] = 0; lrscv_[17] = 4;
  lrscv_[18] = 1; lrscv_[19] = 5;
  lrscv_[20] = 2; lrscv_[21] = 6;
  lrscv_[22] = 3; lrscv_[23] = 7;

    // boundary integration point ip node mapping (ip on an ordinal to local node number)
  ipNodeMap_.resize(24); // 4 ips * 6 faces
  // face 0;
  ipNodeMap_[0] = 0;  ipNodeMap_[1] = 1;  ipNodeMap_[2] = 5;  ipNodeMap_[3] = 4;
  // face 1;
  ipNodeMap_[4] = 1;  ipNodeMap_[5] = 2;  ipNodeMap_[6] = 6;  ipNodeMap_[7] = 5;
  // face 2;
  ipNodeMap_[8] = 2;  ipNodeMap_[9] = 3;  ipNodeMap_[10] = 7; ipNodeMap_[11] = 6;
  // face 3;
  ipNodeMap_[12] = 0; ipNodeMap_[13] = 4; ipNodeMap_[14] = 7; ipNodeMap_[15] = 3;
  // face 4;
  ipNodeMap_[16] = 0; ipNodeMap_[17] = 3; ipNodeMap_[18] = 2; ipNodeMap_[19] = 1;
  // face 5;
  ipNodeMap_[20] = 4; ipNodeMap_[21] = 5; ipNodeMap_[22] = 6; ipNodeMap_[23] = 7;

  // standard integration location
  intgLoc_.resize(36);
  intgLoc_[0]  =  0.00; intgLoc_[1]  = -0.25; intgLoc_[2]  = -0.25; // surf 1    1->2
  intgLoc_[3]  =  0.25; intgLoc_[4]  =  0.00; intgLoc_[5]  = -0.25; // surf 2    2->3
  intgLoc_[6]  =  0.00; intgLoc_[7]  =  0.25; intgLoc_[8]  = -0.25; // surf 3    3->4
  intgLoc_[9]  = -0.25; intgLoc_[10] =  0.00; intgLoc_[11] = -0.25; // surf 4    1->4
  intgLoc_[12] =  0.00; intgLoc_[13] = -0.25; intgLoc_[14] =  0.25; // surf 5    5->6
  intgLoc_[15] =  0.25; intgLoc_[16] =  0.00; intgLoc_[17] =  0.25; // surf 6    6->7
  intgLoc_[18] =  0.00; intgLoc_[19] =  0.25; intgLoc_[20] =  0.25; // surf 7    7->8
  intgLoc_[21] = -0.25; intgLoc_[22] =  0.00; intgLoc_[23] =  0.25; // surf 8    5->8
  intgLoc_[24] = -0.25; intgLoc_[25] = -0.25; intgLoc_[26] =  0.00; // surf 9    1->5
  intgLoc_[27] =  0.25; intgLoc_[28] = -0.25; intgLoc_[29] =  0.00; // surf 10   2->6
  intgLoc_[30] =  0.25; intgLoc_[31] =  0.25; intgLoc_[32] =  0.00; // surf 11   3->7
  intgLoc_[33] = -0.25; intgLoc_[34] =  0.25; intgLoc_[35] =  0.00; // surf 12   4->8

  // nodes for collocation calculations
  nodeLoc_.resize(24);
  // node 0
  nodeLoc_[0] = -0.5; nodeLoc_[1] = -0.5; nodeLoc_[2] = -0.5;
  // node 1
  nodeLoc_[3] =  0.5; nodeLoc_[4] = -0.5; nodeLoc_[5] = -0.5;
  // node 2
  nodeLoc_[6] =  0.5; nodeLoc_[7] =  0.5; nodeLoc_[8] = -0.5;
  // node 3
  nodeLoc_[9] = -0.5; nodeLoc_[10] =  0.5; nodeLoc_[11] = -0.5;
  // node 4
  nodeLoc_[12] = -0.5; nodeLoc_[13] = -0.5; nodeLoc_[14] =  0.5;
  // node 5
  nodeLoc_[15] =  0.5; nodeLoc_[16] = -0.5; nodeLoc_[17] =  0.5;
  // node 6
  nodeLoc_[18] =  0.5; nodeLoc_[19] =  0.5; nodeLoc_[20] =  0.5;
  // node 7
  nodeLoc_[21] = -0.5; nodeLoc_[22] =  0.5; nodeLoc_[23] =  0.5;

  // mapping from a side ordinal to the node ordinals on that side
  sideNodeOrdinals_ = {
      0, 1, 5, 4, // ordinal 0
      1, 2, 6, 5, // ordinal 1
      2, 3, 7, 6, // ordinal 2
      0, 4, 7, 3, // ordinal 3
      0, 3, 2, 1, // ordinal 4
      4, 5, 6, 7  // ordinal 5
  };
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
HexSCS::~HexSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
HexSCS::ipNodeMap(
  int ordinal)
{
  // define ip->node mappings for each face (ordinal);
  return &ipNodeMap_[ordinal*4];
}


//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
HexSCS::shape_fcn(SharedMemView<DoubleType**> &shpfc)
{
  hex8_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- hex8_shape_fcn --------------------------------------------------
//--------------------------------------------------------------------------
void
HexSCS::hex8_shape_fcn(
  const int  &npts,
  const double *isoParCoord,
  SharedMemView<DoubleType**> &shape_fcn)
{
  const DoubleType half = 0.50;
  const DoubleType one4th = 0.25;
  const DoubleType one8th = 0.125;
  for ( int j = 0; j < numIntPoints_; ++j ) {

    const DoubleType s1 = isoParCoord[j*3];
    const DoubleType s2 = isoParCoord[j*3+1];
    const DoubleType s3 = isoParCoord[j*3+2];

    shape_fcn(j,0) = one8th + one4th*(-s1 - s2 - s3)
      + half*( s2*s3 + s3*s1 + s1*s2 ) - s1*s2*s3;
    shape_fcn(j,1) = one8th + one4th*( s1 - s2 - s3)
      + half*( s2*s3 - s3*s1 - s1*s2 ) + s1*s2*s3;
    shape_fcn(j,2) = one8th + one4th*( s1 + s2 - s3)
      + half*(-s2*s3 - s3*s1 + s1*s2 ) - s1*s2*s3;
    shape_fcn(j,3) = one8th + one4th*(-s1 + s2 - s3)
      + half*(-s2*s3 + s3*s1 - s1*s2 ) + s1*s2*s3;
    shape_fcn(j,4) = one8th + one4th*(-s1 - s2 + s3)
      + half*(-s2*s3 - s3*s1 + s1*s2 ) + s1*s2*s3;
    shape_fcn(j,5) = one8th + one4th*( s1 - s2 + s3)
      + half*(-s2*s3 + s3*s1 - s1*s2 ) - s1*s2*s3;
    shape_fcn(j,6) = one8th + one4th*( s1 + s2 + s3)
      + half*( s2*s3 + s3*s1 + s1*s2 ) + s1*s2*s3;
    shape_fcn(j,7) = one8th + one4th*(-s1 + s2 + s3)
      + half*( s2*s3 - s3*s1 - s1*s2 ) - s1*s2*s3;
  }
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void HexSCS::grad_op(
  SharedMemView<DoubleType**>&coords,
  SharedMemView<DoubleType***>&gradop,
  SharedMemView<DoubleType***>&deriv)
{
  hex8_derivative(numIntPoints_, &intgLoc_[0], deriv);
  generic_grad_op<AlgTraitsHex8>(deriv, coords, gradop);
 }

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void HexSCS::determinant(
  SharedMemView<DoubleType**>&coords,
  SharedMemView<DoubleType**>&areav)
{
  constexpr int hex_edge_facet_table[12][4] = {
      { 20,  8, 12, 26 },
      { 24,  9, 12, 26 },
      { 10, 12, 26, 23 },
      { 11, 25, 26, 12 },
      { 13, 20, 26, 17 },
      { 17, 14, 24, 26 },
      { 17, 15, 23, 26 },
      { 16, 17, 26, 25 },
      { 19, 20, 26, 25 },
      { 20, 18, 24, 26 },
      { 22, 23, 26, 24 },
      { 21, 25, 26, 23 }
  };

  DoubleType coordv[27][3];
  subdivide_hex_8(coords, coordv);

  constexpr int npf = 4;
  constexpr int nscs = 12;
  for (int ics=0; ics < nscs; ++ics) {
    DoubleType scscoords[4][3];
    for (int inode = 0; inode < npf; ++inode) {
      const int itrianglenode = hex_edge_facet_table[ics][inode];
      for (int d=0; d < 3; ++d) {
        scscoords[inode][d] = coordv[itrianglenode][d];
      }
    }
    quad_area_by_triangulation(ics, scscoords, areav);
  }
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void HexSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  SIERRA_FORTRAN(hex_scs_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords, areav );

  // all is always well; no error checking
  *error = 0;
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void HexSCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;
  
  SIERRA_FORTRAN(hex_derivative)
    ( &numIntPoints_,
      &intgLoc_[0], deriv );

  SIERRA_FORTRAN(hex_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    NaluEnv::self().naluOutput() << "sorry, negative HexSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- adjacentNodes ---------------------------------------------------
//--------------------------------------------------------------------------
const int *
HexSCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
HexSCS::shape_fcn(double *shpfc)
{
  SIERRA_FORTRAN(hex_shape_fcn)
    (&numIntPoints_,&intgLoc_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
HexSCS::isInElement(
    const double * elem_nodal_coor,     // (8,3)
    const double * point_coor,          // (3)
    double * par_coor )
{
  const int maxNonlinearIter = 20;
  const double isInElemConverged = 1.0e-16;
  // Translate element so that (x,y,z) coordinates of the first node are (0,0,0)

  double x[] = {0.,
        0.125*(elem_nodal_coor[1] - elem_nodal_coor[0]),
        0.125*(elem_nodal_coor[2] - elem_nodal_coor[0]),
        0.125*(elem_nodal_coor[3] - elem_nodal_coor[0]),
        0.125*(elem_nodal_coor[4] - elem_nodal_coor[0]),
        0.125*(elem_nodal_coor[5] - elem_nodal_coor[0]),
        0.125*(elem_nodal_coor[6] - elem_nodal_coor[0]),
        0.125*(elem_nodal_coor[7] - elem_nodal_coor[0]) };
  double y[] = {0.,
        0.125*(elem_nodal_coor[9 ] - elem_nodal_coor[8]),
        0.125*(elem_nodal_coor[10] - elem_nodal_coor[8]),
        0.125*(elem_nodal_coor[11] - elem_nodal_coor[8]),
        0.125*(elem_nodal_coor[12] - elem_nodal_coor[8]),
        0.125*(elem_nodal_coor[13] - elem_nodal_coor[8]),
        0.125*(elem_nodal_coor[14] - elem_nodal_coor[8]),
        0.125*(elem_nodal_coor[15] - elem_nodal_coor[8]) };
  double z[] = {0.,
        0.125*(elem_nodal_coor[17] - elem_nodal_coor[16]),
        0.125*(elem_nodal_coor[18] - elem_nodal_coor[16]),
        0.125*(elem_nodal_coor[19] - elem_nodal_coor[16]),
        0.125*(elem_nodal_coor[20] - elem_nodal_coor[16]),
        0.125*(elem_nodal_coor[21] - elem_nodal_coor[16]),
        0.125*(elem_nodal_coor[22] - elem_nodal_coor[16]),
        0.125*(elem_nodal_coor[23] - elem_nodal_coor[16]) };

  // (xp,yp,zp) is the point at which we're searching for (xi,eta,zeta)
  // (must translate this also)

  double xp = point_coor[0] - elem_nodal_coor[0];
  double yp = point_coor[1] - elem_nodal_coor[8];
  double zp = point_coor[2] - elem_nodal_coor[16];

  // Newton-Raphson iteration for (xi,eta,zeta)
  double j[9];
  double f[3];
  double shapefct[8];
  double xinew = 0.5;     // initial guess
  double etanew = 0.5;
  double zetanew = 0.5;
  double xicur = 0.5;
  double etacur = 0.5;
  double zetacur = 0.5;
  double xidiff[] = { 1.0, 1.0, 1.0 };
  int i = 0;

  do
  {
    j[0]=
      -(1.0-etacur)*(1.0-zetacur)*x[1]
      -(1.0+etacur)*(1.0-zetacur)*x[2]
      +(1.0+etacur)*(1.0-zetacur)*x[3]
      +(1.0-etacur)*(1.0+zetacur)*x[4]
      -(1.0-etacur)*(1.0+zetacur)*x[5]
      -(1.0+etacur)*(1.0+zetacur)*x[6]
      +(1.0+etacur)*(1.0+zetacur)*x[7];

    j[1]=
       (1.0+xicur)*(1.0-zetacur)*x[1]
      -(1.0+xicur)*(1.0-zetacur)*x[2]
      -(1.0-xicur)*(1.0-zetacur)*x[3]
      +(1.0-xicur)*(1.0+zetacur)*x[4]
      +(1.0+xicur)*(1.0+zetacur)*x[5]
      -(1.0+xicur)*(1.0+zetacur)*x[6]
      -(1.0-xicur)*(1.0+zetacur)*x[7];

    j[2]=
       (1.0-etacur)*(1.0+xicur)*x[1]
      +(1.0+etacur)*(1.0+xicur)*x[2]
      +(1.0+etacur)*(1.0-xicur)*x[3]
      -(1.0-etacur)*(1.0-xicur)*x[4]
      -(1.0-etacur)*(1.0+xicur)*x[5]
      -(1.0+etacur)*(1.0+xicur)*x[6]
      -(1.0+etacur)*(1.0-xicur)*x[7];

    j[3]=
      -(1.0-etacur)*(1.0-zetacur)*y[1]
      -(1.0+etacur)*(1.0-zetacur)*y[2]
      +(1.0+etacur)*(1.0-zetacur)*y[3]
      +(1.0-etacur)*(1.0+zetacur)*y[4]
      -(1.0-etacur)*(1.0+zetacur)*y[5]
      -(1.0+etacur)*(1.0+zetacur)*y[6]
      +(1.0+etacur)*(1.0+zetacur)*y[7];

    j[4]=
       (1.0+xicur)*(1.0-zetacur)*y[1]
      -(1.0+xicur)*(1.0-zetacur)*y[2]
      -(1.0-xicur)*(1.0-zetacur)*y[3]
      +(1.0-xicur)*(1.0+zetacur)*y[4]
      +(1.0+xicur)*(1.0+zetacur)*y[5]
      -(1.0+xicur)*(1.0+zetacur)*y[6]
      -(1.0-xicur)*(1.0+zetacur)*y[7];

    j[5]=
       (1.0-etacur)*(1.0+xicur)*y[1]
      +(1.0+etacur)*(1.0+xicur)*y[2]
      +(1.0+etacur)*(1.0-xicur)*y[3]
      -(1.0-etacur)*(1.0-xicur)*y[4]
      -(1.0-etacur)*(1.0+xicur)*y[5]
      -(1.0+etacur)*(1.0+xicur)*y[6]
      -(1.0+etacur)*(1.0-xicur)*y[7];

    j[6]=
      -(1.0-etacur)*(1.0-zetacur)*z[1]
      -(1.0+etacur)*(1.0-zetacur)*z[2]
      +(1.0+etacur)*(1.0-zetacur)*z[3]
      +(1.0-etacur)*(1.0+zetacur)*z[4]
      -(1.0-etacur)*(1.0+zetacur)*z[5]
      -(1.0+etacur)*(1.0+zetacur)*z[6]
      +(1.0+etacur)*(1.0+zetacur)*z[7];

    j[7]=
       (1.0+xicur)*(1.0-zetacur)*z[1]
      -(1.0+xicur)*(1.0-zetacur)*z[2]
      -(1.0-xicur)*(1.0-zetacur)*z[3]
      +(1.0-xicur)*(1.0+zetacur)*z[4]
      +(1.0+xicur)*(1.0+zetacur)*z[5]
      -(1.0+xicur)*(1.0+zetacur)*z[6]
      -(1.0-xicur)*(1.0+zetacur)*z[7];

    j[8]=
       (1.0-etacur)*(1.0+xicur)*z[1]
      +(1.0+etacur)*(1.0+xicur)*z[2]
      +(1.0+etacur)*(1.0-xicur)*z[3]
      -(1.0-etacur)*(1.0-xicur)*z[4]
      -(1.0-etacur)*(1.0+xicur)*z[5]
      -(1.0+etacur)*(1.0+xicur)*z[6]
      -(1.0+etacur)*(1.0-xicur)*z[7];

    double jdet=-(j[2]*j[4]*j[6])+j[1]*j[5]*j[6]+j[2]*j[3]*j[7]-
      j[0]*j[5]*j[7]-j[1]*j[3]*j[8]+j[0]*j[4]*j[8];

    if (!jdet) {
      i = maxNonlinearIter;
      break;
    }
    shapefct[0]=(1.0-etacur)*(1.0-xicur)*(1.0-zetacur);

    shapefct[1]=(1.0-etacur)*(1.0+xicur)*(1.0-zetacur);

    shapefct[2]=(1.0+etacur)*(1.0+xicur)*(1.0-zetacur);

    shapefct[3]=(1.0+etacur)*(1.0-xicur)*(1.0-zetacur);

    shapefct[4]=(1.0-etacur)*(1.0-xicur)*(1.0+zetacur);

    shapefct[5]=(1.0-etacur)*(1.0+xicur)*(1.0+zetacur);

    shapefct[6]=(1.0+etacur)*(1.0+xicur)*(1.0+zetacur);

    shapefct[7]=(1.0+etacur)*(1.0-xicur)*(1.0+zetacur);

    f[0]=xp-shapefct[1]*x[1]-shapefct[2]*x[2]-shapefct[3]*x[3]-shapefct[4]*x[4]-\
      shapefct[5]*x[5]-shapefct[6]*x[6]-shapefct[7]*x[7];

    f[1]=yp-shapefct[1]*y[1]-shapefct[2]*y[2]-shapefct[3]*y[3]-shapefct[4]*y[4]-\
      shapefct[5]*y[5]-shapefct[6]*y[6]-shapefct[7]*y[7];

    f[2]=zp-shapefct[1]*z[1]-shapefct[2]*z[2]-shapefct[3]*z[3]-shapefct[4]*z[4]-\
      shapefct[5]*z[5]-shapefct[6]*z[6]-shapefct[7]*z[7];

    xinew = (jdet*xicur+f[2]*(j[2]*j[4]-j[1]*j[5])-f[1]*j[2]*j[7]+f[0]*j[5]*j[7]+
       f[1]*j[1]*j[8]-f[0]*j[4]*j[8])/jdet;

    etanew = (etacur*jdet+f[2]*(-(j[2]*j[3])+j[0]*j[5])+f[1]*j[2]*j[6]-f[0]*j[5]*j[6]-
        f[1]*j[0]*j[8]+f[0]*j[3]*j[8])/jdet;

    zetanew = (jdet*zetacur+f[2]*(j[1]*j[3]-j[0]*j[4])-f[1]*j[1]*j[6]+
         f[0]*j[4]*j[6]+f[1]*j[0]*j[7]-f[0]*j[3]*j[7])/jdet;

    xidiff[0] = xinew - xicur;
    xidiff[1] = etanew - etacur;
    xidiff[2] = zetanew - zetacur;
    xicur = xinew;
    etacur = etanew;
    zetacur = zetanew;

  }
  while ( !within_tolerance( vector_norm_sq(xidiff,3), isInElemConverged) && ++i < maxNonlinearIter);

  par_coor[0] = par_coor[1] = par_coor[2] = std::numeric_limits<double>::max();
  double dist = std::numeric_limits<double>::max();

  if (i <maxNonlinearIter) {
    par_coor[0] = xinew;
    par_coor[1] = etanew;
    par_coor[2] = zetanew;

    std::vector<double> xtmp(3);
    xtmp[0] = par_coor[0];
    xtmp[1] = par_coor[1];
    xtmp[2] = par_coor[2];
    dist = parametric_distance(xtmp);
  }
  return dist;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
HexSCS::interpolatePoint(
    const int  & ncomp_field,
    const double * par_coord,           // (3)
    const double * field,               // (8,ncomp_field)
    double * result ) // (ncomp_field)
{
  // 'field' is a flat array of dimension (8,ncomp_field) (Fortran ordering);
  double xi   = par_coord[0];
  double eta  = par_coord[1];
  double zeta = par_coord[2];

  // NOTE: this uses a [-1,1] definition of the reference element,
  // contrary to the rest of the code

  for ( int i = 0; i < ncomp_field; i++ )
  {
    // Base 'field array' index for ith component
    int b = 8*i;

    result[i] = 0.125 * (
        (1 - xi) * (1 - eta) * (1 - zeta) * field[b + 0]
      + (1 + xi) * (1 - eta) * (1 - zeta) * field[b + 1]
      + (1 + xi) * (1 + eta) * (1 - zeta) * field[b + 2]
      + (1 - xi) * (1 + eta) * (1 - zeta) * field[b + 3]
      + (1 - xi) * (1 - eta) * (1 + zeta) * field[b + 4]
      + (1 + xi) * (1 - eta) * (1 + zeta) * field[b + 5]
      + (1 + xi) * (1 + eta) * (1 + zeta) * field[b + 6]
      + (1 - xi) * (1 + eta) * (1 + zeta) * field[b + 7]
    );
  }
}

//--------------------------------------------------------------------------
//-------- parametric_distance ---------------------------------------------
//--------------------------------------------------------------------------
double HexSCS::parametric_distance(const std::vector<double> &x)
{
  std::vector<double> y(3);
  for (int i=0; i<3; ++i) {
    y[i] = std::fabs(x[i]);
  }

  double d = 0;
  for (int i=0; i<3; ++i) {
    if (d < y[i]) {
      d = y[i];
    }
  }
  return d;
}

}
}
