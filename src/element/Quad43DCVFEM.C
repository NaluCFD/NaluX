/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "element/Element.h"
#include "element/Quad43DCVFEM.h"

#include <element/FORTRAN_Proto.h>

namespace sierra{
namespace nalu{


//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Quad3DSCS::Quad3DSCS()  
  : Element(),
    elemThickness_(0.1)
{
  nDim_ = 3;
  nodesPerElement_ = 4;
  numIntPoints_ = 4;
  scaleToStandardIsoFac_ = 2.0;

  // define ip node mappings; ordinal size = 1
  ipNodeMap_.resize(4);
  ipNodeMap_[0] = 0;
  ipNodeMap_[1] = 1;
  ipNodeMap_[2] = 2;
  ipNodeMap_[3] = 3;

  // standard integration location
  intgLoc_.resize(8);    
  intgLoc_[0]  = -0.25; intgLoc_[1] = -0.25; // surf 1
  intgLoc_[2]  =  0.25; intgLoc_[3] = -0.25; // surf 2
  intgLoc_[4]  =  0.25; intgLoc_[5] =  0.25; // surf 3
  intgLoc_[6]  = -0.25; intgLoc_[7] =  0.25; // surf 4
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Quad3DSCS::~Quad3DSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad3DSCS::ipNodeMap(
  int /*ordinal*/)
{
  // define ip->node mappings for each face (single ordinal); 
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Quad3DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(quad3d_scs_det)
    ( &nelem, coords, areav );

  // fake check
  *error = (lerr == 0) ? 0.0 : 1.0;
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Quad3DSCS::shape_fcn(double *shpfc)
{
  SIERRA_FORTRAN(quad3d_shape_fcn)
    (&numIntPoints_,&intgLoc_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Quad3DSCS::shape_fcn(SharedMemView<DoubleType**> &shpfc)
{
  quad4_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- quad4_shape_fcn --------------------------------------------------
//--------------------------------------------------------------------------
void
Quad3DSCS::quad4_shape_fcn(
  const int  &npts,
  const double *isoParCoord,
  SharedMemView<DoubleType**> &shape_fcn)
{
  // -1/2:1/2 isoparametric range
  const DoubleType half = 0.50;
  const DoubleType one4th = 0.25;
  for ( int j = 0; j < numIntPoints_; ++j ) {

    const DoubleType s1 = isoParCoord[j*2];
    const DoubleType s2 = isoParCoord[j*2+1];

    shape_fcn(j,0) = one4th + half*(-s1 - s2) + s1*s2;
    shape_fcn(j,1) = one4th + half*( s1 - s2) - s1*s2;
    shape_fcn(j,2) = one4th + half*( s1 + s2) + s1*s2;
    shape_fcn(j,3) = one4th + half*(-s1 + s2) - s1*s2;
  }
}

} // namespace nalu
} // namespace sierra
