/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Quad43DSCS_h
#define Quad43DSCS_h

#include "element/Element.h"

namespace sierra{
namespace nalu{

// 3D Quad 4
class Quad3DSCS : public Element
{
public:

  Quad3DSCS();
  virtual ~Quad3DSCS();

  const int * ipNodeMap(int ordinal = 0);
 
  // NGP-ready methods first
  void shape_fcn(
    SharedMemView<DoubleType**> &shpfc);

  void quad4_shape_fcn(
    const int  &numIp,
    const double *isoParCoord,
    SharedMemView<DoubleType**> &shpfc);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  const double elemThickness_;
};
    
} // namespace nalu
} // namespace Sierra

#endif
