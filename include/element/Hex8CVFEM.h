/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Hex8CVFEM_h
#define Hex8CVFEM_h

#include<AlgTraits.h>
#include<element/Element.h>

namespace sierra{
namespace nalu{

// Hex 8 subcontrol volume
class HexSCV : public Element
{
public:
  using AlgTraits = AlgTraitsHex8;

  HexSCV();
  virtual ~HexSCV();

  const int * ipNodeMap(int ordinal = 0);

  // NGP-ready methods first
  void determinant(
    SharedMemView<DoubleType**>& coords,
    SharedMemView<DoubleType*>& volume);

  // Non-NGPS-ready methods second
  void determinant(
    const int nelem,
    const double *coords,
    double *volume,
    double * error );
};

// Hex 8 subcontrol surface
class HexSCS : public Element
{
public:
  using AlgTraits = AlgTraitsHex8;
  using AlgTraitsFace = AlgTraitsQuad4;

  HexSCS();
  virtual ~HexSCS();

  const int * ipNodeMap(int ordinal = 0);
  const int * adjacentNodes();

  // NGP-ready methods first
  void shape_fcn(
    SharedMemView<DoubleType**> &shpfc);

  void grad_op(
    SharedMemView<DoubleType**>&coords,
    SharedMemView<DoubleType***>&gradop,
    SharedMemView<DoubleType***>&deriv);

  void determinant(
    SharedMemView<DoubleType**>&coords,
    SharedMemView<DoubleType**>&areav);

  void hex8_gradient_operator(
    const int nodesPerElem,
    const int numIntgPts,
    SharedMemView<DoubleType***> &deriv,
    SharedMemView<DoubleType**> &cordel,
    SharedMemView<DoubleType***> &gradop,
    SharedMemView<DoubleType*> &det_j,
    DoubleType &error,
    int &lerr);

  // NGP helper
  void hex8_shape_fcn(
    const int  &numIp,
    const double *isoParCoord,
    SharedMemView<DoubleType**> &shpfc);

  // non NGP-ready methods second
  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shape_fcn(
    double *shpfc);

  // for transfer
  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  double parametric_distance(const std::vector<double> &x);
};
    
} // namespace nalu
} // namespace Sierra

#endif
