/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Element_h
#define Element_h

#include <element/ElementFactory.h>

// NGP-based includes
#include "SimdInterface.h"
#include "KokkosInterface.h"

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <array>

namespace stk {
  struct topology;
}

namespace sierra{
namespace nalu{

namespace MEconstants {
  static const double realmin = std::numeric_limits<double>::min();
}

namespace Jacobian{
enum Direction
{
  S_DIRECTION = 0,
  T_DIRECTION = 1,
  U_DIRECTION = 2
};
}

class Element
{
public:
  Element();
  virtual ~Element();

  // NGP-ready methods first
  virtual void determinant(
    SharedMemView<DoubleType**>&coords,
    SharedMemView<DoubleType**>&areav) {
    throw std::runtime_error("determinant using SharedMemView is not implemented");}

  virtual void grad_op(
    SharedMemView<DoubleType**>&coords,
    SharedMemView<DoubleType***>&gradop,
    SharedMemView<DoubleType***>&deriv) {
    throw std::runtime_error("grad_op using SharedMemView is not implemented");}

  virtual void shape_fcn(
    SharedMemView<DoubleType**> &shpfc) {
    throw std::runtime_error("shape_fcn using SharedMemView is not implemented");}

  // non-NGP-ready methods second
  virtual void determinant(
    const int nelem,
    const double *coords,
    double *volume,
    double * error ) {
    throw std::runtime_error("determinant not implemented");}

  virtual void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error ) {
    throw std::runtime_error("grad_op not implemented");}

  virtual void shape_fcn(
    double *shpfc) {
    throw std::runtime_error("shape_fcn not implemented"); }

  virtual const int * adjacentNodes() {
    throw std::runtime_error("adjacentNodes not implemented");
    return NULL;}

  virtual const int * ipNodeMap(int ordinal = 0) {
      throw std::runtime_error("ipNodeMap not implemented");
      return NULL;}

  // for transfer
  virtual double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord) {
    throw std::runtime_error("isInElement not implemented"); 
    return 1.0e6; }

  virtual void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result) {
    throw std::runtime_error("interpolatePoint not implemented"); }

  //double parametric_distance(const std::array<double, 2>& x);
  bool within_tolerance(const double & val, const double & tol);
  double vector_norm_sq(const double * vect, int len);

  int nDim_;
  int nodesPerElement_;
  int numIntPoints_;
  double scaleToStandardIsoFac_ = 2.0;

  std::vector<int> lrscv_;
  std::vector<int> ipNodeMap_;
  std::vector<double> intgLoc_;
  std::vector<double> nodeLoc_;
  std::vector<int> sideNodeOrdinals_;
  std::vector<int> sideOffset_;

  // FEM
  std::vector<double>weights_;
};


} // namespace nalu
} // namespace Sierra

#endif
