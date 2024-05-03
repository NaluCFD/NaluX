/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ElementFunctions_h
#define ElementFunctions_h

#include <AlgTraits.h>

#include <element/Element.h>

#include <SimdInterface.h>
#include <Kokkos_Core.hpp>

#include <stk_util/util/ReportHandler.hpp>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <array>
#include <type_traits>

namespace sierra {
namespace nalu {

  template<typename ftype> inline void cofactorMatrix(ftype adjJac[][3], const ftype jact[][3]) {
    adjJac[0][0] = jact[1][1] * jact[2][2] - jact[2][1] * jact[1][2];
    adjJac[0][1] = jact[1][2] * jact[2][0] - jact[2][2] * jact[1][0];
    adjJac[0][2] = jact[1][0] * jact[2][1] - jact[2][0] * jact[1][1];

    adjJac[1][0] = jact[0][2] * jact[2][1] - jact[2][2] * jact[0][1];
    adjJac[1][1] = jact[0][0] * jact[2][2] - jact[2][0] * jact[0][2];
    adjJac[1][2] = jact[0][1] * jact[2][0] - jact[2][1] * jact[0][0];

    adjJac[2][0] = jact[0][1] * jact[1][2] - jact[1][1] * jact[0][2];
    adjJac[2][1] = jact[0][2] * jact[1][0] - jact[1][2] * jact[0][0];
    adjJac[2][2] = jact[0][0] * jact[1][1] - jact[1][0] * jact[0][1];
  }
  template<typename ftype> inline void cofactorMatrix(ftype adjJac[][2], const ftype jact[][2]) {
    adjJac[0][0] =  jact[1][1];
    adjJac[0][1] = -jact[1][0];
    adjJac[1][0] = -jact[0][1];
    adjJac[1][1] =  jact[0][0];
  }

  template <typename AlgTraits, typename GradViewType, typename CoordViewType, typename OutputViewType>
  void generic_grad_op(const GradViewType& referenceGradWeights, const CoordViewType& coords, OutputViewType& weights)
  {
    constexpr int dim = AlgTraits::nDim_;
    
    using ftype = typename CoordViewType::value_type;
    static_assert(std::is_same<ftype, typename GradViewType::value_type>::value,  "Incompatiable value type for views");
    static_assert(std::is_same<ftype, typename OutputViewType::value_type>::value,  "Incompatiable value type for views");
    static_assert(GradViewType::rank   ==   3, "grad view assumed to be rank 3");
    static_assert(CoordViewType::rank  ==   2, "Coordinate view assumed to be rank 2");
    static_assert(OutputViewType::rank ==   3, "Weight view assumed to be rank 3");

    STK_ThrowAssert(AlgTraits::nodesPerElement_ == referenceGradWeights.extent(1));
    STK_ThrowAssert(AlgTraits::nDim_            == referenceGradWeights.extent(2));
    for (int i=0; i<dim; ++i) 
      STK_ThrowAssert(weights.extent(i) == referenceGradWeights.extent(i));

    for (unsigned ip = 0; ip < referenceGradWeights.extent(0); ++ip) {
      NALU_ALIGNED ftype jact[dim][dim];
      for (int i=0; i<dim; ++i) 
        for (int j=0; j<dim; ++j) 
          jact[i][j] = ftype(0.0);

      NALU_ALIGNED ftype refGrad[AlgTraits::nodesPerElement_][dim];
      for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
        for (int i=0; i<dim; ++i) {
          refGrad[n][i] = referenceGradWeights(ip, n, i);
        }
        for (int i=0; i<dim; ++i) {
          for (int j=0; j<dim; ++j) {
            jact[i][j] += refGrad[n][j] * coords(n, i);
          }
        }
      }

      NALU_ALIGNED ftype adjJac[dim][dim];
      cofactorMatrix(adjJac, jact);

      NALU_ALIGNED ftype det = ftype(0.0);
      for (int i=0; i<dim; ++i) det += jact[i][0] * adjJac[i][0];
      STK_ThrowAssertMsg(
        stk::simd::are_any(det > tiny_positive_value()),
        "Problem with Jacobian determinant"
      );

      NALU_ALIGNED const ftype inv_detj = ftype(1.0) / det;

      for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
        for (int i=0; i<dim; ++i) {
          weights(ip, n, i) = ftype(0.0);
          for (int j=0; j<dim; ++j) {
            weights(ip, n, i) += adjJac[i][j] * refGrad[n][j];
          }
          weights(ip, n, i) *= inv_detj;
        }
      }
    }
  }

} // namespace nalu
} // namespace Sierra

#endif
