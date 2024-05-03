/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AlgTraits_h
#define AlgTraits_h

#include <stk_topology/topology.hpp>

namespace sierra {
namespace nalu {

// limited supported now (P=1 3D elements)
struct AlgTraitsHex8 {
  static constexpr int nDim_ = 3;
  static constexpr int nodesPerElement_ = 8;
  static constexpr int numScsIp_ = 12;
  static constexpr int numScvIp_ = 8;
  static constexpr int numGp_ = 8; // for FEM
  static constexpr stk::topology::topology_t topo_ = stk::topology::HEX_8;
};

struct AlgTraitsQuad4
{
  static constexpr int nDim_ = 3;
  static constexpr int nodesPerElement_ = 4;
  static constexpr int nodesPerFace_ = nodesPerElement_;
  static constexpr int numScsIp_ = 4;
  static constexpr int numFaceIp_ = numScsIp_;
  static constexpr stk::topology::topology_t topo_ = stk::topology::QUAD_4;
};

//-------------------------------------------------------------------------------------------
template <typename AlgTraitsFace, typename AlgTraitsElem>
struct AlgTraitsFaceElem
{
  using FaceTraits = AlgTraitsFace;
  using ElemTraits = AlgTraitsElem;

  static constexpr int nDim_ = ElemTraits::nDim_;
  static_assert( nDim_ == FaceTraits::nDim_, "inconsistent dimension specification");

  static constexpr int nodesPerElement_ = ElemTraits::nodesPerElement_;
  static constexpr int nodesPerFace_ = FaceTraits::nodesPerElement_;
  
  static constexpr int numScsIp_ = ElemTraits::numScsIp_;
  static constexpr int numScvIp_ = ElemTraits::numScvIp_;
  static constexpr int numGp_ = ElemTraits::numGp_;
  
  static constexpr int numFaceIp_ = FaceTraits::numFaceIp_;
  
  static constexpr stk::topology::topology_t elemTopo_ = ElemTraits::topo_;
  static constexpr stk::topology::topology_t faceTopo_ = FaceTraits::topo_;
};

using AlgTraitsQuad4Hex8 = AlgTraitsFaceElem<AlgTraitsQuad4, AlgTraitsHex8>;

} // namespace nalu
} // namespace Sierra

#endif
