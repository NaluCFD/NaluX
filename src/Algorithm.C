/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <Algorithm.h>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// Algorithm - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Algorithm::Algorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : realm_(realm)
{
  // push back on partVec
  partVec_.push_back(part);
}

// alternative; provide full partVec
Algorithm::Algorithm(
  Realm &realm,
  stk::mesh::PartVector &partVec)
  : realm_(realm),
    partVec_(partVec)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Algorithm::~Algorithm()
{
  // nothing to do
}

} // namespace nalu
} // namespace Sierra
