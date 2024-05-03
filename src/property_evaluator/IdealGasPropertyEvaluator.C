/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <property_evaluator/PropertyEvaluator.h>
#include <property_evaluator/IdealGasPropertyEvaluator.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <vector>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// IdealGasPTYkrefPropertyEvaluator - evaluates density as a function of 
//                                    P, T and Ykref
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPTYkrefPropertyEvaluator::IdealGasPTYkrefPropertyEvaluator(
  double universalR,
  std::vector<std::pair<double, double> > mwMassFracVec,
  stk::mesh::MetaData &metaData)
  : PropertyEvaluator(),
    R_(universalR),
    mwRef_(0.0),
    pressure_(NULL)
{

  // save off mass fraction field
  pressure_ = metaData.get_field<double>(stk::topology::NODE_RANK, "pressure");

  // compute mixture mw
  double sum = 0.0;
  for (std::size_t k = 0; k < mwMassFracVec.size(); ++k ){
    sum += mwMassFracVec[k].second/mwMassFracVec[k].first;
  }
  mwRef_ = 1.0/sum;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
IdealGasPTYkrefPropertyEvaluator::~IdealGasPTYkrefPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
IdealGasPTYkrefPropertyEvaluator::execute(
  double *indVarList,
  stk::mesh::Entity node)
{
  const double T = indVarList[0];
  const double P = *stk::mesh::field_data(*pressure_, node);
  return P*mwRef_/R_/T;
}

} // namespace nalu
} // namespace Sierra
