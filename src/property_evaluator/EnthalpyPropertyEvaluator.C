/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <property_evaluator/EnthalpyPropertyEvaluator.h>

#include <FieldTypeDef.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <map>
#include <stdexcept>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// EnthalpyConstSpecHeatPropertyEvaluator - evaluates H based on Cp and Tref
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyConstSpecHeatPropertyEvaluator::EnthalpyConstSpecHeatPropertyEvaluator(
  const double & specificHeat,
  const double & referenceTemperature)
  : specificHeat_(specificHeat),
    referenceTemperature_(referenceTemperature)
{
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyConstSpecHeatPropertyEvaluator::~EnthalpyConstSpecHeatPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
EnthalpyConstSpecHeatPropertyEvaluator::execute(
    double *indVarList,
    stk::mesh::Entity /*node*/)
{
  const double T = indVarList[0];
  return specificHeat_ * (T - referenceTemperature_);
}

} // namespace nalu
} // namespace Sierra
