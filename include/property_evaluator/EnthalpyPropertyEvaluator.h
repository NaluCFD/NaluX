/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EnthalpyPropertyEvaluator_h
#define EnthalpyPropertyEvaluator_h

#include <property_evaluator/PropertyEvaluator.h>

#include <string>
#include <map>
#include <vector>

namespace stk {
namespace mesh {
  struct Entity;
}
}

namespace sierra{
namespace nalu{

class EnthalpyConstSpecHeatPropertyEvaluator : public PropertyEvaluator
{
public:

  EnthalpyConstSpecHeatPropertyEvaluator(
    const double & specificHeat,
    const double & referenceTemperature);

  virtual ~EnthalpyConstSpecHeatPropertyEvaluator();
  
  double execute(
    double *indVarList,
    stk::mesh::Entity node);

  double specificHeat_;
  double referenceTemperature_;
};

} // namespace nalu
} // namespace Sierra

#endif
