/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <Realm.h>
#include <BoundaryConditions.h>
#include <NaluEnv.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>
#include <NaluParsing.h>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// BoundaryCondition - do some stuff
//==========================================================================
BoundaryCondition * BoundaryCondition::load(const YAML::Node & node) 
{
  if ( node["symmetry_boundary_condition"]) {
    SymmetryBoundaryConditionData& symmetryBC = *new SymmetryBoundaryConditionData(*parent());
    node >> symmetryBC;
    NaluEnv::self().naluOutputP0() << "Symmetry BC name:    " << symmetryBC.bcName_
                    << " on " << symmetryBC.targetName_ << std::endl;
    return &symmetryBC;
  }
  else if (node["periodic_boundary_condition"]) {
    PeriodicBoundaryConditionData& periodicBC = *new PeriodicBoundaryConditionData(*parent());
    node >> periodicBC;
    NaluEnv::self().naluOutputP0() << "Periodic BC name:    " << periodicBC.bcName_
                    << " between " << periodicBC.monarchSubject_.monarch_
                    << " and "<< periodicBC.monarchSubject_.subject_ << std::endl;
    return &periodicBC;
  }

  else {
    throw std::runtime_error("parser error BoundaryConditions::load: no such bc type");
  }
  return 0;
}

  Simulation* BoundaryCondition::root() { return parent()->root(); }
  BoundaryConditions *BoundaryCondition::parent() { return &boundaryConditions_; }

  Simulation* BoundaryConditions::root() { return parent()->root(); }
  Realm *BoundaryConditions::parent() { return &realm_; }

} // namespace nalu
} // namespace Sierra
