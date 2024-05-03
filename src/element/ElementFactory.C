/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "element/ElementFactory.h"
#include "element/Element.h"

// CVFEM-based
#include "element/Hex8CVFEM.h"
#include "element/Quad43DCVFEM.h"

#include "NaluEnv.h"
#include "nalu_make_unique.h"

#include <stk_util/util/ReportHandler.hpp>
#include <stk_topology/topology.hpp>

#include <cmath>
#include <iostream>
#include <memory>

namespace sierra{
namespace nalu{

  std::unique_ptr<Element>
  create_surface_element(stk::topology topo)
  {

    switch ( topo.value() ) {

      case stk::topology::HEX_8:
        return make_unique<HexSCS>();

      case stk::topology::QUAD_4:
        return make_unique<Quad3DSCS>();

      default:
        NaluEnv::self().naluOutputP0() << "sorry, we only support hex8,"
                                          " quad3d surface elements" << std::endl;
        NaluEnv::self().naluOutputP0() << "your type is " << topo.value() << std::endl;
        break;

    }
    return nullptr;
  }
  //--------------------------------------------------------------------------
  std::unique_ptr<Element>
  create_volume_element(stk::topology topo)
  {

    switch ( topo.value() ) {

      case stk::topology::HEX_8:
        return make_unique<HexSCV>();

      default:
        NaluEnv::self().naluOutputP0() << "sorry, we only support hex8 " << std::endl;
        NaluEnv::self().naluOutputP0() << "your type is " << topo.value() << std::endl;
        break;
    }
    return nullptr;
  }

  std::map<stk::topology, std::unique_ptr<Element>> ElementRepo::surfaceMeMap_;

  Element* ElementRepo::get_surface_element(
    const stk::topology& theTopo)
  {
    auto it = surfaceMeMap_.find(theTopo);
    if (it == surfaceMeMap_.end()) {
      surfaceMeMap_[theTopo] = create_surface_element(theTopo);
    }
    Element* theElem = surfaceMeMap_.at(theTopo).get();
    STK_ThrowRequire(theElem != nullptr);
    return theElem;
  }

  std::map<stk::topology, std::unique_ptr<Element>> ElementRepo::volumeMeMap_;

  Element* ElementRepo::get_volume_element(
    const stk::topology& theTopo)
  {
    auto it = volumeMeMap_.find(theTopo);
    if (it == volumeMeMap_.end()) {
      volumeMeMap_[theTopo] = create_volume_element(theTopo);
    }
    Element* theElem = volumeMeMap_.at(theTopo).get();
    STK_ThrowRequire(theElem != nullptr);
    return theElem;
  }

  void ElementRepo::clear()
  {
    surfaceMeMap_.clear();
    volumeMeMap_.clear();
  }

}
}
