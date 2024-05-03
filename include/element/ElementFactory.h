/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ElementFactory_h
#define ElementFactory_h

#include <string>
#include <map>
#include <memory>

namespace stk { struct topology; }

namespace sierra{
namespace nalu{
  class Element;

  struct ElementRepo
  {
  public:
    static Element*
    get_surface_element(
      const stk::topology& theTopo);

    static Element*
    get_volume_element(
      const stk::topology& theTopo);

    static Element*
    get_fem_element(
      const stk::topology& theTopo);

    static void clear();
  private:
    ElementRepo() = default;
    static std::map<stk::topology, std::unique_ptr<Element>> surfaceMeMap_;
    static std::map<stk::topology, std::unique_ptr<Element>> volumeMeMap_;
  };

} // namespace nalu
} // namespace sierra

#endif
