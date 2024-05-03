#include <gtest/gtest.h>
#include "NaluEnv.h"

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_topology/topology.hpp>

#include "nalu_make_unique.h"

#include "UnitTestUtils.h"

#include <algorithm>
#include <string>
#include <array>
#include <random>

namespace unit_test_utils {


stk::mesh::Entity create_one_element(
  stk::mesh::BulkData& bulk,
  stk::topology topo, const
  std::vector<std::vector<double>>& nodeLocations)
{
  // create just one element

   auto& meta = bulk.mesh_meta_data();
   stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", topo);
   stk::io::put_io_part_attribute(block_1);
   stk::mesh::PartVector allSurfaces = { &meta.declare_part("all_surfaces", meta.side_rank()) };
   stk::io::put_io_part_attribute(*allSurfaces.front());

   stk::mesh::PartVector individualSurfaces(topo.num_sides());
   for (unsigned k = 0u; k < topo.num_sides(); ++k) {
     individualSurfaces[k] = &meta.declare_part_with_topology("surface_" + std::to_string(k), topo.side_topology(k));
   }

   // set a coordinate field
   auto& coordField = meta.declare_field<double>(stk::topology::NODE_RANK, "coordinates");
   stk::mesh::put_field_on_mesh(coordField, block_1, meta.spatial_dimension(), nullptr);
   stk::mesh::put_field_on_mesh(coordField, stk::mesh::selectUnion(allSurfaces), meta.spatial_dimension(), nullptr);
   meta.set_coordinate_field(&coordField);
   meta.commit();

   stk::mesh::EntityIdVector nodeIds(topo.num_nodes());
   std::iota(nodeIds.begin(), nodeIds.end(), 1);

   
   bulk.modification_begin();

   for (auto id : nodeIds) {
     bulk.declare_entity(stk::topology::NODE_RANK, id, stk::mesh::PartVector{});
   }
   auto elem = stk::mesh::declare_element(bulk, block_1, bulk.parallel_rank()+1, nodeIds);

   bulk.modification_end();

   stk::mesh::create_all_sides(bulk, block_1, allSurfaces, false);

   auto surfaceSelector = stk::mesh::selectUnion(allSurfaces);
   stk::mesh::EntityVector all_faces;
   stk::mesh::get_selected_entities(surfaceSelector, bulk.get_buckets(meta.side_rank(), surfaceSelector), all_faces);
   STK_ThrowRequire(all_faces.size() == topo.num_sides());

   bulk.modification_begin();
   for (unsigned k = 0u; k < all_faces.size(); ++k) {
     const int ordinal = bulk.begin_element_ordinals(all_faces[k])[0];
     bulk.change_entity_parts(all_faces[k], {individualSurfaces[ordinal]}, stk::mesh::PartVector{});
   }
   bulk.modification_end();

   const auto* nodes = bulk.begin_nodes(elem);
   for (unsigned j = 0; j  < bulk.num_nodes(elem); ++j) {
     for (unsigned i = 0; i < topo.dimension(); ++i) {
       stk::mesh::field_data(coordField, nodes[j])[i] = nodeLocations.at(j).at(i);
     }
   }

   return elem;
}

stk::mesh::Entity create_one_reference_hex8_element(stk::mesh::BulkData& bulk)
{
   std::vector<std::vector<double>> nodeLocations =
   {
       {-0.5,-0.5,-0.5}, {+0.5,-0.5,-0.5}, {+0.5,+0.5,-0.5}, {-0.5,+0.5,-0.5},
       {-0.5,-0.5,+0.5}, {+0.5,-0.5,+0.5}, {+0.5,+0.5,+0.5}, {-0.5,+0.5,+0.5}
   };
   return create_one_element(bulk, stk::topology::HEX_8, nodeLocations);
}

stk::mesh::Entity create_one_reference_element(stk::mesh::BulkData& bulk, stk::topology topo)
{
  switch (topo.value())
  {
    case stk::topology::HEXAHEDRON_8:
      return create_one_reference_hex8_element(bulk);

    default:
      EXPECT_TRUE(false) << " Element type " + topo.name() + " not implemented";
      return stk::mesh::Entity(stk::mesh::Entity::InvalidEntity);
  }
}

}
