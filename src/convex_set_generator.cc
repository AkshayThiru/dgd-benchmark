#include "dgd_benchmark/convex_set_generator.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "dgd/geometry/3d/cone.h"
#include "dgd/geometry/3d/cuboid.h"
#include "dgd/geometry/3d/cylinder.h"
#include "dgd/geometry/3d/ellipsoid.h"
#include "dgd/geometry/3d/frustum.h"
#include "dgd/geometry/3d/mesh.h"
#include "dgd/geometry/3d/polytope.h"
#include "dgd/geometry/xd/capsule.h"
#include "dgd/geometry/xd/sphere.h"
#include "dgd/mesh_loader.h"
#include "dgd/utils.h"
#include "dgd_benchmark/data_types.h"
#include "dgd_benchmark/dsf/dsf.h"
#include "dgd_benchmark/inc/polyhedron.h"

/**
 * Mesh functions.
 */

void SetMeshFromObjFile(const std::string& filename, MeshProperties& mp) {
  dgd::MeshLoader ml{};

  try {
    ml.LoadObj(filename);
    if (!ml.MakeVertexGraph(mp.vert, mp.vgraph) ||
        !ml.MakeFacetGraph(mp.normal, mp.offset, mp.fgraph,
                           mp.interior_point)) {
      throw std::runtime_error("Qhull error: Failed to parse the file");
    }
    if ((mp.inradius = ml.ComputeInradius(mp.interior_point, true)) <= 0.0) {
      throw std::runtime_error("Nonpositive inradius");
    }
  } catch (const std::runtime_error& e) {
    std::cerr << "Error loading mesh from file '" << filename
              << "': " << e.what() << std::endl;
    throw;  // Re-throw the exception to be handled by the caller.
  }

  mp.nvert = static_cast<int>(mp.vert.size());
  mp.nfacet = static_cast<int>(mp.normal.size());
}

/**
 * Random convex set generator.
 */

namespace {

inline double Random(Range<double> range) {
  return dgd::Random(range.low, range.high);
}

inline int Random(Range<int> range) {
  std::uniform_int_distribution<int> dis(range.low, range.high);
  return dis(dgd::generator);
}

}  // namespace

ConvexSetGenerator::ConvexSetGenerator(const ConvexSetFeatureRange& fr)
    : fr_(fr),
      ccount_(static_cast<int>(CurvedPrimitive::Count_)),
      fcount_(static_cast<int>(FlatPrimitive::Count_)) {
  dgd::SetDefaultSeed();
  if (fr.polytope.nvert.low < 4) {
    throw std::invalid_argument("Polytopes must have at least 4 vertices.");
  }
  polytope_vert_.reserve(fr.polytope.nvert.high);
  meshes_.clear();
  polyhedra_.clear();
  vdsfs_.clear();
  nmeshes_ = 0;
}

void ConvexSetGenerator::LoadMeshesFromObjFiles(
    const std::vector<std::string>& obj_filenames) {
  meshes_.clear();
  polyhedra_.clear();
  vdsfs_.clear();

  MeshProperties mp;
  for (const auto& filename : obj_filenames) {
    try {
      SetMeshFromObjFile(filename, mp);
    } catch (const std::runtime_error& e) {
      continue;
    }
    SetZeroVertexCenter(mp);
    SetZeroFacetCenter(mp);

    meshes_.push_back(
        std::make_shared<dgd::Mesh>(mp.vert, mp.vgraph, 0.0, mp.inradius));
    polyhedra_.push_back(
        std::make_shared<inc::Polyhedron>(mp.normal, mp.offset, mp.fgraph));
    vdsfs_.push_back(std::make_shared<dsf::VDSF<kDsfExp>>(mp.vert));
  }
  nmeshes_ = static_cast<int>(meshes_.size());
}

ConvexSetGenerator::ConvexSetPtr ConvexSetGenerator::GeneratePrimitiveSet(
    CurvedPrimitive type) const {
  switch (type) {
    case CurvedPrimitive::Capsule:
      return std::make_shared<dgd::Capsule>(Random(fr_.capsule.hlx),
                                            Random(fr_.capsule.radius), 0.0);

    case CurvedPrimitive::Cone:
      return std::make_shared<dgd::Cone>(Random(fr_.cone.radius),
                                         Random(fr_.cone.height), 0.0);

    case CurvedPrimitive::Cylinder:
      return std::make_shared<dgd::Cylinder>(Random(fr_.cylinder.hlx),
                                             Random(fr_.cylinder.radius), 0.0);

    case CurvedPrimitive::Ellipsoid:
      return std::make_shared<dgd::Ellipsoid>(Random(fr_.ellipsoid.hlx),
                                              Random(fr_.ellipsoid.hly),
                                              Random(fr_.ellipsoid.hlz), 0.0);

    case CurvedPrimitive::Frustum:
      return std::make_shared<dgd::Frustum>(Random(fr_.frustum.base_radius),
                                            Random(fr_.frustum.top_radius),
                                            Random(fr_.frustum.height), 0.0);

    case CurvedPrimitive::Sphere:
      return std::make_shared<dgd::Sphere>(Random(fr_.sphere.radius));

    default:
      throw std::invalid_argument("Invalid curved primitive type");
  }
}

ConvexSetGenerator::ConvexSetPtr ConvexSetGenerator::GeneratePrimitiveSet(
    FlatPrimitive type) {
  switch (type) {
    case FlatPrimitive::Cuboid:
      return std::make_shared<dgd::Cuboid>(Random(fr_.cuboid.hlx),
                                           Random(fr_.cuboid.hly),
                                           Random(fr_.cuboid.hlz), 0.0);

    case FlatPrimitive::Polytope:
      polytope_vert_.clear();

      if (meshes_.empty()) {
        for (int i = 0; i < Random(fr_.polytope.nvert); ++i) {
          Vec3 v{dgd::Random(), dgd::Random(), dgd::Random()};
          v /= (v.lpNorm<4>() + dgd::kEps);  //
          polytope_vert_.push_back(fr_.polytope.size * v);
        }
      } else {
        const int mesh_idx = Random(Range<int>{0, nmeshes_ - 1});
        const int nvert_m = meshes_[mesh_idx]->nvertices();
        for (int i = 0; i < 4; ++i) {
          polytope_vert_.push_back(
              meshes_[mesh_idx]->vertices()[(i * nvert_m) / 3]);
        }
        for (int i = 4; i < fr_.polytope.nvert.high; ++i) {
          const int idx = Random(Range<int>{0, nvert_m - 1});
          polytope_vert_.push_back(meshes_[mesh_idx]->vertices()[idx]);
        }
      }

      {
        Vec3 center;
        for (const auto& v : polytope_vert_) {
          center += v;
        }
        center /= static_cast<double>(polytope_vert_.size());
        for (auto& v : polytope_vert_) {
          v -= center;
        }
      }
      return std::make_shared<dgd::Polytope>(polytope_vert_, 0.0, dgd::kEps);

    default:
      throw std::invalid_argument("Invalid flat primitive type");
  }
}

ConvexSetGenerator::ConvexSetPtr
ConvexSetGenerator::GetRandomCurvedPrimitiveSet() const {
  const int idx = Random(Range<int>{0, ccount_ - 1});
  return GeneratePrimitiveSet(static_cast<CurvedPrimitive>(idx));
}

ConvexSetGenerator::ConvexSetPtr ConvexSetGenerator::GetRandomPrimitiveSet() {
  const int idx = Random(Range<int>{0, ccount_ + fcount_ - 1});
  if (idx < ccount_) {
    return GeneratePrimitiveSet(static_cast<CurvedPrimitive>(idx));
  } else {
    return GeneratePrimitiveSet(static_cast<FlatPrimitive>(idx - ccount_));
  }
}

ConvexSetGenerator::ConvexSetPtr ConvexSetGenerator::GetRandomMeshSet(
    int* idx) const {
  if (meshes_.empty()) {
    throw std::runtime_error(
        "No meshes available to retrieve a random mesh set");
  }

  const int mesh_idx = Random(Range<int>{0, nmeshes_ - 1});
  if (idx) *idx = mesh_idx;
  return meshes_[mesh_idx];
}

ConvexSetGenerator::ConvexSetPtr ConvexSetGenerator::GetRandomSet() {
  const int idx = Random(Range<int>{0, ccount_ + fcount_ + nmeshes_ - 1});
  if (idx < ccount_) {
    return GetRandomCurvedPrimitiveSet();
  } else if (idx < ccount_ + fcount_) {
    return GetRandomPrimitiveSet();
  } else {
    return GetRandomMeshSet();
  }
}
