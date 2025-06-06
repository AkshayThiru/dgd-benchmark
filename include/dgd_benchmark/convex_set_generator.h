#ifndef DGD_BENCHMARK_CONVEX_SET_GENERATOR_H_
#define DGD_BENCHMARK_CONVEX_SET_GENERATOR_H_

#include <memory>
#include <string>
#include <vector>

#include "dgd_benchmark/data_types.h"

/**
 * Mesh functions.
 */

// Mesh properties.
struct MeshProperties {
  std::vector<Vec3> vert, normal;
  std::vector<double> offset;
  std::vector<int> vgraph, fgraph;
  Vec3 interior_point;
  double inradius = 0.0;
  int nvert = 0;
  int nfacet = 0;
};

// Sets vertex mesh from .obj file.
void SetMeshFromObjFile(const std::string& filename, MeshProperties& mp);

// Sets the center of the vertex mesh as the origin.
inline void SetZeroVertexCenter(MeshProperties& mp) {
  for (auto& v : mp.vert) {
    v -= mp.interior_point;
  }
}

// Sets the center of the facet mesh as the origin.
inline void SetZeroFacetCenter(MeshProperties& mp) {
  for (int i = 0; i < mp.nfacet; ++i) {
    mp.offset[i] += mp.normal[i].dot(mp.interior_point);
  }
}

/**
 * Random convex set generator.
 */

// Forward declarations.
namespace dgd {

template <int dim>
class ConvexSet;

class Mesh;

}  // namespace dgd

namespace dsf {
template <int exp>
class VDSF;
}

namespace inc {
class Polyhedron;
}

// Feature ranges for convex sets.
template <typename T>
struct Range {
  T low{};
  T high{};
};

struct ConvexSetFeatureRange {
  struct {
    Range<double> hlx{{0.25 * 1e-2}, {0.25}};
    Range<double> radius{{0.25 * 1e-2}, {0.25}};
  } capsule;

  struct {
    Range<double> radius{{0.25 * 1e-2}, {0.25}};
    Range<double> height{{0.5 * 1e-2}, {0.5}};
  } cone;

  struct {
    Range<double> hlx{{0.4 * 1e-2}, {0.4}};
    Range<double> radius{{0.25 * 1e-2}, {0.25}};
  } cylinder;

  struct {
    Range<double> hlx{{0.25 * 1e-2}, {0.25}};
    Range<double> hly{{0.25 * 1e-2}, {0.25}};
    Range<double> hlz{{0.25 * 1e-2}, {0.25}};
  } ellipsoid;

  struct {
    Range<double> base_radius{{0.25 * 1e-2}, {0.25}};
    Range<double> top_radius{{0.25 * 1e-2}, {0.25}};
    Range<double> height{{0.5 * 1e-2}, {0.5}};
  } frustum;

  struct {
    Range<double> radius{{0.25 * 1e-2}, {0.25}};
  } sphere;

  struct {
    Range<double> hlx{{0.25 * 1e-2}, {0.25}};
    Range<double> hly{{0.25 * 1e-2}, {0.25}};
    Range<double> hlz{{0.25 * 1e-2}, {0.25}};
  } cuboid;

  struct {
    Range<int> nvert{{4}, {32}};
    double size = 0.4;
  } polytope;
};

// Convex set generator.
class ConvexSetGenerator {
 public:
  using ConvexSetPtr = std::shared_ptr<dgd::ConvexSet<3>>;
  using MeshPtr = std::shared_ptr<dgd::Mesh>;
  using PolyhedronPtr = std::shared_ptr<inc::Polyhedron>;
  using VdsfPtr = std::shared_ptr<dsf::VDSF<kDsfExp>>;

  enum class CurvedPrimitive {
    Capsule,
    Cone,
    Cylinder,
    Ellipsoid,
    Frustum,
    Sphere,
    Count_,
  };

  enum class FlatPrimitive {
    Cuboid,
    Polytope,
    Count_,
  };

  explicit ConvexSetGenerator(const ConvexSetFeatureRange& fr);

  ~ConvexSetGenerator() = default;

  // Load meshes from .obj files.
  void LoadMeshesFromObjFiles(const std::vector<std::string>& obj_filenames);

  // Generate a random convex set of a primitive type.
  ConvexSetPtr GeneratePrimitiveSet(CurvedPrimitive type) const;

  ConvexSetPtr GeneratePrimitiveSet(FlatPrimitive type);

  // Get random primitive convex sets.
  ConvexSetPtr GetRandomCurvedPrimitiveSet() const;

  ConvexSetPtr GetRandomPrimitiveSet();

  // Retrieve a random Mesh set with index.
  ConvexSetPtr GetRandomMeshSet(int* idx = nullptr) const;

  // Get a random convex set.
  ConvexSetPtr GetRandomSet();

  const std::vector<MeshPtr>& meshes() const { return meshes_; }

  const std::vector<PolyhedronPtr>& polyhedra() const { return polyhedra_; }

  const std::vector<VdsfPtr>& vdsfs() const { return vdsfs_; }

 private:
  std::vector<MeshPtr> meshes_;
  std::vector<PolyhedronPtr> polyhedra_;
  std::vector<VdsfPtr> vdsfs_;
  std::vector<Vec3> polytope_vert_;
  const ConvexSetFeatureRange fr_;
  const int ccount_, fcount_;
  int nmeshes_;
};

#endif  // DGD_BENCHMARK_CONVEX_SET_GENERATOR_H_
