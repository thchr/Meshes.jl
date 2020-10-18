module Meshes

using StaticArrays
using IterTools
using LinearAlgebra
using EarCut_jll

using Base: @propagate_inbounds

import Base: +, -

# basic concepts
include("vectors.jl")
include("matrices.jl")
include("points.jl")

include("geometries.jl")

include("interfaces.jl")
include("viewtypes.jl")
include("primitives.jl")
include("meshes.jl")
include("triangulation.jl")
include("lines.jl")
include("boundingboxes.jl")

# points
export Point, Point2, Point3, Point2f, Point3f

# vectors
export Vec, Vec2, Vec3, Vec2f, Vec3f
export vunit, vfill

# geometries
export Geometry, Primitive, Polytope
export Rectangle, Cylinder, Pyramid, Sphere
export coordtype

export LineFace, Line, NgonFace
export LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon
export Simplex, connect, Triangle, NSimplex, Tetrahedron
export coordinates, TetrahedronFace
export AbstractFace, TriangleFace, QuadFace
export FaceView, SimpleFaceView, TupleView
export decompose, coordinates, faces, convert_simplex
export Tesselation

# TODO: review these
export Cylinder2, Cylinder3

export height, origin, radius, width, widths, xwidth, yheight
export before, during, isinside, isoutside, overlaps, intersects, finishes
export direction, area, update, extremity
export self_intersections, split_intersections

# meshes
export AbstractMesh, Mesh, TriangleMesh, PlainMesh, PlainTriangleMesh

# mesh creation
export triangle_mesh

# bounding boxes
export boundingbox

end # module
