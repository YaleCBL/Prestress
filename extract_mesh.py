#!/usr/bin/env python3
"""
Extract mesh data from eva.vtu for the prestress algorithm.

Geometry: four-chamber patient-specific cardiac mesh from Pfaller et al. (2019,
2020).  282,288 nodes, 167,232 quadratic tetrahedral (tet10) elements, 846,864
structural DOFs.  Imaging data acquired at King's College London (Philips
Achieva 1.5T), meshed with Gmsh at 2 mm resolution.

Boundary conditions for prestress:
  LV endocardium (dsurf2)   — Neumann pressure
  RV endocardium (dsurf5)   — Neumann pressure
  LA endocardium (dsurf6)   — Neumann pressure
  RA endocardium (dsurf7)   — Neumann pressure
  External tissue (dsurf20) — homogeneous Dirichlet (zero displacement)
  All other surfaces        — homogeneous Neumann (traction-free)

Outputs (matching ventricle_preStress_main.m variables):
  VT             — node coordinates                      (n x 3)
  ET             — element connectivity                  (m x 10, tet10)
  F_pressure_lv  — LV endocardial surface faces          (f1 x 6, tri6)
  F_pressure_rv  — RV endocardial surface faces          (f2 x 6, tri6)
  F_pressure_la  — LA endocardial surface faces          (f3 x 6, tri6)
  F_pressure_ra  — RA endocardial surface faces          (f4 x 6, tri6)
  bcSupportList  — Dirichlet BC node indices             (b x 1)
  fiber_dir_IP   — fiber directions at Gauss points      (m x 4 x 3)
  element_mat    — element material IDs                  (m x 1)
"""

import argparse
import sys
from collections import Counter

import meshio
import numpy as np


# ---------------------------------------------------------------------------
# Mesh reading
# ---------------------------------------------------------------------------

def read_vtu(vtu_path):
    """Read VTU and return the meshio Mesh object."""
    print(f"Reading {vtu_path} ...")
    mesh = meshio.read(vtu_path)
    print(f"  Nodes:    {mesh.points.shape[0]}")
    for i, blk in enumerate(mesh.cells):
        print(f"  Cells[{i}]: {blk.type}, {blk.data.shape[0]} elements")
    return mesh


# ---------------------------------------------------------------------------
# Boundary face extraction (tet10)
# ---------------------------------------------------------------------------

# Local face definitions for a tet10 element.
#   tet10 node order: 0-3 corners, 4=mid(0-1), 5=mid(1-2), 6=mid(0-2),
#                     7=mid(0-3), 8=mid(1-3), 9=mid(2-3)
FACE_CORNER = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
FACE_FULL = [
    (0, 2, 1, 6, 5, 4),
    (0, 1, 3, 4, 8, 7),
    (1, 2, 3, 5, 9, 8),
    (0, 3, 2, 7, 9, 6),
]


def extract_boundary_faces(tet10):
    """Return dict mapping boundary face (sorted-corner-tuple) -> (elem_idx, local_face_idx)."""
    face_count = Counter()
    face_info = {}
    for ei in range(tet10.shape[0]):
        for fi, (a, b, c) in enumerate(FACE_CORNER):
            key = tuple(sorted((tet10[ei, a], tet10[ei, b], tet10[ei, c])))
            face_count[key] += 1
            face_info[key] = (ei, fi)
    return {k: face_info[k] for k, cnt in face_count.items() if cnt == 1}


def faces_on_surface(boundary, tet10, surf_nodes):
    """Return (f x 6) array of tri6 boundary faces whose nodes all lie in *surf_nodes*."""
    faces = []
    for fkey, (ei, fi) in boundary.items():
        local = FACE_FULL[fi]
        global_nodes = [tet10[ei, n] for n in local]
        if all(n in surf_nodes for n in global_nodes):
            faces.append(global_nodes)
    return np.array(faces, dtype=np.int64) if faces else np.empty((0, 6), dtype=np.int64)


# ---------------------------------------------------------------------------
# Surface / BC identification
# ---------------------------------------------------------------------------

def dsurf_nodes(mesh, name):
    """Return the set of node indices flagged by a dsurf array."""
    return set(np.where(mesh.point_data[name] != 0)[0])


def get_surface_faces(mesh, tet10, boundary, surf_name):
    """Extract boundary faces for a single named surface."""
    nodes = dsurf_nodes(mesh, surf_name)
    f = faces_on_surface(boundary, tet10, nodes)
    print(f"  {surf_name}: {f.shape[0]} tri6 faces ({len(nodes)} nodes)")
    return f


def get_bc_nodes(mesh, surf_name):
    """Return Dirichlet BC node indices from a named dsurf."""
    nodes = np.where(mesh.point_data[surf_name] != 0)[0]
    print(f"  {surf_name}: {len(nodes)} nodes (Dirichlet BC)")
    return nodes


# ---------------------------------------------------------------------------
# Fiber interpolation to integration points (tet10)
# ---------------------------------------------------------------------------

# 4-point Gauss quadrature for tetrahedra (exact for quadratic polynomials)
_A = (5.0 - np.sqrt(5.0)) / 20.0
_B = (5.0 + 3.0 * np.sqrt(5.0)) / 20.0
GAUSS_TET4 = np.array([
    [_A, _A, _A],
    [_B, _A, _A],
    [_A, _B, _A],
    [_A, _A, _B],
])
GAUSS_W4 = np.full(4, 1.0 / 24.0)


def tet10_shape(xi, eta, zeta):
    """Evaluate tet10 shape functions at natural coordinates (xi, eta, zeta).

    Barycentric coordinates: L1=1-xi-eta-zeta, L2=xi, L3=eta, L4=zeta.
    """
    L1 = 1.0 - xi - eta - zeta
    L2, L3, L4 = xi, eta, zeta
    return np.array([
        L1 * (2 * L1 - 1),
        L2 * (2 * L2 - 1),
        L3 * (2 * L3 - 1),
        L4 * (2 * L4 - 1),
        4 * L1 * L2,
        4 * L2 * L3,
        4 * L1 * L3,
        4 * L1 * L4,
        4 * L2 * L4,
        4 * L3 * L4,
    ])


def interpolate_to_gauss(tet10, node_data, gauss_pts=GAUSS_TET4):
    """Interpolate any nodal field to Gauss points (vectorized over elements)."""
    nElem = tet10.shape[0]
    nGP = gauss_pts.shape[0]
    d = node_data.shape[1] if node_data.ndim > 1 else 1
    if node_data.ndim == 1:
        node_data = node_data[:, None]

    N_all = np.array([tet10_shape(*gp) for gp in gauss_pts])  # (nGP, 10)

    ip_data = np.zeros((nElem, nGP, d))
    for g in range(nGP):
        N = N_all[g]
        for a in range(10):
            ip_data[:, g, :] += N[a] * node_data[tet10[:, a]]

    return ip_data


def compute_fiber_at_ip(mesh, tet10, fiber_key="node-fiber1"):
    """Interpolate nodal fiber vectors to Gauss points and normalise."""
    node_fibers = mesh.point_data[fiber_key]
    VT = mesh.points

    print(f"  Interpolating '{fiber_key}' to {GAUSS_TET4.shape[0]} Gauss points ...")
    fiber_ip = interpolate_to_gauss(tet10, node_fibers)
    pos_ip = interpolate_to_gauss(tet10, VT)

    norms = np.linalg.norm(fiber_ip, axis=2, keepdims=True)
    norms[norms < 1e-12] = 1.0
    fiber_ip = fiber_ip / norms

    return fiber_ip, pos_ip


# ---------------------------------------------------------------------------
# Default surface assignments for eva.vtu
# ---------------------------------------------------------------------------
ENDO_LV = "dsurf2"
ENDO_RV = "dsurf5"
ENDO_LA = "dsurf6"
ENDO_RA = "dsurf7"
EPICARDIUM = "dsurf1"
VALVE_MITRAL = "dsurf9"
VALVE_TRICUSPID = "dsurf10"
VALVE_AORTIC = "dsurf11"
VALVE_PULMONARY = "dsurf12"
DIRICHLET = "dsurf20"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("vtu", nargs="?", default="eva.vtu",
                        help="Input VTU file (default: eva.vtu)")
    parser.add_argument("--endo-lv", default=ENDO_LV,
                        help=f"LV endocardial surface (default: {ENDO_LV})")
    parser.add_argument("--endo-rv", default=ENDO_RV,
                        help=f"RV endocardial surface (default: {ENDO_RV})")
    parser.add_argument("--endo-la", default=ENDO_LA,
                        help=f"LA endocardial surface (default: {ENDO_LA})")
    parser.add_argument("--endo-ra", default=ENDO_RA,
                        help=f"RA endocardial surface (default: {ENDO_RA})")
    parser.add_argument("--epicardium", default=EPICARDIUM,
                        help=f"Epicardial surface (default: {EPICARDIUM})")
    parser.add_argument("--valve-mitral", default=VALVE_MITRAL,
                        help=f"Mitral valve ring surface (default: {VALVE_MITRAL})")
    parser.add_argument("--valve-tricuspid", default=VALVE_TRICUSPID,
                        help=f"Tricuspid valve ring surface (default: {VALVE_TRICUSPID})")
    parser.add_argument("--valve-aortic", default=VALVE_AORTIC,
                        help=f"Aortic valve ring surface (default: {VALVE_AORTIC})")
    parser.add_argument("--valve-pulmonary", default=VALVE_PULMONARY,
                        help=f"Pulmonary valve ring surface (default: {VALVE_PULMONARY})")
    parser.add_argument("--dirichlet", default=DIRICHLET,
                        help=f"Dirichlet BC surface (default: {DIRICHLET})")
    parser.add_argument("--fiber-key", default="node-fiber1",
                        help="Point data key for fiber direction (default: node-fiber1)")
    parser.add_argument("--fiber-key2", default="node-fiber2",
                        help="Point data key for second fiber family (default: node-fiber2)")
    parser.add_argument("--material", type=int, nargs="*", default=[1],
                        help="Restrict to specific volume material IDs (default: [1] = LV). "
                             "Pass --material with no values to keep all materials.")
    parser.add_argument("--output", default="mesh_data",
                        help="Output file base name (default: mesh_data)")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Read mesh
    # ------------------------------------------------------------------
    mesh = read_vtu(args.vtu)
    VT = mesh.points

    # ------------------------------------------------------------------
    # Extract tet10 elements
    # ------------------------------------------------------------------
    tet_block = None
    for blk in mesh.cells:
        if blk.type == "tetra10":
            tet_block = blk
            break
    if tet_block is None:
        sys.exit("ERROR: no tetra10 cells found in VTU")

    ET = tet_block.data
    tet_mat = mesh.cell_data["element-material"][0].flatten()

    if args.material:
        mask = np.isin(tet_mat, args.material)
        ET = ET[mask]
        tet_mat = tet_mat[mask]
        print(f"\nFiltering to materials {args.material}: {ET.shape[0]} elements")

        # Subset nodes to those used by the filtered elements and remap
        # connectivity so all downstream extraction (boundary faces, dsurf
        # surfaces, BC nodes, fibers) is restricted to this geometry.
        used_nodes = np.unique(ET)
        node_map = -np.ones(VT.shape[0], dtype=np.int64)
        node_map[used_nodes] = np.arange(used_nodes.size)
        ET = node_map[ET]
        VT = VT[used_nodes]
        mesh.points = VT
        mesh.point_data = {k: v[used_nodes] for k, v in mesh.point_data.items()}
        print(f"  Retained {VT.shape[0]} of {node_map.size} nodes")

    print(f"\nNodes (VT): {VT.shape}")
    print(f"Elements (ET): {ET.shape}")

    # ------------------------------------------------------------------
    # Boundary faces
    # ------------------------------------------------------------------
    print("\nExtracting boundary faces ...")
    boundary = extract_boundary_faces(ET)
    print(f"  Total boundary faces: {len(boundary)}")

    # ------------------------------------------------------------------
    # Pressure surfaces (one per chamber). Empty surfaces are dropped.
    # ------------------------------------------------------------------
    print("\nIdentifying endocardial (pressure) surfaces ...")
    pressure_specs = [
        ("F_pressure_lv", args.endo_lv, "LV"),
        ("F_pressure_rv", args.endo_rv, "RV"),
        ("F_pressure_la", args.endo_la, "LA"),
        ("F_pressure_ra", args.endo_ra, "RA"),
    ]
    pressure_faces = {}
    for var_name, surf_name, label in pressure_specs:
        f = get_surface_faces(mesh, ET, boundary, surf_name)
        if f.shape[0] > 0:
            pressure_faces[var_name] = (f, surf_name, label)
        else:
            print(f"    -> dropped (no faces in this geometry)")

    # ------------------------------------------------------------------
    # Other named surfaces: epicardium and valve rings. Dropped if empty.
    # ------------------------------------------------------------------
    print("\nIdentifying epicardium and valve-ring surfaces ...")
    extra_specs = [
        ("F_epicardium",      args.epicardium,      "epicardium"),
        ("F_valve_mitral",    args.valve_mitral,    "mitral valve ring"),
        ("F_valve_tricuspid", args.valve_tricuspid, "tricuspid valve ring"),
        ("F_valve_aortic",    args.valve_aortic,    "aortic valve ring"),
        ("F_valve_pulmonary", args.valve_pulmonary, "pulmonary valve ring"),
    ]

    extra_faces = {}
    for var_name, surf_name, label in extra_specs:
        f = get_surface_faces(mesh, ET, boundary, surf_name)
        if f.shape[0] > 0:
            extra_faces[var_name] = (f, surf_name, label)
        else:
            print(f"    -> dropped (no faces in this geometry)")

    # ------------------------------------------------------------------
    # Dirichlet BC nodes (external tissue). Dropped if empty.
    # ------------------------------------------------------------------
    print("\nIdentifying Dirichlet BC nodes ...")
    bcSupportList = get_bc_nodes(mesh, args.dirichlet)
    if bcSupportList.size == 0:
        print(f"    -> dropped (no Dirichlet nodes in this geometry)")
        bcSupportList = None

    # ------------------------------------------------------------------
    # Fiber directions at integration points
    # ------------------------------------------------------------------
    print("\nComputing fiber directions at integration points ...")
    fiber_dir_IP, pos_IP = compute_fiber_at_ip(mesh, ET, args.fiber_key)
    print(f"  fiber_dir_IP: {fiber_dir_IP.shape}")

    fiber_dir_IP2 = None
    if args.fiber_key2 in mesh.point_data:
        print(f"\nComputing second fiber family ...")
        fiber_dir_IP2, _ = compute_fiber_at_ip(mesh, ET, args.fiber_key2)
        print(f"  fiber_dir_IP2: {fiber_dir_IP2.shape}")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("SUMMARY")
    print("=" * 65)
    print(f"  VT             : {VT.shape}           (node coordinates)")
    print(f"  ET             : {ET.shape}      (tet10 connectivity)")
    for var_name, (f, surf_name, label) in pressure_faces.items():
        print(f"  {var_name:18s} : {f.shape}        ({label} endocardial faces, {surf_name})")
    for var_name, (f, surf_name, label) in extra_faces.items():
        print(f"  {var_name:18s} : {f.shape}        ({label} faces, {surf_name})")
    if bcSupportList is not None:
        print(f"  bcSupportList  : {bcSupportList.shape}        (Dirichlet BC nodes, {args.dirichlet})")
    print(f"  fiber_dir_IP   : {fiber_dir_IP.shape}  (fiber1 at 4 GP)")
    if fiber_dir_IP2 is not None:
        print(f"  fiber_dir_IP2  : {fiber_dir_IP2.shape}  (fiber2 at 4 GP)")
    print(f"  pos_IP         : {pos_IP.shape}  (Gauss point positions)")
    print(f"  element_mat    : {tet_mat.shape}        (material IDs)")

    # ------------------------------------------------------------------
    # Save .npz (0-indexed)
    # ------------------------------------------------------------------
    npz_path = args.output + ".npz"
    save_dict = dict(
        VT=VT,
        ET=ET,
        fiber_dir_IP=fiber_dir_IP,
        pos_IP=pos_IP,
        element_material=tet_mat,
        gauss_points=GAUSS_TET4,
        gauss_weights=GAUSS_W4,
    )
    for var_name, (f, _, _) in pressure_faces.items():
        save_dict[var_name] = f
    for var_name, (f, _, _) in extra_faces.items():
        save_dict[var_name] = f
    if bcSupportList is not None:
        save_dict["bcSupportList"] = bcSupportList
    if fiber_dir_IP2 is not None:
        save_dict["fiber_dir_IP2"] = fiber_dir_IP2
    np.savez(npz_path, **save_dict)
    print(f"\nSaved (0-indexed): {npz_path}")

    # ------------------------------------------------------------------
    # Save .mat (1-indexed for MATLAB)
    # ------------------------------------------------------------------
    try:
        from scipy.io import savemat
        mat_dict = dict(
            VT=VT,
            ET=ET + 1,
            fiber_dir_IP=fiber_dir_IP,
            pos_IP=pos_IP,
            element_material=tet_mat,
            gauss_points=GAUSS_TET4,
            gauss_weights=GAUSS_W4,
        )
        for var_name, (f, _, _) in pressure_faces.items():
            mat_dict[var_name] = f + 1
        for var_name, (f, _, _) in extra_faces.items():
            mat_dict[var_name] = f + 1
        if bcSupportList is not None:
            mat_dict["bcSupportList"] = bcSupportList + 1
        if fiber_dir_IP2 is not None:
            mat_dict["fiber_dir_IP2"] = fiber_dir_IP2
        mat_path = args.output + ".mat"
        savemat(mat_path, mat_dict)
        print(f"Saved (1-indexed): {mat_path}")
    except ImportError:
        print("scipy not available - .mat file not saved")

    print(f"""
NOTES FOR COLLABORATOR:
  Geometry: four-chamber heart from Pfaller et al. (2019, 2020).

  Element type: tet10 (10-node quadratic tetrahedra), not hex8.
    -> FEBio: Elements type='tet10', 10 nodes per element.
  Surface faces: tri6 (6-node quadratic triangles), not quad4.
    -> FEBio: Surface type='tri6', 6 nodes per face.

  Boundary conditions (prestress):
    LV endocardium ({args.endo_lv}): Neumann pressure (one value)
    RV endocardium ({args.endo_rv}): Neumann pressure (one value)
    LA endocardium ({args.endo_la}): Neumann pressure (one value)
    RA endocardium ({args.endo_ra}): Neumann pressure (one value)
    External tissue ({args.dirichlet}): homogeneous Dirichlet (zero displacement)
    All other surfaces: homogeneous Neumann (traction-free)

  Fiber directions:
    node-fiber1 = myofiber direction (f_0)
    node-fiber2 = sheet direction (s_0)
    Interpolated to 4 Gauss points per tet10 element.

  The .mat file uses 1-indexed connectivity (MATLAB convention).
  The .npz file uses 0-indexed connectivity (Python convention).
""")


if __name__ == "__main__":
    main()
