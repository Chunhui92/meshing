import numpy as np

from src.mir_mesh import UcdSiloMesh
from src.mir_utils import (
    CLAMP_EPS,
    EPS,
    HEX_FACES,
    VTK_HEXAHEDRON,
    VTK_TETRA,
    choose_split_normal,
    clip_tets_volume,
    hex_volume,
    split_rectilinear_hex,
    tet_volume,
)


class InMemoryLegacyVtkMesh:
    def __init__(self, original_points):
        original_points = np.asarray(original_points, dtype=np.float64)
        self.points = [original_points[idx].copy() for idx in range(len(original_points))]
        self.cells = []
        self.cell_types = []
        self.material_ids = []
        self.zone_ids = []
        self.material_volumes = {}
        self.edge_point_cache = {}
        self.center_point_ids = {}
        self.degenerate_cells = 0

    @property
    def num_cells(self):
        return len(self.cells)

    @property
    def num_points(self):
        return len(self.points)

    def add_point(self, point):
        self.points.append(np.asarray(point, dtype=np.float64))
        return len(self.points) - 1

    def get_or_add_center(self, zone_id, center):
        if zone_id not in self.center_point_ids:
            self.center_point_ids[zone_id] = self.add_point(center)
        return self.center_point_ids[zone_id]

    def get_or_add_edge_point(self, v0, v1, mat_pair, phi0, phi1):
        key = (min(v0, v1), max(v0, v1), mat_pair)
        point_id = self.edge_point_cache.get(key)
        if point_id is not None:
            return point_id

        denom = phi0 - phi1
        t = 0.5 if abs(denom) <= EPS else phi0 / denom
        t = min(max(t, CLAMP_EPS), 1.0 - CLAMP_EPS)
        point = self.points[v0] + t * (self.points[v1] - self.points[v0])
        point_id = self.add_point(point)
        self.edge_point_cache[key] = point_id
        return point_id

    def add_cell_ids(self, point_ids, cell_type, material_id, source_zone_id):
        if len(set(point_ids)) != len(point_ids):
            self.degenerate_cells += 1
            return False
        coords = np.asarray([self.points[idx] for idx in point_ids], dtype=np.float64)
        volume = tet_volume(coords) if cell_type == VTK_TETRA else hex_volume(coords)
        if volume <= 1.0e-12:
            self.degenerate_cells += 1
            return False
        self.cells.append(list(point_ids))
        self.cell_types.append(int(cell_type))
        self.material_ids.append(int(material_id))
        self.zone_ids.append(int(source_zone_id))
        self.material_volumes[material_id] = self.material_volumes.get(material_id, 0.0) + volume
        return True

    def add_cell_points(self, points, cell_type, material_id, source_zone_id):
        point_ids = [self.add_point(point) for point in np.asarray(points, dtype=np.float64)]
        return self.add_cell_ids(point_ids, cell_type, material_id, source_zone_id)

    def write(self, output_path):
        with open(output_path, "w", encoding="utf-8") as out:
            out.write("# vtk DataFile Version 3.0\n")
            out.write("multi-material MIR reconstruction\n")
            out.write("ASCII\n")
            out.write("DATASET UNSTRUCTURED_GRID\n")
            out.write("POINTS {} float\n".format(self.num_points))
            for point in self.points:
                out.write("{:.9g} {:.9g} {:.9g}\n".format(point[0], point[1], point[2]))

            connectivity_size = sum(len(cell) + 1 for cell in self.cells)
            out.write("CELLS {} {}\n".format(self.num_cells, connectivity_size))
            for cell in self.cells:
                out.write("{} {}\n".format(len(cell), " ".join(str(idx) for idx in cell)))

            out.write("CELL_TYPES {}\n".format(self.num_cells))
            for cell_type in self.cell_types:
                out.write("{}\n".format(cell_type))

            out.write("CELL_DATA {}\n".format(self.num_cells))
            out.write("SCALARS material_id int 1\n")
            out.write("LOOKUP_TABLE default\n")
            for material_id in self.material_ids:
                out.write("{}\n".format(material_id))
            out.write("SCALARS source_zone_id int 1\n")
            out.write("LOOKUP_TABLE default\n")
            for zone_id in self.zone_ids:
                out.write("{}\n".format(zone_id))


def iter_selected_zones(mesh, max_zones=None):
    for index, zone in enumerate(mesh.iter_zones()):
        if max_zones is not None and index >= max_zones:
            break
        yield zone


def compute_reference_material_volumes(mesh, max_zones=None):
    volumes = {}
    for zone in iter_selected_zones(mesh, max_zones=max_zones):
        for material_id, vf in zone.materials:
            volumes[material_id] = volumes.get(material_id, 0.0) + zone.volume * vf
    return volumes


def compute_material_centroids(mesh, max_zones=None):
    weighted_sum = {}
    weights = {}
    for zone in iter_selected_zones(mesh, max_zones=max_zones):
        for material_id, vf in zone.materials:
            weight = zone.volume * vf
            weighted_sum[material_id] = weighted_sum.get(material_id, np.zeros(3, dtype=np.float64)) + zone.center * weight
            weights[material_id] = weights.get(material_id, 0.0) + weight
    return {
        material_id: weighted_sum[material_id] / weights[material_id]
        for material_id in weights
        if weights[material_id] > 0.0
    }


def compute_node_vfs(mesh, max_zones=None):
    sums = [dict() for _ in range(mesh.num_nodes)]
    adjacency = np.zeros(mesh.num_nodes, dtype=np.int32)
    for zone in iter_selected_zones(mesh, max_zones=max_zones):
        for node_id in zone.node_ids:
            adjacency[node_id] += 1
            node_sum = sums[node_id]
            for material_id, vf in zone.materials:
                node_sum[material_id] = node_sum.get(material_id, 0.0) + vf

    node_vfs = []
    for node_id in range(mesh.num_nodes):
        if adjacency[node_id] == 0:
            node_vfs.append({})
            continue
        scale = 1.0 / float(adjacency[node_id])
        values = {
            material_id: value * scale
            for material_id, value in sums[node_id].items()
            if value * scale > EPS
        }
        total = sum(values.values())
        if total > 0.0:
            inv_total = 1.0 / total
            values = {material_id: value * inv_total for material_id, value in values.items()}
        node_vfs.append(values)
    return node_vfs


def normalize_materials(materials, min_fraction):
    active = [(int(material_id), float(vf)) for material_id, vf in materials if vf > min_fraction]
    if not active:
        return []
    total = sum(vf for _, vf in active)
    return [] if total <= 0.0 else [(material_id, vf / total) for material_id, vf in active]


def canonical_material_pair(material_a, material_b):
    return (int(min(material_a, material_b)), int(max(material_a, material_b)))


def build_consistent_hex_tets(zone_node_ids, center_id):
    tets = []
    for face in HEX_FACES:
        face_ids = [int(zone_node_ids[idx]) for idx in face]
        diag_02 = (min(face_ids[0], face_ids[2]), max(face_ids[0], face_ids[2]))
        diag_13 = (min(face_ids[1], face_ids[3]), max(face_ids[1], face_ids[3]))
        tris = ((face[0], face[1], face[2]), (face[0], face[2], face[3])) if diag_02 <= diag_13 else (
            (face[1], face[2], face[3]),
            (face[1], face[3], face[0]),
        )
        for tri in tris:
            tets.append((center_id, zone_node_ids[tri[0]], zone_node_ids[tri[1]], zone_node_ids[tri[2]]))
    return tets


def _edge_point(builder, tet_ids, phi_values, mat_pair, i, j):
    return builder.get_or_add_edge_point(
        int(tet_ids[i]),
        int(tet_ids[j]),
        mat_pair,
        float(phi_values[i]),
        float(phi_values[j]),
    )


def _negative_tets_for_scalar(builder, tet_ids, phi_values, mat_pair):
    inside = [idx for idx, value in enumerate(phi_values) if value <= EPS]
    if not inside:
        return []
    if len(inside) == 4:
        return [list(map(int, tet_ids))]

    outside = [idx for idx in range(4) if idx not in inside]
    if len(inside) == 1:
        a = inside[0]
        b, c, d = outside
        return [[
            int(tet_ids[a]),
            _edge_point(builder, tet_ids, phi_values, mat_pair, a, b),
            _edge_point(builder, tet_ids, phi_values, mat_pair, a, c),
            _edge_point(builder, tet_ids, phi_values, mat_pair, a, d),
        ]]

    if len(inside) == 2:
        a, b = inside
        c, d = outside
        pac = _edge_point(builder, tet_ids, phi_values, mat_pair, a, c)
        pad = _edge_point(builder, tet_ids, phi_values, mat_pair, a, d)
        pbc = _edge_point(builder, tet_ids, phi_values, mat_pair, b, c)
        pbd = _edge_point(builder, tet_ids, phi_values, mat_pair, b, d)
        return [
            [int(tet_ids[a]), pac, pad, int(tet_ids[b])],
            [pac, pad, int(tet_ids[b]), pbc],
            [pad, int(tet_ids[b]), pbc, pbd],
        ]

    a, b, c = inside
    d = outside[0]
    pad = _edge_point(builder, tet_ids, phi_values, mat_pair, a, d)
    pbd = _edge_point(builder, tet_ids, phi_values, mat_pair, b, d)
    pcd = _edge_point(builder, tet_ids, phi_values, mat_pair, c, d)
    return [
        [int(tet_ids[a]), int(tet_ids[b]), int(tet_ids[c]), pad],
        [int(tet_ids[b]), int(tet_ids[c]), pad, pbd],
        [int(tet_ids[c]), pad, pbd, pcd],
    ]


def emit_scalar_clipped_tet(builder, tet_ids, phi_values, material_a, material_b, zone_id):
    min_phi = float(np.min(phi_values))
    max_phi = float(np.max(phi_values))
    if min_phi >= -EPS:
        return int(builder.add_cell_ids(list(map(int, tet_ids)), VTK_TETRA, material_a, zone_id))
    if max_phi <= EPS:
        return int(builder.add_cell_ids(list(map(int, tet_ids)), VTK_TETRA, material_b, zone_id))

    mat_pair = canonical_material_pair(material_a, material_b)
    emitted = 0
    for tet in _negative_tets_for_scalar(builder, tet_ids, phi_values, mat_pair):
        emitted += int(builder.add_cell_ids(tet, VTK_TETRA, material_b, zone_id))
    for tet in _negative_tets_for_scalar(builder, tet_ids, [-value for value in phi_values], mat_pair):
        emitted += int(builder.add_cell_ids(tet, VTK_TETRA, material_a, zone_id))
    return emitted


def reconstruct_zone_zoo(builder, zone, node_vfs, materials):
    material_a, _ = materials[0]
    material_b, _ = materials[1]
    center_id = builder.get_or_add_center(zone.zone_id, zone.center)
    material_map = dict(materials)
    center_phi = material_map.get(material_a, 0.0) - material_map.get(material_b, 0.0)
    emitted = 0

    for tet_ids in build_consistent_hex_tets(zone.node_ids, center_id):
        phi_values = []
        for point_id in tet_ids:
            if point_id == center_id:
                phi_values.append(center_phi)
            else:
                values = node_vfs[int(point_id)]
                phi_values.append(values.get(material_a, 0.0) - values.get(material_b, 0.0))
        emitted += emit_scalar_clipped_tet(builder, tet_ids, phi_values, material_a, material_b, zone.zone_id)
    return emitted


def reconstruct_zone_fast(builder, zone, materials, material_centroids):
    emitted = 0
    residual_points = zone.points.copy()
    remaining_fraction = 1.0
    for material_id, vf in materials[:-1]:
        local_fraction = vf / remaining_fraction
        normal = choose_split_normal(zone, material_id, material_centroids)
        assigned_points, residual_points = split_rectilinear_hex(residual_points, normal, local_fraction)
        emitted += int(builder.add_cell_points(assigned_points, VTK_HEXAHEDRON, material_id, zone.zone_id))
        remaining_fraction -= vf
    emitted += int(builder.add_cell_points(residual_points, VTK_HEXAHEDRON, materials[-1][0], zone.zone_id))
    return emitted


def _split_tets_by_fraction(tets, normal, fraction):
    if fraction <= EPS:
        return [], tets
    if fraction >= 1.0 - EPS:
        return tets, []

    projections = []
    total_volume = sum(tet_volume(tet) for tet in tets)
    target_volume = fraction * total_volume
    for tet in tets:
        projections.extend(np.dot(tet, normal))
    low = min(projections) - EPS
    high = max(projections) + EPS

    negative = []
    for _ in range(20):
        mid = 0.5 * (low + high)
        negative, volume = clip_tets_volume(tets, normal, mid, keep_negative=True)
        if volume < target_volume:
            low = mid
        else:
            high = mid

    negative, _ = clip_tets_volume(tets, normal, high, keep_negative=True)
    positive, _ = clip_tets_volume(tets, normal, high, keep_negative=False)
    return negative, positive


def reconstruct_zone_plane(builder, mesh, zone, materials, material_centroids):
    emitted = 0
    residual = mesh.zone_tets(zone)
    remaining_fraction = 1.0
    for material_id, vf in materials[:-1]:
        local_fraction = vf / remaining_fraction
        normal = choose_split_normal(zone, material_id, material_centroids)
        assigned, residual = _split_tets_by_fraction(residual, normal, local_fraction)
        for tet in assigned:
            emitted += int(builder.add_cell_points(tet, VTK_TETRA, material_id, zone.zone_id))
        remaining_fraction -= vf
    for tet in residual:
        emitted += int(builder.add_cell_points(tet, VTK_TETRA, materials[-1][0], zone.zone_id))
    return emitted


def reconstruct(mesh, output_path, min_fraction=1.0e-4, max_zones=None, mode="zoo", top_k=2):
    if mode == "fast" and isinstance(mesh, UcdSiloMesh):
        raise ValueError("mode=fast is only supported for rectilinear meshes")

    builder = InMemoryLegacyVtkMesh(mesh.original_points())
    node_vfs = compute_node_vfs(mesh, max_zones=max_zones) if mode == "zoo" else None
    material_centroids = compute_material_centroids(mesh, max_zones=max_zones)
    stats = {"pure_hex_zones": 0, "zoo_zones": 0, "plane_zones": 0, "plane_fallback_zones": 0, "fast_zones": 0}

    for zone in iter_selected_zones(mesh, max_zones=max_zones):
        materials = normalize_materials(zone.materials, min_fraction=min_fraction)
        if not materials:
            continue
        if len(materials) == 1:
            stats["pure_hex_zones"] += 1
            builder.add_cell_ids(zone.node_ids.tolist(), VTK_HEXAHEDRON, materials[0][0], zone.zone_id)
            continue

        if mode == "zoo":
            if len(materials) == 2 and top_k >= 2:
                stats["zoo_zones"] += 1
                reconstruct_zone_zoo(builder, zone, node_vfs, materials)
            else:
                stats["plane_fallback_zones"] += 1
                reconstruct_zone_plane(builder, mesh, zone, materials, material_centroids)
        elif mode == "fast":
            stats["fast_zones"] += 1
            reconstruct_zone_fast(builder, zone, materials, material_centroids)
        else:
            stats["plane_zones"] += 1
            reconstruct_zone_plane(builder, mesh, zone, materials, material_centroids)

    builder.write(output_path)
    return {
        "num_cells": builder.num_cells,
        "num_points": builder.num_points,
        "material_volumes": builder.material_volumes,
        "degenerate_cells": builder.degenerate_cells,
        "stats": stats,
    }


def print_volume_report(reference_volumes, reconstructed_volumes):
    total_reference = sum(reference_volumes.values())
    print("Material volume report:")
    for material_id in sorted(set(reference_volumes) | set(reconstructed_volumes)):
        reference = reference_volumes.get(material_id, 0.0)
        reconstructed = reconstructed_volumes.get(material_id, 0.0)
        rel_error = 0.0 if total_reference <= 0.0 else abs(reconstructed - reference) / total_reference
        print(
            "  material {:>3d}: ref={:.9f} recon={:.9f} rel_total_err={:.6e}".format(
                material_id,
                reference,
                reconstructed,
                rel_error,
            )
        )
