from dataclasses import dataclass

import numpy as np

from src.mir_utils import CLAMP_EPS, EPS, HEX_FACES, VTK_HEXAHEDRON, VTK_TETRA, hex_volume, hull_triangles, tet_volume


@dataclass
class TetFragment:
    point_ids: tuple


@dataclass
class ZonePoint:
    coord: np.ndarray
    values: np.ndarray
    kind: str
    output_point_id: int | None = None


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

    def get_or_add_edge_point(self, v0, v1, clip_context, t):
        key = (min(v0, v1), max(v0, v1), clip_context)
        point_id = self.edge_point_cache.get(key)
        if point_id is not None:
            return point_id

        t = min(max(float(t), CLAMP_EPS), 1.0 - CLAMP_EPS)
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


class ZoneGeometryContext:
    def __init__(self, builder, zone, node_vfs, materials):
        self.builder = builder
        self.zone = zone
        self.material_ids = [int(material_id) for material_id, _ in materials]
        self.material_index = {material_id: idx for idx, material_id in enumerate(self.material_ids)}
        self.zone_material_values = np.asarray([float(vf) for _, vf in materials], dtype=np.float64)
        self.points = {}
        self.next_point_id = 0
        self.local_edge_cache = {}
        self.helper_cache = {}

        self.corner_local_ids = []
        for node_id in zone.node_ids:
            values = self._project_node_values(node_vfs[int(node_id)])
            local_id = self._add_point(zone.points[len(self.corner_local_ids)], values, "corner", output_point_id=int(node_id))
            self.corner_local_ids.append(local_id)

        center_output_id = builder.get_or_add_center(zone.zone_id, zone.center)
        center_id = self._add_point(zone.center, self.zone_material_values.copy(), "helper", output_point_id=center_output_id)
        self.center_local_id = center_id

    def _project_node_values(self, values):
        projected = np.zeros(len(self.material_ids), dtype=np.float64)
        for material_id, idx in self.material_index.items():
            projected[idx] = float(values.get(material_id, 0.0))
        total = float(np.sum(projected))
        if total <= EPS:
            return self.zone_material_values.copy()
        return projected / total

    def _add_point(self, coord, values, kind, output_point_id=None):
        point_id = self.next_point_id
        self.next_point_id += 1
        self.points[point_id] = ZonePoint(
            coord=np.asarray(coord, dtype=np.float64),
            values=np.asarray(values, dtype=np.float64),
            kind=kind,
            output_point_id=output_point_id,
        )
        return point_id

    def point_coord(self, point_id):
        return self.points[int(point_id)].coord

    def point_values(self, point_id):
        return self.points[int(point_id)].values

    def materialize_point(self, point_id):
        point = self.points[int(point_id)]
        if point.output_point_id is None:
            point.output_point_id = self.builder.add_point(point.coord)
        return point.output_point_id

    def get_or_add_intersection(self, point_a, point_b, phi_a, phi_b, clip_context):
        edge = (min(int(point_a), int(point_b)), max(int(point_a), int(point_b)), clip_context)
        cached = self.local_edge_cache.get(edge)
        if cached is not None:
            return cached

        denom = float(phi_a - phi_b)
        if abs(denom) <= EPS:
            t = 0.5
        else:
            t = float(phi_a / denom)
        t = min(max(t, CLAMP_EPS), 1.0 - CLAMP_EPS)

        point0 = self.points[int(point_a)]
        point1 = self.points[int(point_b)]
        coord = point0.coord + t * (point1.coord - point0.coord)
        values = point0.values + t * (point1.values - point0.values)

        output_point_id = None
        if point0.output_point_id is not None and point1.output_point_id is not None:
            output_point_id = self.builder.get_or_add_edge_point(point0.output_point_id, point1.output_point_id, clip_context, t)

        point_id = self._add_point(coord, values, "edge", output_point_id=output_point_id)
        self.local_edge_cache[edge] = point_id
        return point_id

    def get_or_add_helper(self, point_ids, clip_context):
        key = (clip_context, tuple(sorted(int(point_id) for point_id in point_ids)))
        cached = self.helper_cache.get(key)
        if cached is not None:
            return cached

        coords = np.asarray([self.point_coord(point_id) for point_id in point_ids], dtype=np.float64)
        values = np.asarray([self.point_values(point_id) for point_id in point_ids], dtype=np.float64)
        helper_id = self._add_point(np.mean(coords, axis=0), np.mean(values, axis=0), "helper")
        self.helper_cache[key] = helper_id
        return helper_id


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
        return [], len(materials)
    total = sum(vf for _, vf in active)
    if total <= 0.0:
        return [], len(materials)
    normalized = [(material_id, vf / total) for material_id, vf in active]
    normalized.sort(key=lambda item: (-item[1], item[0]))
    return normalized, len(materials) - len(normalized)


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
            tets.append(TetFragment(point_ids=(center_id, zone_node_ids[tri[0]], zone_node_ids[tri[1]], zone_node_ids[tri[2]])))
    return tets


def fragment_volume(context, fragment):
    coords = np.asarray([context.point_coord(point_id) for point_id in fragment.point_ids], dtype=np.float64)
    return tet_volume(coords)


def total_fragment_volume(context, fragments):
    return sum(fragment_volume(context, fragment) for fragment in fragments)


def phi_for_material_split(context, point_id, target_material, remaining_materials):
    values = context.point_values(point_id)
    target = values[context.material_index[target_material]]
    residual = sum(values[context.material_index[material_id]] for material_id in remaining_materials if material_id != target_material)
    return float(target - residual)


def _classify_vertex(phi_value, offset, keep_positive):
    signed = float(phi_value) - float(offset)
    return signed >= -EPS if keep_positive else signed <= EPS


def _crosses_level(phi_a, phi_b, offset):
    return (phi_a - offset < -EPS and phi_b - offset > EPS) or (phi_a - offset > EPS and phi_b - offset < -EPS)


def _unique_point_ids(context, point_ids):
    unique = []
    coords = []
    for point_id in point_ids:
        coord = context.point_coord(point_id)
        duplicate = False
        for other in coords:
            if np.linalg.norm(coord - other) <= 1.0e-10:
                duplicate = True
                break
        if not duplicate:
            unique.append(int(point_id))
            coords.append(coord)
    return unique


def clip_tet_fragment(context, fragment, phi_values, offset, keep_positive, clip_context):
    kept = []
    tet_ids = fragment.point_ids
    for local_idx, point_id in enumerate(tet_ids):
        if _classify_vertex(phi_values[local_idx], offset, keep_positive):
            kept.append(int(point_id))

    for i in range(4):
        for j in range(i + 1, 4):
            phi_i = float(phi_values[i])
            phi_j = float(phi_values[j])
            if _crosses_level(phi_i, phi_j, offset):
                kept.append(
                    context.get_or_add_intersection(
                        tet_ids[i],
                        tet_ids[j],
                        phi_i - offset,
                        phi_j - offset,
                        clip_context,
                    )
                )

    kept = _unique_point_ids(context, kept)
    if len(kept) < 4:
        return []

    coords = np.asarray([context.point_coord(point_id) for point_id in kept], dtype=np.float64)
    triangles = hull_triangles(coords)
    if not triangles:
        return []

    centroid_id = context.get_or_add_helper(kept, (clip_context, "centroid"))
    clipped = []
    for tri in sorted(triangles):
        tet_ids_out = (centroid_id, kept[tri[0]], kept[tri[1]], kept[tri[2]])
        if len(set(tet_ids_out)) < 4:
            continue
        tet_coords = np.asarray([context.point_coord(point_id) for point_id in tet_ids_out], dtype=np.float64)
        if tet_volume(tet_coords) <= 1.0e-12:
            continue
        clipped.append(TetFragment(point_ids=tet_ids_out))
    return clipped


def clip_fragments(context, fragments, target_material, remaining_materials, offset, clip_context):
    assigned = []
    residual = []
    assigned_volume = 0.0
    for fragment in fragments:
        phi_values = [phi_for_material_split(context, point_id, target_material, remaining_materials) for point_id in fragment.point_ids]
        target_pieces = clip_tet_fragment(context, fragment, phi_values, offset, True, (clip_context, "target"))
        residual_pieces = clip_tet_fragment(context, fragment, phi_values, offset, False, (clip_context, "residual"))
        assigned.extend(target_pieces)
        residual.extend(residual_pieces)
        assigned_volume += sum(fragment_volume(context, piece) for piece in target_pieces)
    return assigned, residual, assigned_volume


def solve_split_offset(context, fragments, target_material, remaining_materials, target_fraction, max_iters, volume_tol, offset_tol, clip_context):
    total_volume = total_fragment_volume(context, fragments)
    if total_volume <= EPS:
        return 0.0, [], [], 0.0, 0, True

    target_volume = target_fraction * total_volume
    phi_values = [
        phi_for_material_split(context, point_id, target_material, remaining_materials)
        for fragment in fragments
        for point_id in fragment.point_ids
    ]
    low = min(phi_values) - 1.0
    high = max(phi_values) + 1.0

    best = None
    best_error = float("inf")
    iterations = 0
    for iterations in range(1, max_iters + 1):
        mid = 0.5 * (low + high)
        assigned, residual, assigned_volume = clip_fragments(
            context,
            fragments,
            target_material,
            remaining_materials,
            mid,
            (clip_context, iterations),
        )
        error = assigned_volume - target_volume
        abs_error = abs(error)
        if abs_error < best_error:
            best = (mid, assigned, residual, assigned_volume)
            best_error = abs_error
        if abs_error <= max(volume_tol * total_volume, 1.0e-12) or abs(high - low) <= offset_tol:
            break
        if error > 0.0:
            low = mid
        else:
            high = mid

    success = best_error <= max(volume_tol * total_volume, 1.0e-12)
    offset, assigned, residual, assigned_volume = best
    return offset, assigned, residual, assigned_volume, iterations, success


def emit_fragments(builder, context, fragments, material_id, zone_id):
    emitted = 0
    for fragment in fragments:
        point_ids = [context.materialize_point(point_id) for point_id in fragment.point_ids]
        emitted += int(builder.add_cell_ids(point_ids, VTK_TETRA, material_id, zone_id))
    return emitted


def reconstruct_zone_zoo(builder, zone, node_vfs, materials, correction_iters, correction_tol, stats):
    context = ZoneGeometryContext(builder, zone, node_vfs, materials)
    fragments = build_consistent_hex_tets(context.corner_local_ids, context.center_local_id)
    remaining_materials = [material_id for material_id, _ in materials]
    emitted = 0

    for split_index, (material_id, vf) in enumerate(materials[:-1]):
        current_total = sum(weight for current_material, weight in materials[split_index:])
        local_fraction = vf / current_total if current_total > EPS else 0.0
        offset, assigned, fragments, assigned_volume, iterations, success = solve_split_offset(
            context,
            fragments,
            material_id,
            remaining_materials,
            local_fraction,
            correction_iters,
            correction_tol,
            1.0e-10,
            (zone.zone_id, material_id, tuple(remaining_materials)),
        )
        stats["volume_correction_iterations"] += iterations
        if not success:
            stats["volume_correction_failed_zones"] += 1
        if abs(offset) > EPS or not success:
            stats["volume_correction_adjusted_splits"] += 1
        emitted += emit_fragments(builder, context, assigned, material_id, zone.zone_id)
        remaining_materials = remaining_materials[1:]

    if remaining_materials and fragments:
        emitted += emit_fragments(builder, context, fragments, remaining_materials[0], zone.zone_id)
    return emitted


def reconstruct(mesh, output_path, min_fraction=1.0e-4, max_zones=None, mode="zoo", top_k=None, volume_correction_iters=24, volume_correction_tol=1.0e-4):
    if mode != "zoo":
        raise ValueError("Only mode=zoo is supported in the current implementation")
    if top_k is not None:
        raise ValueError("--top-k is no longer supported in zoo-only mode")

    builder = InMemoryLegacyVtkMesh(mesh.original_points())
    node_vfs = compute_node_vfs(mesh, max_zones=max_zones)
    stats = {
        "pure_zones": 0,
        "zoo_2mat_zones": 0,
        "zoo_nmat_zones": 0,
        "filtered_small_material_zones": 0,
        "volume_correction_iterations": 0,
        "volume_correction_failed_zones": 0,
        "volume_correction_adjusted_splits": 0,
        "degenerate_fragments": 0,
    }

    for zone in iter_selected_zones(mesh, max_zones=max_zones):
        materials, filtered_count = normalize_materials(zone.materials, min_fraction=min_fraction)
        if filtered_count > 0:
            stats["filtered_small_material_zones"] += 1
        if not materials:
            continue
        if len(materials) == 1:
            stats["pure_zones"] += 1
            builder.add_cell_ids(zone.node_ids.tolist(), VTK_HEXAHEDRON, materials[0][0], zone.zone_id)
            continue

        if len(materials) == 2:
            stats["zoo_2mat_zones"] += 1
        else:
            stats["zoo_nmat_zones"] += 1
        reconstruct_zone_zoo(
            builder,
            zone,
            node_vfs,
            materials,
            correction_iters=volume_correction_iters,
            correction_tol=volume_correction_tol,
            stats=stats,
        )

    stats["degenerate_fragments"] = builder.degenerate_cells
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
