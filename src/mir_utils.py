import itertools
import math

import numpy as np


VTK_TETRA = 10
VTK_HEXAHEDRON = 12
EPS = 1.0e-9
CLAMP_EPS = 1.0e-8

HEX_TO_TETS = (
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
)
HEX_FACES = (
    (0, 1, 2, 3),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)
TET_EDGES = (
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 2),
    (1, 3),
    (2, 3),
)


def tet_volume(tet):
    a, b, c, d = np.asarray(tet, dtype=np.float64)
    return abs(np.dot(b - a, np.cross(c - a, d - a))) / 6.0


def tetrahedralize_hex(points):
    points = np.asarray(points, dtype=np.float64)
    return [points[np.array(ids, dtype=np.int64)] for ids in HEX_TO_TETS]


def hex_volume(points):
    return sum(tet_volume(tet) for tet in tetrahedralize_hex(points))


def unique_points(points, eps=EPS):
    kept = []
    for point in points:
        if not any(np.linalg.norm(point - other) <= eps for other in kept):
            kept.append(point)
    return kept


def plane_basis(normal):
    normal = np.asarray(normal, dtype=np.float64)
    normal /= np.linalg.norm(normal)
    helper = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    if abs(np.dot(helper, normal)) > 0.9:
        helper = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    u = np.cross(normal, helper)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    return u, v


def hull_triangles(points, eps=EPS):
    points = [np.asarray(point, dtype=np.float64) for point in points]
    if len(points) < 4:
        return []

    planes = {}
    for i, j, k in itertools.combinations(range(len(points)), 3):
        pi, pj, pk = points[i], points[j], points[k]
        normal = np.cross(pj - pi, pk - pi)
        norm = np.linalg.norm(normal)
        if norm <= eps:
            continue
        normal /= norm
        distances = [np.dot(normal, point - pi) for point in points]
        if any(value > eps for value in distances) and any(value < -eps for value in distances):
            continue
        outward = normal if max(distances) <= eps else -normal
        offset = np.dot(outward, pi)
        key = tuple(np.round(np.concatenate((outward, [offset])), 8))
        planes[key] = (outward, offset)

    triangles = []
    for normal, offset in planes.values():
        coplanar = [
            idx for idx, point in enumerate(points)
            if abs(np.dot(normal, point) - offset) <= 1.0e-7
        ]
        if len(coplanar) < 3:
            continue
        plane_center = np.mean([points[idx] for idx in coplanar], axis=0)
        u, v = plane_basis(normal)
        ordered = []
        for idx in coplanar:
            rel = points[idx] - plane_center
            ordered.append((math.atan2(np.dot(rel, v), np.dot(rel, u)), idx))
        ordered.sort()
        face = [idx for _, idx in ordered]
        for n in range(1, len(face) - 1):
            triangles.append((face[0], face[n], face[n + 1]))
    return triangles


def clip_tet_with_plane(tet, normal, offset, keep_negative=True):
    tet = np.asarray(tet, dtype=np.float64)
    signed = np.dot(tet, normal) - offset
    if not keep_negative:
        signed = -signed

    kept = []
    for idx, value in enumerate(signed):
        if value <= EPS:
            kept.append(tet[idx])
    for i, j in TET_EDGES:
        di, dj = signed[i], signed[j]
        if (di < -EPS and dj > EPS) or (di > EPS and dj < -EPS):
            t = di / (di - dj)
            kept.append(tet[i] + t * (tet[j] - tet[i]))

    kept = unique_points(kept)
    if len(kept) < 4:
        return []

    triangles = hull_triangles(kept)
    if not triangles:
        return []

    centroid = np.mean(np.asarray(kept), axis=0)
    result = []
    for a, b, c in triangles:
        poly_tet = np.asarray([centroid, kept[a], kept[b], kept[c]], dtype=np.float64)
        if tet_volume(poly_tet) > 1.0e-12:
            result.append(poly_tet)
    return result


def clip_tets_volume(tets, normal, offset, keep_negative=True):
    clipped = []
    total = 0.0
    for tet in tets:
        pieces = clip_tet_with_plane(tet, normal, offset, keep_negative=keep_negative)
        clipped.extend(pieces)
        total += sum(tet_volume(piece) for piece in pieces)
    return clipped, total


def choose_split_normal(zone, material_id, material_centroids):
    target = material_centroids.get(material_id)
    if target is not None:
        normal = zone.center - target
        norm = np.linalg.norm(normal)
        if norm > 1.0e-12:
            return normal / norm

    spans = np.ptp(zone.points, axis=0)
    axis = int(np.argmax(spans))
    normal = np.zeros(3, dtype=np.float64)
    normal[axis] = 1.0
    return normal


def split_rectilinear_hex(points, normal, fraction):
    points = np.asarray(points, dtype=np.float64)
    axis = int(np.argmax(np.abs(normal)))
    coord_min = np.min(points[:, axis])
    coord_max = np.max(points[:, axis])
    span = coord_max - coord_min
    if span <= 1.0e-12:
        return points.copy(), points.copy()

    take_low_side = normal[axis] >= 0.0
    cut = coord_min + fraction * span if take_low_side else coord_max - fraction * span
    assigned = points.copy()
    residual = points.copy()

    low_ids = [idx for idx, point in enumerate(points) if abs(point[axis] - coord_min) <= 1.0e-12]
    high_ids = [idx for idx, point in enumerate(points) if abs(point[axis] - coord_max) <= 1.0e-12]

    if take_low_side:
        for idx in high_ids:
            assigned[idx, axis] = cut
        for idx in low_ids:
            residual[idx, axis] = cut
    else:
        for idx in low_ids:
            assigned[idx, axis] = cut
        for idx in high_ids:
            residual[idx, axis] = cut
    return assigned, residual
