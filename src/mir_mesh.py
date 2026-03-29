import re
from dataclasses import dataclass

import numpy as np

from src.mir_utils import tet_volume, tetrahedralize_hex


@dataclass
class Symbol:
    name: str
    dtype_name: str
    count: int
    offset: int


@dataclass
class Zone:
    zone_id: int
    center: np.ndarray
    volume: float
    materials: list
    points: np.ndarray
    node_ids: np.ndarray


class PDBLiteReader:
    TYPE_MAP = {
        "float": np.float32,
        "integer": np.int32,
        "double": np.float64,
        "short": np.int16,
    }

    def __init__(self, path):
        self.path = path
        with open(path, "rb") as handle:
            self.bytes = handle.read()
        self.symbols = self._parse_symbols()

    def _parse_symbols(self):
        text = re.sub(rb"[^ -~]+", b"\n", self.bytes).decode("latin1", errors="ignore")
        tokens = [token.strip() for token in text.splitlines() if token.strip()]
        symbols = {}
        for idx in range(len(tokens) - 5):
            name = tokens[idx]
            dtype_name = tokens[idx + 1]
            if not name.startswith("/") or dtype_name not in self.TYPE_MAP:
                continue
            try:
                count = int(tokens[idx + 2])
                offset = int(tokens[idx + 3])
                marker = int(tokens[idx + 4])
                repeat = int(tokens[idx + 5])
            except ValueError:
                continue
            if marker != 0 or repeat != count:
                continue
            symbols[name] = Symbol(name=name, dtype_name=dtype_name, count=count, offset=offset)
        return symbols

    def has(self, name):
        return name in self.symbols

    def read_array(self, name):
        symbol = self.symbols[name]
        dtype = self.TYPE_MAP[symbol.dtype_name]
        return np.frombuffer(self.bytes, dtype=dtype, count=symbol.count, offset=symbol.offset).copy()


def decode_materials(matlist, mix_vf, mix_next, mix_mat, zone_index):
    code = int(matlist[zone_index])
    if code > 0:
        return [(code, 1.0)]

    materials = []
    cursor = -code - 1
    while True:
        materials.append((int(mix_mat[cursor]), float(mix_vf[cursor])))
        nxt = int(mix_next[cursor])
        if nxt <= 0:
            break
        cursor = nxt - 1
    materials.sort(key=lambda item: item[1], reverse=True)
    return materials


class RectilinearSiloMesh:
    def __init__(self, reader, mesh_prefix, material_prefix):
        self.reader = reader
        self.mesh_prefix = mesh_prefix
        self.material_prefix = material_prefix
        self.x = reader.read_array("/{}_coord0".format(mesh_prefix)).astype(np.float64)
        self.y = reader.read_array("/{}_coord1".format(mesh_prefix)).astype(np.float64)
        self.z = reader.read_array("/{}_coord2".format(mesh_prefix)).astype(np.float64)
        self.matlist = reader.read_array("/{}_matlist".format(material_prefix))
        self.mix_vf = reader.read_array("/{}_mix_vf".format(material_prefix)).astype(np.float64)
        self.mix_next = reader.read_array("/{}_mix_next".format(material_prefix))
        self.mix_mat = reader.read_array("/{}_mix_mat".format(material_prefix))
        self.nz = len(self.z) - 1
        self.ny = len(self.y) - 1
        self.nx = len(self.x) - 1
        self.num_zones = self.nx * self.ny * self.nz
        self.num_nodes = (self.nx + 1) * (self.ny + 1) * (self.nz + 1)
        self.node_points = self._build_node_points()

    def _node_id(self, i, j, k):
        return int(i + (self.nx + 1) * (j + (self.ny + 1) * k))

    def _build_node_points(self):
        points = np.zeros((self.num_nodes, 3), dtype=np.float64)
        for k, z_value in enumerate(self.z):
            for j, y_value in enumerate(self.y):
                for i, x_value in enumerate(self.x):
                    points[self._node_id(i, j, k)] = (x_value, y_value, z_value)
        return points

    def original_points(self):
        return self.node_points

    def iter_zones(self):
        zone_id = 0
        for k in range(self.nz):
            z0, z1 = self.z[k], self.z[k + 1]
            for j in range(self.ny):
                y0, y1 = self.y[j], self.y[j + 1]
                for i in range(self.nx):
                    x0, x1 = self.x[i], self.x[i + 1]
                    points = np.asarray([
                        [x0, y0, z0],
                        [x1, y0, z0],
                        [x1, y1, z0],
                        [x0, y1, z0],
                        [x0, y0, z1],
                        [x1, y0, z1],
                        [x1, y1, z1],
                        [x0, y1, z1],
                    ], dtype=np.float64)
                    node_ids = np.asarray([
                        self._node_id(i, j, k),
                        self._node_id(i + 1, j, k),
                        self._node_id(i + 1, j + 1, k),
                        self._node_id(i, j + 1, k),
                        self._node_id(i, j, k + 1),
                        self._node_id(i + 1, j, k + 1),
                        self._node_id(i + 1, j + 1, k + 1),
                        self._node_id(i, j + 1, k + 1),
                    ], dtype=np.int64)
                    volume = abs((x1 - x0) * (y1 - y0) * (z1 - z0))
                    center = np.mean(points, axis=0)
                    yield Zone(
                        zone_id=zone_id,
                        center=center,
                        volume=volume,
                        materials=decode_materials(self.matlist, self.mix_vf, self.mix_next, self.mix_mat, zone_id),
                        points=points,
                        node_ids=node_ids,
                    )
                    zone_id += 1

    def zone_tets(self, zone):
        return tetrahedralize_hex(zone.points)


class UcdSiloMesh:
    def __init__(self, reader, mesh_prefix, material_prefix, zonelist_prefix):
        self.reader = reader
        self.mesh_prefix = mesh_prefix
        self.material_prefix = material_prefix
        self.zonelist_prefix = zonelist_prefix
        self.x = reader.read_array("/{}_coord0".format(mesh_prefix)).astype(np.float64)
        self.y = reader.read_array("/{}_coord1".format(mesh_prefix)).astype(np.float64)
        self.z = reader.read_array("/{}_coord2".format(mesh_prefix)).astype(np.float64)
        self.points = np.column_stack((self.x, self.y, self.z))
        self.num_nodes = len(self.points)
        self.matlist = reader.read_array("/{}_matlist".format(material_prefix))
        self.mix_vf = reader.read_array("/{}_mix_vf".format(material_prefix)).astype(np.float64)
        self.mix_next = reader.read_array("/{}_mix_next".format(material_prefix))
        self.mix_mat = reader.read_array("/{}_mix_mat".format(material_prefix))
        self.shapesize = int(reader.read_array("/{}_shapesize".format(zonelist_prefix))[0])
        self.shapecnt = int(reader.read_array("/{}_shapecnt".format(zonelist_prefix))[0])
        self.nodelist = reader.read_array("/{}_nodelist".format(zonelist_prefix))
        self.num_zones = self.shapecnt
        if self.shapesize != 8:
            raise ValueError("Only 8-node hexahedral UCD zones are supported, got {}".format(self.shapesize))

    def original_points(self):
        return self.points

    def iter_zones(self):
        for zone_id in range(self.num_zones):
            start = zone_id * self.shapesize
            node_ids = self.nodelist[start:start + self.shapesize].astype(np.int64)
            points = self.points[node_ids]
            volume = sum(tet_volume(tet) for tet in tetrahedralize_hex(points))
            center = np.mean(points, axis=0)
            yield Zone(
                zone_id=zone_id,
                center=center,
                volume=volume,
                materials=decode_materials(self.matlist, self.mix_vf, self.mix_next, self.mix_mat, zone_id),
                points=points,
                node_ids=node_ids,
            )

    def zone_tets(self, zone):
        return tetrahedralize_hex(zone.points)


class MultiBlockUcdSiloMesh:
    def __init__(self, reader, block_specs):
        self.reader = reader
        self.block_specs = list(block_specs)
        self.blocks = []
        self.global_points = {}
        self.num_zones = 0

        for block_index, spec in enumerate(self.block_specs):
            mesh_prefix = spec["mesh_prefix"]
            material_prefix = spec["material_prefix"]
            zonelist_prefix = spec["zonelist_prefix"]

            x = reader.read_array("/{}_coord0".format(mesh_prefix)).astype(np.float64)
            y = reader.read_array("/{}_coord1".format(mesh_prefix)).astype(np.float64)
            z = reader.read_array("/{}_coord2".format(mesh_prefix)).astype(np.float64)
            points = np.column_stack((x, y, z))
            gnodeno_name = "/{}_gnodeno".format(mesh_prefix)
            if reader.has(gnodeno_name):
                global_node_ids = reader.read_array(gnodeno_name).astype(np.int64)
            else:
                base = len(self.global_points)
                global_node_ids = np.arange(base, base + len(points), dtype=np.int64)

            for local_id, global_id in enumerate(global_node_ids):
                global_id = int(global_id)
                point = points[local_id]
                existing = self.global_points.get(global_id)
                if existing is None:
                    self.global_points[global_id] = point.copy()
                elif not np.allclose(existing, point, atol=1.0e-10, rtol=0.0):
                    raise ValueError(
                        "Global node {} has inconsistent coordinates across blocks".format(global_id)
                    )

            self.blocks.append(
                {
                    "block_index": block_index,
                    "mesh_prefix": mesh_prefix,
                    "material_prefix": material_prefix,
                    "zonelist_prefix": zonelist_prefix,
                    "points": points,
                    "global_node_ids": global_node_ids,
                    "matlist": reader.read_array("/{}_matlist".format(material_prefix)),
                    "mix_vf": reader.read_array("/{}_mix_vf".format(material_prefix)).astype(np.float64),
                    "mix_next": reader.read_array("/{}_mix_next".format(material_prefix)),
                    "mix_mat": reader.read_array("/{}_mix_mat".format(material_prefix)),
                    "shapesize": int(reader.read_array("/{}_shapesize".format(zonelist_prefix))[0]),
                    "shapecnt": int(reader.read_array("/{}_shapecnt".format(zonelist_prefix))[0]),
                    "nodelist": reader.read_array("/{}_nodelist".format(zonelist_prefix)).astype(np.int64),
                }
            )
            self.num_zones += self.blocks[-1]["shapecnt"]

        if not self.global_points:
            raise ValueError("No multiblock UCD points were detected")

        self.global_id_to_index = {global_id: idx for idx, global_id in enumerate(sorted(self.global_points))}
        self.points = np.zeros((len(self.global_id_to_index), 3), dtype=np.float64)
        for global_id, idx in self.global_id_to_index.items():
            self.points[idx] = self.global_points[global_id]
        self.num_nodes = len(self.points)

        for block in self.blocks:
            if block["shapesize"] != 8:
                raise ValueError(
                    "Only 8-node hexahedral multiblock UCD zones are supported, got {}".format(block["shapesize"])
                )
            block["node_ids"] = np.asarray(
                [self.global_id_to_index[int(global_id)] for global_id in block["global_node_ids"]],
                dtype=np.int64,
            )

    def original_points(self):
        return self.points

    def iter_zones(self):
        zone_id = 0
        for block in self.blocks:
            shapesize = block["shapesize"]
            for local_zone_id in range(block["shapecnt"]):
                start = local_zone_id * shapesize
                local_ids = block["nodelist"][start:start + shapesize]
                node_ids = block["node_ids"][local_ids]
                points = self.points[node_ids]
                volume = sum(tet_volume(tet) for tet in tetrahedralize_hex(points))
                center = np.mean(points, axis=0)
                yield Zone(
                    zone_id=zone_id,
                    center=center,
                    volume=volume,
                    materials=decode_materials(
                        block["matlist"],
                        block["mix_vf"],
                        block["mix_next"],
                        block["mix_mat"],
                        local_zone_id,
                    ),
                    points=points,
                    node_ids=node_ids,
                )
                zone_id += 1

    def zone_tets(self, zone):
        return tetrahedralize_hex(zone.points)


def detect_rectilinear_mesh(reader):
    for symbol in sorted(reader.symbols):
        if symbol.endswith("_coord0"):
            prefix = symbol[1:-7]
            if reader.has("/{}_coord1".format(prefix)) and reader.has("/{}_coord2".format(prefix)):
                return prefix
    return None


def detect_material_prefix(reader):
    for symbol in sorted(reader.symbols):
        if symbol.endswith("_matlist"):
            prefix = symbol[1:-8]
            needed = (
                "/{}_mix_vf".format(prefix),
                "/{}_mix_next".format(prefix),
                "/{}_mix_mat".format(prefix),
            )
            if all(reader.has(name) for name in needed):
                return prefix
    return None


def detect_zonelist_prefix(reader):
    for symbol in sorted(reader.symbols):
        if symbol.endswith("_nodelist"):
            prefix = symbol[1:-9]
            needed = (
                "/{}_shapesize".format(prefix),
                "/{}_shapecnt".format(prefix),
            )
            if all(reader.has(name) for name in needed) and prefix.startswith("zl"):
                return prefix
    return None


def detect_multiblock_specs(reader):
    block_names = sorted(
        {
            match.group(1)
            for symbol in reader.symbols
            for match in [re.match(r"^/(block\d+)/", symbol)]
            if match is not None
        }
    )
    block_specs = []
    for block_name in block_names:
        mesh_prefix = "{}/mesh1".format(block_name)
        material_prefix = "{}/mat1".format(block_name)
        zonelist_prefix = "{}/zl1".format(block_name)
        needed = (
            "/{}_coord0".format(mesh_prefix),
            "/{}_coord1".format(mesh_prefix),
            "/{}_coord2".format(mesh_prefix),
            "/{}_matlist".format(material_prefix),
            "/{}_mix_vf".format(material_prefix),
            "/{}_mix_next".format(material_prefix),
            "/{}_mix_mat".format(material_prefix),
            "/{}_nodelist".format(zonelist_prefix),
            "/{}_shapesize".format(zonelist_prefix),
            "/{}_shapecnt".format(zonelist_prefix),
        )
        if all(reader.has(name) for name in needed):
            block_specs.append(
                {
                    "mesh_prefix": mesh_prefix,
                    "material_prefix": material_prefix,
                    "zonelist_prefix": zonelist_prefix,
                }
            )
    return block_specs


def load_mesh(path):
    reader = PDBLiteReader(path)
    multiblock_specs = detect_multiblock_specs(reader)
    if multiblock_specs:
        return MultiBlockUcdSiloMesh(reader, multiblock_specs)

    material_prefix = detect_material_prefix(reader)
    if material_prefix is None:
        raise ValueError("No Silo material block was found in {}".format(path))

    zonelist_prefix = detect_zonelist_prefix(reader)
    if zonelist_prefix is not None:
        mesh_prefix = None
        for symbol in sorted(reader.symbols):
            if symbol.endswith("_coord0") and not symbol.startswith("/exterior_faces"):
                mesh_prefix = symbol[1:-7]
                if reader.has("/{}_coord1".format(mesh_prefix)) and reader.has("/{}_coord2".format(mesh_prefix)):
                    break
        if mesh_prefix is None:
            raise ValueError("Could not find a UCD mesh coordinate block")
        return UcdSiloMesh(reader, mesh_prefix, material_prefix, zonelist_prefix)

    mesh_prefix = detect_rectilinear_mesh(reader)
    if mesh_prefix is None:
        raise ValueError("Could not detect a supported mesh in {}".format(path))
    return RectilinearSiloMesh(reader, mesh_prefix, material_prefix)
