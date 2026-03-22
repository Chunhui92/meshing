import argparse
import sys
from pathlib import Path

import pyvista as pv


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.mir_mesh import load_mesh
from src.mir_reconstruct import compute_reference_material_volumes
from src.mir_utils import VTK_HEXAHEDRON, VTK_TETRA, hex_volume, tet_volume


def build_parser():
    parser = argparse.ArgumentParser(description="Quantitatively analyze generated VTK files.")
    parser.add_argument("input_vtk", help="Generated VTK file")
    parser.add_argument("--reference-silo", default=None, help="Optional source Silo file for material-volume comparison")
    parser.add_argument("--max-zones", type=int, default=None, help="Restrict the reference Silo comparison to the first N zones")
    parser.add_argument("--expect-cell-types", nargs="*", type=int, default=[VTK_TETRA, VTK_HEXAHEDRON], help="Allowed VTK cell types")
    return parser


def vtk_cell_volume(cell):
    if cell.type == VTK_TETRA:
        return tet_volume(cell.points)
    if cell.type == VTK_HEXAHEDRON:
        return hex_volume(cell.points)
    raise ValueError("Unsupported VTK cell type {}".format(cell.type))


def count_feature_edges(surface, *, boundary=False, non_manifold=False):
    edges = surface.extract_feature_edges(
        boundary_edges=boundary,
        feature_edges=False,
        manifold_edges=False,
        non_manifold_edges=non_manifold,
    )
    return int(edges.n_cells)


def main():
    args = build_parser().parse_args()
    vtk_path = Path(args.input_vtk)
    if not vtk_path.exists():
        raise FileNotFoundError("VTK file not found: {}".format(vtk_path))

    mesh = pv.read(vtk_path)
    print("Readable:", True)
    print("Cells:", mesh.n_cells)
    print("Points:", mesh.n_points)
    print("Arrays:", mesh.array_names)

    cell_types = sorted(int(cell_type) for cell_type in set(mesh.celltypes))
    print("Cell types:", cell_types)
    unexpected = sorted(set(cell_types) - set(args.expect_cell_types))
    if unexpected:
        print("Unexpected cell types:", unexpected)

    if "material_id" not in mesh.cell_data:
        raise ValueError("VTK file does not contain CELL_DATA material_id")

    output_volumes = {}
    for cell_id in range(mesh.n_cells):
        cell = mesh.get_cell(cell_id)
        material_id = int(mesh.cell_data["material_id"][cell_id])
        output_volumes[material_id] = output_volumes.get(material_id, 0.0) + vtk_cell_volume(cell)

    print("Output material volumes:")
    for material_id in sorted(output_volumes):
        print("  {} -> {:.9f}".format(material_id, output_volumes[material_id]))

    if args.reference_silo:
        reference = compute_reference_material_volumes(load_mesh(args.reference_silo), max_zones=args.max_zones)
        total_reference = sum(reference.values())
        print("Reference material volumes:")
        for material_id in sorted(reference):
            print("  {} -> {:.9f}".format(material_id, reference[material_id]))

        print("Relative total-domain errors:")
        for material_id in sorted(set(reference) | set(output_volumes)):
            reference_volume = reference.get(material_id, 0.0)
            output_volume = output_volumes.get(material_id, 0.0)
            rel_error = abs(output_volume - reference_volume) / total_reference if total_reference > 0.0 else 0.0
            print("  {} -> {:.6e}".format(material_id, rel_error))

    print("Surface topology by material:")
    for material_id in sorted(set(int(value) for value in mesh.cell_data["material_id"])):
        block = mesh.threshold(
            value=(material_id - 0.5, material_id + 0.5),
            scalars="material_id",
            preference="cell",
        )
        surface = block.extract_surface(algorithm="dataset_surface").triangulate()
        print(
            "  material {} -> boundary_edges={} non_manifold_edges={}".format(
                material_id,
                count_feature_edges(surface, boundary=True, non_manifold=False),
                count_feature_edges(surface, boundary=False, non_manifold=True),
            )
        )


if __name__ == "__main__":
    main()
