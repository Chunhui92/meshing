import argparse

from src.mir_mesh import load_mesh
from src.mir_reconstruct import compute_reference_material_volumes, print_volume_report, reconstruct


def build_argument_parser():
    parser = argparse.ArgumentParser(description="Zoo-only shared-point MIR reconstruction from PDB Silo files.")
    parser.add_argument("input", help="Input Silo file")
    parser.add_argument("output", help="Output legacy VTK unstructured grid path")
    parser.add_argument("--min-fraction", type=float, default=1.0e-4, help="Ignore smaller per-zone material fractions")
    parser.add_argument("--max-zones", type=int, default=None, help="Only process the first N zones")
    parser.add_argument("--mode", choices=("zoo",), default="zoo", help="Reconstruction mode. Only zoo is supported.")
    parser.add_argument(
        "--volume-correction-iters",
        type=int,
        default=24,
        help="Maximum number of per-split bisection iterations used to match target material volume fractions",
    )
    parser.add_argument(
        "--volume-correction-tol",
        type=float,
        default=1.0e-4,
        help="Relative volume tolerance for per-split zoo volume correction",
    )
    return parser


def main():
    args = build_argument_parser().parse_args()
    mesh = load_mesh(args.input)
    reference_volumes = compute_reference_material_volumes(mesh, max_zones=args.max_zones)
    summary = reconstruct(
        mesh,
        args.output,
        min_fraction=args.min_fraction,
        max_zones=args.max_zones,
        mode=args.mode,
        volume_correction_iters=args.volume_correction_iters,
        volume_correction_tol=args.volume_correction_tol,
    )

    print("Wrote {} cells and {} points to {}".format(summary["num_cells"], summary["num_points"], args.output))
    print("Reconstruction stats:", summary["stats"])
    print("Degenerate cells skipped:", summary["degenerate_cells"])
    print_volume_report(reference_volumes, summary["material_volumes"])


if __name__ == "__main__":
    main()
