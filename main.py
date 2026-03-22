import argparse

from src.mir_mesh import load_mesh
from src.mir_reconstruct import compute_reference_material_volumes, print_volume_report, reconstruct


def build_argument_parser():
    parser = argparse.ArgumentParser(description="Shared-point MIR reconstruction from PDB Silo files.")
    parser.add_argument("input", help="Input Silo file")
    parser.add_argument("output", help="Output legacy VTK unstructured grid path")
    parser.add_argument("--min-fraction", type=float, default=1.0e-4, help="Ignore smaller per-zone material fractions")
    parser.add_argument("--max-zones", type=int, default=None, help="Only process the first N zones")
    parser.add_argument(
        "--mode",
        choices=("zoo", "plane", "fast"),
        default="zoo",
        help="Reconstruction mode: shared-point zoo, plane clipping baseline, or rectilinear fast split",
    )
    parser.add_argument(
        "--rectilinear-mode",
        choices=("zoo", "plane", "fast"),
        default=None,
        help="Deprecated alias for --mode",
    )
    parser.add_argument("--top-k", type=int, default=2, help="Maximum active materials per zone for mode=zoo")
    return parser


def main():
    args = build_argument_parser().parse_args()
    mode = args.rectilinear_mode or args.mode
    if args.top_k < 2:
        raise ValueError("--top-k must be at least 2 for the current zoo implementation")

    mesh = load_mesh(args.input)
    reference_volumes = compute_reference_material_volumes(mesh, max_zones=args.max_zones)
    summary = reconstruct(
        mesh,
        args.output,
        min_fraction=args.min_fraction,
        max_zones=args.max_zones,
        mode=mode,
        top_k=args.top_k,
    )

    print("Wrote {} cells and {} points to {}".format(summary["num_cells"], summary["num_points"], args.output))
    print("Reconstruction stats:", summary["stats"])
    print("Degenerate cells skipped:", summary["degenerate_cells"])
    print_volume_report(reference_volumes, summary["material_volumes"])


if __name__ == "__main__":
    main()
