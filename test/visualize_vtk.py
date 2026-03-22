import argparse
from pathlib import Path

import pyvista as pv


def build_parser():
    parser = argparse.ArgumentParser(description="Visualize generated VTK files with PyVista.")
    parser.add_argument("input", help="Path to the VTK file")
    parser.add_argument("--surface", action="store_true", help="Show the extracted outer surface")
    parser.add_argument("--interfaces", action="store_true", help="Extract and show per-material interface surfaces")
    parser.add_argument("--opacity", type=float, default=1.0, help="Mesh opacity in [0, 1]")
    parser.add_argument("--show-edges", action="store_true", help="Show mesh edges")
    parser.add_argument("--slice", dest="slice_view", action="store_true", help="Show an orthogonal slice view")
    parser.add_argument("--screenshot", default=None, help="Optional screenshot output path")
    return parser


def main():
    args = build_parser().parse_args()
    vtk_path = Path(args.input)
    if not vtk_path.exists():
        raise FileNotFoundError("VTK file not found: {}".format(vtk_path))

    mesh = pv.read(vtk_path)
    print(mesh)
    print("Arrays:", mesh.array_names)

    plotter = pv.Plotter(off_screen=bool(args.screenshot))
    plotter.set_background("#f5f1e8", top="#d8e1e8")
    plotter.add_axes()

    if args.interfaces and "material_id" in mesh.array_names:
        material_ids = sorted(int(value) for value in set(mesh.cell_data["material_id"]))
        for material_id in material_ids:
            block = mesh.threshold(
                value=(material_id - 0.5, material_id + 0.5),
                scalars="material_id",
                preference="cell",
            )
            if block.n_cells == 0:
                continue
            surface = block.extract_surface(algorithm="dataset_surface").triangulate()
            plotter.add_mesh(
                surface,
                scalars="material_id" if "material_id" in surface.array_names else None,
                show_edges=args.show_edges,
                opacity=args.opacity,
                cmap="glasbey",
                clim=[min(material_ids), max(material_ids)],
                scalar_bar_args={"title": "material_id"},
                lighting=True,
            )
    else:
        if args.surface:
            display_mesh = mesh.extract_surface(algorithm="dataset_surface").triangulate()
        elif args.slice_view:
            display_mesh = mesh.slice_orthogonal()
        else:
            display_mesh = mesh

        plotter.add_mesh(
            display_mesh,
            scalars="material_id" if "material_id" in mesh.array_names else None,
            show_edges=args.show_edges,
            opacity=args.opacity,
            cmap="glasbey",
            lighting=True,
        )
    plotter.add_text(str(vtk_path.name), font_size=10)

    if args.screenshot:
        plotter.show(screenshot=str(Path(args.screenshot)))
        print("Saved screenshot to {}".format(args.screenshot))
    else:
        plotter.show()


if __name__ == "__main__":
    main()
