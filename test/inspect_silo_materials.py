import argparse
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.mir_mesh import PDBLiteReader, decode_materials, detect_material_prefix, load_mesh


def build_parser():
    parser = argparse.ArgumentParser(description="Inspect material-count distribution in a Silo mesh.")
    parser.add_argument("input", help="Input Silo file")
    parser.add_argument("--max-zones", type=int, default=None, help="Only inspect the first N zones")
    parser.add_argument("--show-examples", type=int, default=5, help="How many >2-material examples to print")
    return parser


def main():
    args = build_parser().parse_args()
    path = Path(args.input)
    if not path.exists():
        raise FileNotFoundError("Silo file not found: {}".format(path))

    counts = {}
    examples = []
    total = 0
    mode = "mesh"
    try:
        mesh = load_mesh(str(path))
        for zone in mesh.iter_zones():
            if args.max_zones is not None and total >= args.max_zones:
                break
            total += 1
            num_materials = len(zone.materials)
            counts[num_materials] = counts.get(num_materials, 0) + 1
            if num_materials > 2 and len(examples) < args.show_examples:
                examples.append((zone.zone_id, zone.materials))
    except Exception:
        mode = "material-block-only"
        reader = PDBLiteReader(str(path))
        material_prefix = detect_material_prefix(reader)
        if material_prefix is None:
            raise
        matlist = reader.read_array("/{}_matlist".format(material_prefix))
        mix_vf = reader.read_array("/{}_mix_vf".format(material_prefix))
        mix_next = reader.read_array("/{}_mix_next".format(material_prefix))
        mix_mat = reader.read_array("/{}_mix_mat".format(material_prefix))
        total = len(matlist) if args.max_zones is None else min(len(matlist), args.max_zones)
        for zone_id in range(total):
            materials = decode_materials(matlist, mix_vf, mix_next, mix_mat, zone_id)
            num_materials = len(materials)
            counts[num_materials] = counts.get(num_materials, 0) + 1
            if num_materials > 2 and len(examples) < args.show_examples:
                examples.append((zone_id, materials))

    print("File:", path)
    print("Inspection mode:", mode)
    print("Zones inspected:", total)
    print("Material-count distribution:", counts)
    print("Has >2 material zones:", any(key > 2 and value > 0 for key, value in counts.items()))
    if examples:
        print("Examples:")
        for zone_id, materials in examples:
            print("  zone {} -> {}".format(zone_id, materials))


if __name__ == "__main__":
    main()
