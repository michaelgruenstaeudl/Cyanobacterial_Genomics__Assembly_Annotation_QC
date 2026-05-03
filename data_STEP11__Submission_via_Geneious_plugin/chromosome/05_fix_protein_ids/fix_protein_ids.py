#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

DB = "HT2024"

cds_pattern = re.compile(
    r'product\s+whole\s+local\s+str\s+"([^"]+)".*?locus-tag\s+"([^"]+)"',
    re.DOTALL,
)

geneious_id_pattern = re.compile(
    r"^Limnothrix_sp_HT2024_(?:chromosome|plasmid)\d+$"
)

bad_locus_tag_pattern = re.compile(r"^pgaptmp[_-]?", re.IGNORECASE)


def make_safe_tag(locus_tag: str, counter: int, locus_tag_prefix: str | None) -> str:
    """
    Use the existing locus_tag unless it starts with pgaptmp.
    If locus_tag_prefix is provided, use that prefix instead.
    """
    if locus_tag_prefix:
        return f"{locus_tag_prefix}_{counter:06d}"

    if bad_locus_tag_pattern.match(locus_tag):
        return f"HT2024_{counter:06d}"

    return locus_tag


def fix_file(infile: str, locus_tag_prefix: str | None) -> None:
    in_path = Path(infile)

    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {infile}")

    out_path = in_path.with_name(
        f"{in_path.stem}.protein_ids_fixed{in_path.suffix}"
    )

    text = in_path.read_text()

    id_map = {}
    counter = 1

    for old_id, locus_tag in cds_pattern.findall(text):
        if geneious_id_pattern.match(old_id):
            new_tag = make_safe_tag(locus_tag, counter, locus_tag_prefix)
            id_map[old_id] = new_tag
            counter += 1

    if not id_map:
        raise RuntimeError(
            f"No Geneious protein IDs with locus-tags found in {infile}"
        )

    for old_id in sorted(id_map, key=len, reverse=True):
        new_tag = id_map[old_id]

        replacement = (
            'general {\n'
            f'                    db "{DB}",\n'
            f'                    tag str "{new_tag}"\n'
            '                  }'
        )

        text = re.sub(
            rf'local\s+str\s+"{re.escape(old_id)}"',
            replacement,
            text,
        )

    out_path.write_text(text)

    print(f"Wrote: {out_path}")
    print(f"Replaced {len(id_map)} protein IDs")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Replace Geneious local protein IDs with GenBank-compatible general protein IDs."
    )

    parser.add_argument("infile", help="Input ASN/SQN file")
    parser.add_argument(
        "--locus-tag-prefix",
        help="Optional official locus_tag prefix, for example ABCD01 or HT2024",
        default=None,
    )

    args = parser.parse_args()

    if args.locus_tag_prefix and args.locus_tag_prefix.lower().startswith("pgaptmp"):
        raise ValueError("locus_tag_prefix must not start with pgaptmp")

    fix_file(args.infile, args.locus_tag_prefix)


if __name__ == "__main__":
    main()
