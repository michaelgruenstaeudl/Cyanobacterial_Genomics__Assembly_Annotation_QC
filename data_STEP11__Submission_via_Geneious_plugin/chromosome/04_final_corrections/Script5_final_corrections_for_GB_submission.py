#!/usr/bin/env python3

# Script5_final_corrections_for_GB_submission.py

import argparse
import logging
import re
import sys
from collections import Counter

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, CompoundLocation


PRODUCT_REPLACEMENTS = {
    "Tll0287-like domain-containing protein": "DUF433 domain-containing protein",
    "domain-containing protein": "DUF433 domain-containing protein",
    "slr1658 superfamily regulator": "superfamily regulator",
    "slr1659 superfamily regulator": "superfamily regulator",
    "Npun_R2821/Npun_R2822 family protein": "Npun family protein",
    "Npun R2821/Npun R2822 family protein": "Npun family protein",
}


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def same_location(a, b):
    return str(a.location) == str(b.location)


def first_qualifier(feature, key):
    values = feature.qualifiers.get(key)
    if values:
        return values[0]
    return None


def feature_interval(feature):
    return int(feature.location.start), int(feature.location.end)


def features_overlap(feature_a, feature_b):
    start_a, end_a = feature_interval(feature_a)
    start_b, end_b = feature_interval(feature_b)
    return start_a < end_b and start_b < end_a


def is_joined_feature(feature):
    return isinstance(feature.location, CompoundLocation)


def sanitize_for_protein_id(value):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", value)


def infer_locus_tag_prefix(record):
    prefixes = []

    for feature in record.features:
        locus_tag = first_qualifier(feature, "locus_tag")
        if not locus_tag:
            continue

        match = re.match(r"^(.+?)[_-](\d+)$", locus_tag)
        if match:
            prefixes.append(match.group(1))

    if prefixes:
        return Counter(prefixes).most_common(1)[0][0]

    return "generated"


def existing_locus_tags(record):
    tags = set()

    for feature in record.features:
        locus_tag = first_qualifier(feature, "locus_tag")
        if locus_tag:
            tags.add(locus_tag)

    return tags


def next_generated_locus_tag(record, prefix, used_tags):
    max_number = 0
    pattern = re.compile(rf"^{re.escape(prefix)}[_-](\d+)$")

    for tag in used_tags:
        match = pattern.match(tag)
        if match:
            max_number = max(max_number, int(match.group(1)))

    candidate_number = max_number + 1

    while True:
        candidate = f"{prefix}_{candidate_number:06d}"
        if candidate not in used_tags:
            used_tags.add(candidate)
            return candidate
        candidate_number += 1


def add_missing_gene_features(record):
    existing_genes = [f for f in record.features if f.type == "gene"]
    additions = []
    added_count = 0

    for feature in record.features:
        if feature.type not in {"ncRNA", "tmRNA"}:
            continue

        gene_name = first_qualifier(feature, "gene")
        if not gene_name:
            continue

        found = any(
            same_location(feature, gene_feature)
            for gene_feature in existing_genes
        )

        if not found:
            qualifiers = {"gene": [gene_name]}

            locus_tag = first_qualifier(feature, "locus_tag")
            if locus_tag:
                qualifiers["locus_tag"] = [locus_tag]

            new_gene = SeqFeature(
                location=feature.location,
                type="gene",
                qualifiers=qualifiers,
            )
            additions.append(new_gene)
            existing_genes.append(new_gene)
            added_count += 1
            logging.info(
                "Added gene feature for %s at %s",
                feature.type,
                feature.location,
            )

    record.features.extend(additions)
    return added_count


def remove_pseudo_cds(record):
    kept_features = []
    removed = 0

    for feature in record.features:
        if feature.type == "CDS" and "pseudo" in feature.qualifiers:
            removed += 1
            logging.info(
                "Removed pseudo CDS at %s; locus_tag=%s; product=%s",
                feature.location,
                first_qualifier(feature, "locus_tag"),
                first_qualifier(feature, "product"),
            )
            continue

        kept_features.append(feature)

    record.features = kept_features
    return removed


def remove_all_joined_or_ribosomal_slippage_cds(record):
    kept_features = []
    removed = 0

    for feature in record.features:
        if feature.type == "CDS" and (
            is_joined_feature(feature)
            or "ribosomal_slippage" in feature.qualifiers
        ):
            removed += 1
            logging.info(
                "Removed joined/ribosomal-slippage CDS at %s; locus_tag=%s; product=%s",
                feature.location,
                first_qualifier(feature, "locus_tag"),
                first_qualifier(feature, "product"),
            )
            continue

        kept_features.append(feature)

    record.features = kept_features
    return removed


def remove_trna_features_overlapping_cds(record):
    cds_features = [feature for feature in record.features if feature.type == "CDS"]

    kept_features = []
    removed_trna = 0
    removed_empty_companion_genes = 0
    removed_gene_locations = set()

    for feature in record.features:
        if feature.type == "tRNA" and any(
            features_overlap(feature, cds) for cds in cds_features
        ):
            removed_trna += 1
            removed_gene_locations.add(str(feature.location))
            logging.info(
                "Removed tRNA overlapping CDS at %s; product=%s; inference=%s",
                feature.location,
                first_qualifier(feature, "product"),
                first_qualifier(feature, "inference"),
            )
            continue

        kept_features.append(feature)

    final_features = []

    for feature in kept_features:
        if feature.type == "gene" and str(feature.location) in removed_gene_locations:
            has_gene = bool(feature.qualifiers.get("gene"))
            has_locus_tag = bool(feature.qualifiers.get("locus_tag"))
            has_locus = bool(feature.qualifiers.get("locus"))

            if not has_gene and not has_locus_tag and not has_locus:
                removed_empty_companion_genes += 1
                logging.info(
                    "Removed empty companion gene for removed tRNA at %s",
                    feature.location,
                )
                continue

        final_features.append(feature)

    record.features = final_features
    return removed_trna, removed_empty_companion_genes


def clean_qualifiers(record):
    gene_features = [feature for feature in record.features if feature.type == "gene"]

    removed_dbxref = 0
    added_products = 0
    removed_notes = 0
    added_genes_to_cds = 0

    for feature in record.features:
        if feature.type != "source" and "db_xref" in feature.qualifiers:
            removed_dbxref += len(feature.qualifiers["db_xref"])
            del feature.qualifiers["db_xref"]

        keep_note = (
            feature.type == "misc_feature"
            or (feature.type == "CDS" and "EC_number" in feature.qualifiers)
        )

        if not keep_note and "note" in feature.qualifiers:
            removed_notes += len(feature.qualifiers["note"])
            del feature.qualifiers["note"]

        if feature.type == "CDS":
            if "gene" not in feature.qualifiers or not feature.qualifiers["gene"]:
                for gene_feature in gene_features:
                    if (
                        same_location(feature, gene_feature)
                        and "gene" in gene_feature.qualifiers
                        and gene_feature.qualifiers["gene"]
                    ):
                        feature.qualifiers["gene"] = gene_feature.qualifiers["gene"][:]
                        added_genes_to_cds += 1
                        logging.info("Added /gene= to CDS at %s", feature.location)
                        break

            if (
                "product" not in feature.qualifiers
                or not feature.qualifiers["product"]
                or not feature.qualifiers["product"][0].strip()
            ):
                gene_name = first_qualifier(feature, "gene")

                if gene_name:
                    feature.qualifiers["product"] = [gene_name]
                else:
                    feature.qualifiers["product"] = ["hypothetical protein"]

                added_products += 1
                logging.info("Added /product= to CDS at %s", feature.location)

    return removed_dbxref, added_products, removed_notes, added_genes_to_cds


def ensure_cds_protein_ids(record, protein_id_db):
    changed = 0

    for feature in record.features:
        if feature.type != "CDS":
            continue

        locus_tag = first_qualifier(feature, "locus_tag")

        if locus_tag:
            protein_core = sanitize_for_protein_id(locus_tag)
        else:
            protein_core = sanitize_for_protein_id(
                f"{record.id}_{int(feature.location.start) + 1}_{int(feature.location.end)}"
            )

        protein_id = f"gnl|{protein_id_db}|{protein_core}"

        current = first_qualifier(feature, "protein_id")
        if current != protein_id:
            feature.qualifiers["protein_id"] = [protein_id]
            changed += 1
            logging.info(
                "Set /protein_id=%s for CDS at %s",
                protein_id,
                feature.location,
            )

    return changed


def ensure_gene_locus_tags(record):
    associated_types = {"CDS", "tRNA", "rRNA", "ncRNA", "tmRNA", "misc_RNA"}
    prefix = infer_locus_tag_prefix(record)
    used_tags = existing_locus_tags(record)
    added = 0

    for gene_feature in record.features:
        if gene_feature.type != "gene":
            continue

        existing = first_qualifier(gene_feature, "locus_tag")
        if existing and existing.strip():
            continue

        copied_locus_tag = None

        for feature in record.features:
            if feature.type not in associated_types:
                continue

            if same_location(gene_feature, feature):
                candidate = first_qualifier(feature, "locus_tag")
                if candidate:
                    copied_locus_tag = candidate
                    break

        if copied_locus_tag:
            gene_feature.qualifiers["locus_tag"] = [copied_locus_tag]
            used_tags.add(copied_locus_tag)
            added += 1
            logging.info(
                "Copied /locus_tag=%s to gene at %s",
                copied_locus_tag,
                gene_feature.location,
            )
        else:
            generated = next_generated_locus_tag(record, prefix, used_tags)
            gene_feature.qualifiers["locus_tag"] = [generated]
            added += 1
            logging.info(
                "Generated /locus_tag=%s for gene at %s",
                generated,
                gene_feature.location,
            )

    return added


def clean_suspect_product_names(record):
    changed = 0

    for feature in record.features:
        if feature.type != "CDS":
            continue

        product = first_qualifier(feature, "product")
        if not product:
            continue

        if product in PRODUCT_REPLACEMENTS:
            new_product = PRODUCT_REPLACEMENTS[product]
            feature.qualifiers["product"] = [new_product]
            changed += 1
            logging.info(
                "Changed suspect product name from %s to %s at %s",
                product,
                new_product,
                feature.location,
            )

    return changed


def remove_empty_gene_features(record):
    kept_features = []
    removed = 0

    for feature in record.features:
        if feature.type != "gene":
            kept_features.append(feature)
            continue

        has_gene = bool(feature.qualifiers.get("gene"))
        has_locus_tag = bool(feature.qualifiers.get("locus_tag"))
        has_locus = bool(feature.qualifiers.get("locus"))

        if has_gene or has_locus_tag or has_locus:
            kept_features.append(feature)
        else:
            removed += 1
            logging.info("Removed empty gene feature at %s", feature.location)

    record.features = kept_features
    return removed


def remove_superfluous_gene_features(record):
    associated_types = {"CDS", "tRNA", "rRNA", "ncRNA", "tmRNA", "misc_RNA"}

    associated_locations = {
        str(feature.location)
        for feature in record.features
        if feature.type in associated_types
    }

    kept_features = []
    removed = 0

    for feature in record.features:
        if feature.type != "gene":
            kept_features.append(feature)
            continue

        if "pseudo" in feature.qualifiers:
            kept_features.append(feature)
            continue

        if str(feature.location) in associated_locations:
            kept_features.append(feature)
        else:
            removed += 1
            logging.info(
                "Removed superfluous gene feature at %s; locus_tag=%s; gene=%s",
                feature.location,
                first_qualifier(feature, "locus_tag"),
                first_qualifier(feature, "gene"),
            )

    record.features = kept_features
    return removed


def remove_inversion_end_misc_feature(record):
    kept_features = []
    removed = 0

    for feature in record.features:
        if (
            feature.type == "misc_feature"
            and first_qualifier(feature, "note") == "Inversion_1_end"
        ):
            removed += 1
            logging.info(
                "Removed misc_feature Inversion_1_end at %s",
                feature.location,
            )
            continue

        kept_features.append(feature)

    record.features = kept_features
    return removed


def remove_illegal_ec_numbers(record):
    removed = 0
    legal_ec_pattern = re.compile(r"^(\d+|-)\.(\d+|-)\.(\d+|-)\.(\d+|-)$")

    for feature in record.features:
        if "EC_number" not in feature.qualifiers:
            continue

        valid_values = []

        for value in feature.qualifiers["EC_number"]:
            value = value.strip()

            if value == "2.4.99.28":
                removed += 1
                logging.info(
                    "Removed illegal /EC_number=%s at %s",
                    value,
                    feature.location,
                )
                continue

            if legal_ec_pattern.match(value):
                valid_values.append(value)
            else:
                removed += 1
                logging.info(
                    "Removed malformed /EC_number=%s at %s",
                    value,
                    feature.location,
                )

        if valid_values:
            feature.qualifiers["EC_number"] = valid_values
        else:
            del feature.qualifiers["EC_number"]

    return removed


def remove_cds_gene_equal_to_locus_tag(record):
    removed = 0

    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene_value = first_qualifier(feature, "gene")
        locus_tag_value = first_qualifier(feature, "locus_tag")

        if gene_value and locus_tag_value and gene_value == locus_tag_value:
            del feature.qualifiers["gene"]
            removed += 1
            logging.info(
                "Removed CDS /gene= identical to /locus_tag= at %s",
                feature.location,
            )

    return removed


def force_cds_products_in_flatfile(path):
    with open(path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    output = []
    cds_blocks_fixed = 0
    i = 0

    while i < len(lines):
        line = lines[i]

        if re.match(r"^     CDS\s+", line):
            block = [line]
            i += 1

            while i < len(lines):
                next_line = lines[i]

                if re.match(r"^     \S", next_line):
                    break

                block.append(next_line)
                i += 1

            block_text = "".join(block)
            has_product = re.search(r'/product="[^"]*\S[^"]*"', block_text)

            if not has_product:
                insert_at = len(block)

                for idx, block_line in enumerate(block):
                    if "/translation=" in block_line:
                        insert_at = idx
                        break

                block.insert(
                    insert_at,
                    '                     /product="hypothetical protein"\n',
                )

                cds_blocks_fixed += 1

            output.extend(block)
            continue

        output.append(line)
        i += 1

    with open(path, "w", encoding="utf-8") as handle:
        handle.writelines(output)

    return cds_blocks_fixed


def force_single_line_anticodons(path):
    with open(path, "r", encoding="utf-8") as handle:
        text = handle.read()

    pattern = re.compile(r'(/anticodon="(?:[^"]|\n\s*)*")', flags=re.MULTILINE)
    replacements = 0

    def flatten(match):
        nonlocal replacements
        replacements += 1
        return re.sub(r"\n\s+", "", match.group(1))

    text = pattern.sub(flatten, text)

    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)

    return replacements


def process_file(input_path, output_path, protein_id_db):
    logging.info("Reading GenBank file: %s", input_path)

    records = list(SeqIO.parse(input_path, "genbank"))

    logging.info("Loaded %d record(s)", len(records))

    total_added_genes = 0
    total_removed_pseudo_cds = 0
    total_removed_joined_or_slippage_cds = 0
    total_removed_overlapping_trna = 0
    total_removed_empty_companion_trna_genes = 0
    total_removed_dbxref = 0
    total_added_products = 0
    total_removed_notes = 0
    total_added_cds_genes = 0
    total_set_protein_ids = 0
    total_added_gene_locus_tags = 0
    total_cleaned_suspect_products = 0
    total_removed_empty_genes = 0
    total_removed_bad_ec_numbers = 0
    total_removed_cds_gene_equal_locus_tag = 0
    total_removed_superfluous_genes = 0
    total_removed_inversion_end_misc_features = 0

    for index, record in enumerate(records, start=1):
        logging.info("Processing record %d: %s", index, record.id)

        total_added_genes += add_missing_gene_features(record)
        total_removed_pseudo_cds += remove_pseudo_cds(record)

        total_removed_joined_or_slippage_cds += (
            remove_all_joined_or_ribosomal_slippage_cds(record)
        )

        removed_trna, removed_companion_genes = remove_trna_features_overlapping_cds(
            record
        )
        total_removed_overlapping_trna += removed_trna
        total_removed_empty_companion_trna_genes += removed_companion_genes

        (
            removed_dbxref,
            added_products,
            removed_notes,
            added_cds_genes,
        ) = clean_qualifiers(record)

        total_removed_dbxref += removed_dbxref
        total_added_products += added_products
        total_removed_notes += removed_notes
        total_added_cds_genes += added_cds_genes

        total_set_protein_ids += ensure_cds_protein_ids(record, protein_id_db)
        total_cleaned_suspect_products += clean_suspect_product_names(record)
        total_removed_empty_genes += remove_empty_gene_features(record)
        total_removed_bad_ec_numbers += remove_illegal_ec_numbers(record)

        total_removed_cds_gene_equal_locus_tag += (
            remove_cds_gene_equal_to_locus_tag(record)
        )

        total_removed_superfluous_genes += remove_superfluous_gene_features(record)
        total_added_gene_locus_tags += ensure_gene_locus_tags(record)

        total_removed_inversion_end_misc_features += (
            remove_inversion_end_misc_feature(record)
        )

    logging.info("Writing corrected file: %s", output_path)

    with open(output_path, "w", encoding="utf-8") as handle:
        SeqIO.write(records, handle, "genbank")

    cds_products_fixed_postwrite = force_cds_products_in_flatfile(output_path)
    anticodon_fixed = force_single_line_anticodons(output_path)

    logging.info("Completed processing")
    logging.info("Summary:")
    logging.info("  Added gene features: %d", total_added_genes)
    logging.info("  Removed pseudo CDS annotations: %d", total_removed_pseudo_cds)
    logging.info(
        "  Removed joined/ribosomal-slippage CDS annotations: %d",
        total_removed_joined_or_slippage_cds,
    )
    logging.info(
        "  Removed tRNA features overlapping CDS: %d",
        total_removed_overlapping_trna,
    )
    logging.info(
        "  Removed empty companion genes for removed tRNAs: %d",
        total_removed_empty_companion_trna_genes,
    )
    logging.info("  Removed db_xref tags: %d", total_removed_dbxref)
    logging.info("  Added CDS product tags: %d", total_added_products)
    logging.info(
        "  CDS product tags added during post-write flatfile repair: %d",
        cds_products_fixed_postwrite,
    )
    logging.info("  Removed note tags: %d", total_removed_notes)
    logging.info("  Added CDS gene tags: %d", total_added_cds_genes)
    logging.info(
        "  Set CDS protein_id tags to gnl format: %d",
        total_set_protein_ids,
    )
    logging.info("  Added gene locus_tag qualifiers: %d", total_added_gene_locus_tags)
    logging.info("  Cleaned suspect product names: %d", total_cleaned_suspect_products)
    logging.info("  Removed empty gene features: %d", total_removed_empty_genes)
    logging.info(
        "  Removed illegal or malformed EC_number tags: %d",
        total_removed_bad_ec_numbers,
    )
    logging.info(
        "  Removed CDS /gene= values identical to /locus_tag=: %d",
        total_removed_cds_gene_equal_locus_tag,
    )
    logging.info(
        "  Removed superfluous gene features: %d",
        total_removed_superfluous_genes,
    )
    logging.info(
        "  Removed Inversion_1_end misc_features: %d",
        total_removed_inversion_end_misc_features,
    )
    logging.info("  Flattened anticodon tags: %d", anticodon_fixed)


def main():
    parser = argparse.ArgumentParser(
        description="Final corrections for GenBank flat files before submission via Geneious."
    )
    parser.add_argument("input", help="Input GenBank flat file")
    parser.add_argument("output", help="Output GenBank flat file")
    parser.add_argument(
        "--protein-id-db",
        default="local",
        help="Database prefix used when setting /protein_id= tags. Default: local",
    )

    args = parser.parse_args()

    setup_logging()

    try:
        process_file(args.input, args.output, args.protein_id_db)
    except Exception as exc:
        logging.exception("Fatal error: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
