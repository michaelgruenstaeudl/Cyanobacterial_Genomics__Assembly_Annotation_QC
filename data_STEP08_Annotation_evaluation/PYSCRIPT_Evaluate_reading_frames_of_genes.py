#!/usr/bin/env python3
"""
Evaluates if the reading frames of the genes of a bacterial genome 
(available in GenBank flatfile format) are intact.

Typical usage:
    python PYSCRIPT_Evaluate_reading_frames_of_genes.py genome.gb -o genome_summary

Outputs:
    <prefix>.replicons.tsv   per-record/per-replicon metrics
    <prefix>.genome.tsv      whole-genome aggregated metrics

Notes
-----
- Works with one or many GenBank records in the same file (chromosome + plasmids).
- Uses annotated features already present in the GenBank file.
- "Complete vs partial rRNA" is inferred heuristically from GenBank feature
  locations/qualifiers:
    * partial if location is fuzzy (< or >) or qualifier "partial" is present
    * otherwise counted as complete
- "Intergenic spacers" are computed from merged annotated intervals formed from
  gene/CDS/RNA/pseudogene features, so the result is annotation-centric.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import statistics
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


# ----------------------------- helpers --------------------------------- #

def safe_mean(values: Sequence[float]) -> float:
    return float(statistics.mean(values)) if values else float("nan")


def safe_median(values: Sequence[float]) -> float:
    return float(statistics.median(values)) if values else float("nan")


def gc_percent(seq: str) -> float:
    seq = seq.upper()
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    atgc = a + c + g + t
    if atgc == 0:
        return float("nan")
    return 100.0 * (g + c) / atgc


def format_value(v):
    if isinstance(v, float):
        if math.isnan(v):
            return "nan"
        return f"{v:.1f}"
    return v


def get_feature_span(feature: SeqFeature) -> Tuple[int, int]:
    """
    Return zero-based half-open [start, end) span across the full feature.
    For joined/compound features, this is the min start to max end span.
    """
    loc = feature.location
    if isinstance(loc, CompoundLocation):
        starts = [int(part.start) for part in loc.parts]
        ends = [int(part.end) for part in loc.parts]
        return min(starts), max(ends)
    return int(loc.start), int(loc.end)


def get_feature_length(feature: SeqFeature) -> int:
    try:
        return int(len(feature.location))
    except Exception:
        start, end = get_feature_span(feature)
        return max(0, end - start)


def has_partial_location(feature: SeqFeature) -> bool:
    """
    Heuristic detection of partial/fuzzy locations.
    In Biopython, '<' and '>' usually become BeforePosition / AfterPosition.
    """
    loc = feature.location

    def _pos_is_partial(pos) -> bool:
        cls = pos.__class__.__name__
        return cls in {"BeforePosition", "AfterPosition", "WithinPosition", "OneOfPosition"}

    if isinstance(loc, CompoundLocation):
        for part in loc.parts:
            if _pos_is_partial(part.start) or _pos_is_partial(part.end):
                return True
        return False

    return _pos_is_partial(loc.start) or _pos_is_partial(loc.end)


def qualifier_contains(feature: SeqFeature, key: str, pattern: str) -> bool:
    vals = feature.qualifiers.get(key, [])
    regex = re.compile(pattern, flags=re.IGNORECASE)
    return any(regex.search(v) for v in vals)


def is_pseudogene_feature(feature: SeqFeature) -> bool:
    if feature.type == "pseudogene":
        return True
    if "pseudo" in feature.qualifiers:
        return True
    if "pseudogene" in feature.qualifiers:
        return True
    return False


def infer_rrna_subtype(feature: SeqFeature) -> str:
    """
    Return one of: 5S, 16S, 23S, other, unknown
    """
    texts = []
    for key in ("product", "gene", "note", "standard_name"):
        texts.extend(feature.qualifiers.get(key, []))
    blob = " ".join(texts).lower()

    if "16s" in blob or "small subunit ribosomal rna" in blob or "ssu ribosomal rna" in blob:
        return "16S"
    if "23s" in blob or "large subunit ribosomal rna" in blob or "lsu ribosomal rna" in blob:
        return "23S"
    if "5s" in blob:
        return "5S"
    if "rrna" in blob or "ribosomal rna" in blob:
        return "other"
    return "unknown"


def is_nc_rna(feature: SeqFeature) -> bool:
    if feature.type in {"ncRNA", "misc_RNA"}:
        return True
    if feature.type == "RNA":
        return True
    return False


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    ints = sorted((s, e) for s, e in intervals if e > s)
    if not ints:
        return []
    merged = [ints[0]]
    for s, e in ints[1:]:
        last_s, last_e = merged[-1]
        if s <= last_e:
            merged[-1] = (last_s, max(last_e, e))
        else:
            merged.append((s, e))
    return merged


def interval_gaps_and_overlaps(intervals: Iterable[Tuple[int, int]]) -> Tuple[List[int], List[int]]:
    """
    For sorted merged or unmerged intervals, compute positive gaps and positive overlaps
    between adjacent intervals after sorting by start,end.
    """
    ints = sorted((s, e) for s, e in intervals if e > s)
    gaps: List[int] = []
    overlaps: List[int] = []
    if len(ints) < 2:
        return gaps, overlaps

    prev_s, prev_e = ints[0]
    for s, e in ints[1:]:
        diff = s - prev_e
        if diff > 0:
            gaps.append(diff)
        elif diff < 0:
            overlaps.append(-diff)
        if e > prev_e:
            prev_s, prev_e = s, e
        else:
            prev_s, prev_e = prev_s, prev_e
    return gaps, overlaps


@dataclass
class RepliconSummary:
    replicon_id: str
    replicon_name: str
    record_type: str
    topology: str
    length_bp: int
    gc_percent: float

    genes_total: int
    cds_total: int
    cds_with_translation: int
    pseudogenes_total: int

    trna_total: int
    rrna_total: int
    rrna_5s: int
    rrna_16s: int
    rrna_23s: int
    rrna_other: int
    rrna_complete: int
    rrna_partial: int

    ncrna_total: int
    tmrna_total: int
    other_rna_total: int

    avg_gene_length_bp: float
    median_gene_length_bp: float
    avg_cds_length_bp: float
    median_cds_length_bp: float

    annotated_bases_bp: int
    coding_bases_bp: int
    coding_density_percent: float

    intergenic_spacers_n: int
    intergenic_spacers_total_bp: int
    intergenic_spacers_mean_bp: float
    intergenic_spacers_median_bp: float
    intergenic_spacers_min_bp: float
    intergenic_spacers_max_bp: float

    annotated_overlaps_n: int
    annotated_overlaps_total_bp: int
    annotated_overlaps_mean_bp: float


# -------------------------- core summarization -------------------------- #

def summarize_record(record: SeqRecord) -> RepliconSummary:
    seq = str(record.seq)
    length_bp = len(seq)
    gc = gc_percent(seq)

    genes_total = 0
    cds_total = 0
    cds_with_translation = 0
    pseudogenes_total = 0

    trna_total = 0
    rrna_total = 0
    rrna_5s = 0
    rrna_16s = 0
    rrna_23s = 0
    rrna_other = 0
    rrna_complete = 0
    rrna_partial = 0

    ncrna_total = 0
    tmrna_total = 0
    other_rna_total = 0

    gene_lengths: List[int] = []
    cds_lengths: List[int] = []

    annotation_intervals: List[Tuple[int, int]] = []
    coding_intervals: List[Tuple[int, int]] = []

    for feat in record.features:
        ftype = feat.type

        if is_pseudogene_feature(feat):
            pseudogenes_total += 1

        if ftype == "gene":
            genes_total += 1
            gene_lengths.append(get_feature_length(feat))
            annotation_intervals.append(get_feature_span(feat))

        elif ftype == "CDS":
            cds_total += 1
            cds_len = get_feature_length(feat)
            cds_lengths.append(cds_len)
            annotation_intervals.append(get_feature_span(feat))
            coding_intervals.append(get_feature_span(feat))
            if "translation" in feat.qualifiers:
                cds_with_translation += 1

        elif ftype == "tRNA":
            trna_total += 1
            annotation_intervals.append(get_feature_span(feat))

        elif ftype == "rRNA":
            rrna_total += 1
            annotation_intervals.append(get_feature_span(feat))
            subtype = infer_rrna_subtype(feat)
            if subtype == "5S":
                rrna_5s += 1
            elif subtype == "16S":
                rrna_16s += 1
            elif subtype == "23S":
                rrna_23s += 1
            else:
                rrna_other += 1

            partial = (
                has_partial_location(feat)
                or qualifier_contains(feat, "note", r"\bpartial\b")
                or qualifier_contains(feat, "product", r"\bpartial\b")
            )
            if partial:
                rrna_partial += 1
            else:
                rrna_complete += 1

        elif ftype == "tmRNA":
            tmrna_total += 1
            annotation_intervals.append(get_feature_span(feat))

        elif is_nc_rna(feat):
            ncrna_total += 1
            annotation_intervals.append(get_feature_span(feat))

        elif ftype.endswith("RNA"):
            other_rna_total += 1
            annotation_intervals.append(get_feature_span(feat))

    merged_annotation = merge_intervals(annotation_intervals)
    merged_coding = merge_intervals(coding_intervals)

    annotated_bases_bp = sum(e - s for s, e in merged_annotation)
    coding_bases_bp = sum(e - s for s, e in merged_coding)
    coding_density_percent = (100.0 * coding_bases_bp / length_bp) if length_bp else float("nan")

    gaps, overlaps = interval_gaps_and_overlaps(merged_annotation)

    _, raw_overlaps = interval_gaps_and_overlaps(annotation_intervals)

    source = None
    for feat in record.features:
        if feat.type == "source":
            source = feat
            break

    topology = "unknown"
    if source is not None:
        topology_vals = source.qualifiers.get("topology", [])
        if topology_vals:
            topology = topology_vals[0]

    record_type = "replicon"
    descr = (record.description or "").lower()
    if "plasmid" in descr:
        record_type = "plasmid"
    else:
        source_plasmid = source.qualifiers.get("plasmid", []) if source else []
        if source_plasmid:
            record_type = "plasmid"
        else:
            record_type = "chromosome_or_other"

    return RepliconSummary(
        replicon_id=record.id,
        replicon_name=record.name,
        record_type=record_type,
        topology=topology,
        length_bp=length_bp,
        gc_percent=gc,

        genes_total=genes_total,
        cds_total=cds_total,
        cds_with_translation=cds_with_translation,
        pseudogenes_total=pseudogenes_total,

        trna_total=trna_total,
        rrna_total=rrna_total,
        rrna_5s=rrna_5s,
        rrna_16s=rrna_16s,
        rrna_23s=rrna_23s,
        rrna_other=rrna_other,
        rrna_complete=rrna_complete,
        rrna_partial=rrna_partial,

        ncrna_total=ncrna_total,
        tmrna_total=tmrna_total,
        other_rna_total=other_rna_total,

        avg_gene_length_bp=safe_mean(gene_lengths),
        median_gene_length_bp=safe_median(gene_lengths),
        avg_cds_length_bp=safe_mean(cds_lengths),
        median_cds_length_bp=safe_median(cds_lengths),

        annotated_bases_bp=annotated_bases_bp,
        coding_bases_bp=coding_bases_bp,
        coding_density_percent=coding_density_percent,

        intergenic_spacers_n=len(gaps),
        intergenic_spacers_total_bp=sum(gaps),
        intergenic_spacers_mean_bp=safe_mean(gaps),
        intergenic_spacers_median_bp=safe_median(gaps),
        intergenic_spacers_min_bp=min(gaps) if gaps else float("nan"),
        intergenic_spacers_max_bp=max(gaps) if gaps else float("nan"),

        annotated_overlaps_n=len(raw_overlaps),
        annotated_overlaps_total_bp=sum(raw_overlaps),
        annotated_overlaps_mean_bp=safe_mean(raw_overlaps),
    )


def aggregate_genome(summaries: List[RepliconSummary]) -> Dict[str, float]:
    total_length = sum(s.length_bp for s in summaries)

    if total_length:
        genome_gc = sum(
            s.gc_percent * s.length_bp
            for s in summaries
            if not math.isnan(s.gc_percent)
        ) / total_length
    else:
        genome_gc = float("nan")

    out = {
        "replicons_total": len(summaries),
        "chromosome_or_other_replicons": sum(1 for s in summaries if s.record_type == "chromosome_or_other"),
        "plasmids_total": sum(1 for s in summaries if s.record_type == "plasmid"),
        "length_bp_total": total_length,
        "gc_percent_weighted": genome_gc,

        "genes_total": sum(s.genes_total for s in summaries),
        "cds_total": sum(s.cds_total for s in summaries),
        "cds_with_translation_total": sum(s.cds_with_translation for s in summaries),
        "pseudogenes_total": sum(s.pseudogenes_total for s in summaries),

        "trna_total": sum(s.trna_total for s in summaries),
        "rrna_total": sum(s.rrna_total for s in summaries),
        "rrna_5s_total": sum(s.rrna_5s for s in summaries),
        "rrna_16s_total": sum(s.rrna_16s for s in summaries),
        "rrna_23s_total": sum(s.rrna_23s for s in summaries),
        "rrna_other_total": sum(s.rrna_other for s in summaries),
        "rrna_complete_total": sum(s.rrna_complete for s in summaries),
        "rrna_partial_total": sum(s.rrna_partial for s in summaries),

        "ncrna_total": sum(s.ncrna_total for s in summaries),
        "tmrna_total": sum(s.tmrna_total for s in summaries),
        "other_rna_total": sum(s.other_rna_total for s in summaries),

        "annotated_bases_bp_total": sum(s.annotated_bases_bp for s in summaries),
        "coding_bases_bp_total": sum(s.coding_bases_bp for s in summaries),

        "intergenic_spacers_n_total": sum(s.intergenic_spacers_n for s in summaries),
        "intergenic_spacers_bp_total": sum(s.intergenic_spacers_total_bp for s in summaries),

        "annotated_overlaps_n_total": sum(s.annotated_overlaps_n for s in summaries),
        "annotated_overlaps_bp_total": sum(s.annotated_overlaps_total_bp for s in summaries),
    }

    if total_length:
        out["coding_density_percent"] = 100.0 * out["coding_bases_bp_total"] / total_length
    else:
        out["coding_density_percent"] = float("nan")

    def weighted_avg(attr: str, weight_attr: str = "length_bp") -> float:
        nums = []
        dens = []
        for s in summaries:
            val = getattr(s, attr)
            w = getattr(s, weight_attr)
            if not math.isnan(val):
                nums.append(val * w)
                dens.append(w)
        return sum(nums) / sum(dens) if dens and sum(dens) else float("nan")

    out["avg_gene_length_bp_weighted"] = weighted_avg("avg_gene_length_bp")
    out["avg_cds_length_bp_weighted"] = weighted_avg("avg_cds_length_bp")

    out["intergenic_spacers_mean_bp"] = (
        out["intergenic_spacers_bp_total"] / out["intergenic_spacers_n_total"]
        if out["intergenic_spacers_n_total"] else float("nan")
    )
    out["annotated_overlaps_mean_bp"] = (
        out["annotated_overlaps_bp_total"] / out["annotated_overlaps_n_total"]
        if out["annotated_overlaps_n_total"] else float("nan")
    )

    out["genes_per_mbp"] = (out["genes_total"] / total_length * 1e6) if total_length else float("nan")
    out["cds_per_mbp"] = (out["cds_total"] / total_length * 1e6) if total_length else float("nan")
    out["trna_per_mbp"] = (out["trna_total"] / total_length * 1e6) if total_length else float("nan")
    out["rrna_operon_proxy_min"] = min(
        out["rrna_5s_total"], out["rrna_16s_total"], out["rrna_23s_total"]
    )

    return out


# ------------------------------- output -------------------------------- #

def write_replicon_tsv(summaries: List[RepliconSummary], path: str) -> None:
    if not summaries:
        raise ValueError("No summaries to write")

    fieldnames = list(summaries[0].__dataclass_fields__.keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for s in summaries:
            formatted = {k: format_value(v) for k, v in s.__dict__.items()}
            writer.writerow(formatted)


def write_genome_tsv(genome_summary: Dict[str, float], path: str) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in genome_summary.items():
            writer.writerow([key, format_value(value)])


def print_human_readable(genome_summary: Dict[str, float]) -> None:
    print("Genome summary")
    print("--------------")
    for key, value in genome_summary.items():
        print(f"{key}: {format_value(value)}")


# ------------------------------- main ---------------------------------- #

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Evaluates if the reading frames of the genes of a bacterial genome are intact.")
    p.add_argument("genbank", help="Input GenBank flatfile (.gb, .gbk, .gbff)")
    p.add_argument(
        "-o", "--out-prefix",
        default=None,
        help="Output prefix (default: input filename without extension)"
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    in_path = args.genbank
    out_prefix = args.out_prefix
    if out_prefix is None:
        base = os.path.basename(in_path)
        out_prefix = re.sub(r"\.(gb|gbk|gbff|genbank)$", "", base, flags=re.IGNORECASE)

    records = list(SeqIO.parse(in_path, "genbank"))
    if not records:
        raise SystemExit(f"No GenBank records found in {in_path!r}")

    summaries = [summarize_record(r) for r in records]
    genome_summary = aggregate_genome(summaries)

    replicon_tsv = f"{out_prefix}.replicons.tsv"
    genome_tsv = f"{out_prefix}.genome.tsv"

    write_replicon_tsv(summaries, replicon_tsv)
    write_genome_tsv(genome_summary, genome_tsv)
    print_human_readable(genome_summary)

    print("\nWrote:")
    print(f"  {replicon_tsv}")
    print(f"  {genome_tsv}")


if __name__ == "__main__":
    main()
