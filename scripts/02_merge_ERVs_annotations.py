#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
import re
from bisect import bisect_left
from collections import defaultdict
from tqdm import tqdm

# ==============================
# BED parsing / formatting
# ==============================

def parse_bed_line(line):
    fields = line.strip().split('\t')
    if len(fields) < 6:
        raise ValueError(f"Expected 6 BED columns, got {len(fields)}: {line}")
    chrom, start, end, name, score, strand = fields[:6]

    # e.g. "HERVK9-int_pos_chr1_..." -> "HERVK9-int" -> "HERVK9"
    left = name.split('_pos_')[0]
    subfamily = re.sub(r'-?int$', '', left)

    return {
        'chrom': chrom,
        'start': int(start),
        'end': int(end),
        'name': name,
        'score': score,
        'strand': strand,
        'subfamily': subfamily  # normalized base, no "-int"
    }


def format_bed_line(entry):
    new_name = f"{entry['subfamily']}-int_pos_{entry['chrom']}_{entry['start']+1}_{entry['end']}_strand_{entry['strand']}"
    return f"{entry['chrom']}\t{entry['start']}\t{entry['end']}\t{new_name}\t{entry['score']}\t{entry['strand']}"

# ==============================
# RepeatMasker .out parsing
# (1-based inclusive -> 0-based half-open)
# ==============================

def parse_rmsk_out(rmsk_path):
    def _to_int(x): return int(re.sub(r'[()]', '', x))
    records = []
    with open(rmsk_path, 'r') as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.lstrip().startswith(("SW", "score", "There were")):
                continue
            parts = line.strip().split()
            if len(parts) < 14:
                continue
            try:
                chrom = parts[4]
                g_begin = int(parts[5])  # 1-based inclusive
                g_end   = int(parts[6])  # 1-based inclusive
                orient_token = parts[8]
                rep_name = parts[9]
                classfam = parts[10]

                strand = '-' if orient_token.upper().startswith('C') else '+'
                m1, m2, m3 = parts[11], parts[12], parts[13]

                if strand == '+':
                    rep_begin = _to_int(m1)
                    rep_end   = _to_int(m2)
                    rep_left  = _to_int(m3)
                else:  # minus
                    rep_left  = _to_int(m1)
                    rep_end   = _to_int(m2)
                    rep_begin = _to_int(m3)

                # to BED half-open
                g_start0 = g_begin - 1
                g_end0   = g_end

                records.append({
                    'chrom': chrom,
                    'g_start': g_start0,
                    'g_end':   g_end0,
                    'strand': strand,
                    'repName': rep_name,   # e.g., HERVK9-int
                    'classfam': classfam,
                    'repBegin': rep_begin, # NOTE: for '-' this is the "start-like" monotone coord
                    'repEnd': rep_end,
                    'repLeft': rep_left
                })
            except Exception:
                continue
    return records


# ==============================
# RMSK indices (internal-only)
# ==============================

def build_rmsk_index(rmsk_records):
    """
    Build internal-only indices (repName endswith '-int').
      - by_cs:  key (chrom, strand) -> list of internal hits
      - by_cst: key (chrom, strand, token) -> list of internal hits
    """
    by_cs  = defaultdict(list)
    by_cst = defaultdict(list)
    starts_cs  = {}
    starts_cst = {}

    for r in rmsk_records:
        rep = r['repName']
        if not rep.endswith('-int'):
            continue  # only internal models
        key_cs = (r['chrom'], r['strand'])
        by_cs[key_cs].append(r)

        tokens = re.split(r'[_\-]', rep)
        for t in {t for t in tokens if t}:
            key_cst = (r['chrom'], r['strand'], t)
            by_cst[key_cst].append(r)

    # Sort & cache starts
    for d, starts in ((by_cs, starts_cs), (by_cst, starts_cst)):
        for k, lst in d.items():
            lst.sort(key=lambda x: (x['g_start'], x['g_end']))
            starts[k] = [x['g_start'] for x in lst]

    return (by_cs, starts_cs), (by_cst, starts_cst)

# ==============================
# Queries (internal & same-subfamily only)
# ==============================

def _halfopen_overlap(a0, a1, b0, b1):
    return (a0 < b1) and (a1 > b0)

def query_rmsk_internal_for_subfamily(
    idx_cs, starts_cs, idx_cst, starts_cst,
    chrom, strand, subfamily, g_start, g_end, pad=0
):
    out = []
    subfam_base = re.sub(r'-?int$', '', subfamily)

    def scan_list(lst):
        for rec in lst:
            if rec['g_start'] >= g_end + pad:
                break
            if _halfopen_overlap(rec['g_start'], rec['g_end'], g_start - pad, g_end + pad):
                base = re.sub(r'-?int$', '', rec['repName'])
                tokens = re.split(r'[_\-]', base)
                if base == subfam_base or base.startswith(subfam_base) or (subfam_base in tokens):
                    out.append(rec)

    # token and base lookups (kept as you had them)
    tokens = [t for t in re.split(r'[_\-]', subfamily) if t]
    seen = set()
    for t in tokens:
        k = (chrom, strand, t)
        if k in idx_cst and k not in seen:
            lst = idx_cst[k]
            left = bisect_left(starts_cst[k], max(g_start - pad, 0))
            scan_list(lst[left:])
            seen.add(k)

    base_sub = subfam_base
    k2 = (chrom, strand, base_sub)
    if k2 in idx_cst and k2 not in seen:
        lst = idx_cst[k2]
        left = bisect_left(starts_cst[k2], max(g_start - pad, 0))
        scan_list(lst[left:])
        seen.add(k2)

    if not out:
        kcs = (chrom, strand)
        if kcs in idx_cs:
            lst = idx_cs[kcs]
            left = bisect_left(starts_cs[kcs], max(g_start - pad, 0))
            scan_list(lst[left:])

    return out


def pick_terminal_internal_hits(
    idx_cs, starts_cs, idx_cst, starts_cst,
    chrom, strand, subfamily, start, end, pad=0
):
    hits = query_rmsk_internal_for_subfamily(
        idx_cs, starts_cs, idx_cst, starts_cst,
        chrom, strand, subfamily, start, end, pad=pad
    )
    if not hits:
        return None, None
    hits.sort(key=lambda h: (h['g_start'], h['g_end']))
    leftmost  = hits[0]
    rightmost = max(hits, key=lambda h: h['g_end'])
    return leftmost, rightmost


def model_continuity_ok_boundary_same_subfam(
    rmsk_ix_pack, left_entry, right_entry, pad=0, debug_rows=None
):
    (idx_cs, starts_cs, idx_cst, starts_cst) = rmsk_ix_pack

    same_group = (
        left_entry['chrom'] == right_entry['chrom'] and
        left_entry['strand'] == right_entry['strand'] and
        left_entry['subfamily'] == right_entry['subfamily']
    )
    if not same_group:
        if debug_rows is not None:
            debug_rows.append({
                "chrom": left_entry['chrom'],
                "strand": left_entry['strand'],
                "subfamily": left_entry['subfamily'],
                "left_start": left_entry['start'], "left_end": left_entry['end'],
                "right_start": right_entry['start'], "right_end": right_entry['end'],
                "pad": pad,
                "n_hits_left": 0, "n_hits_right": 0,
                "L_repBegin": "", "L_repEnd": "", "L_repLeft": "",
                "R_repBegin": "", "R_repEnd": "", "R_repLeft": "",
                "decision": False, "reason": "different_group"
            })
        return False

    chrom  = left_entry['chrom']
    strand = left_entry['strand']
    subfam = left_entry['subfamily']

    hits_L = query_rmsk_internal_for_subfamily(
        idx_cs, starts_cs, idx_cst, starts_cst,
        chrom, strand, subfam, left_entry['start'], left_entry['end'], pad=pad
    )
    hits_R = query_rmsk_internal_for_subfamily(
        idx_cs, starts_cs, idx_cst, starts_cst,
        chrom, strand, subfam, right_entry['start'], right_entry['end'], pad=pad
    )
    nL, nR = len(hits_L), len(hits_R)

    if nL == 0 or nR == 0:
        if debug_rows is not None:
            debug_rows.append({
                "chrom": chrom, "strand": strand, "subfamily": subfam,
                "left_start": left_entry['start'], "left_end": left_entry['end'],
                "right_start": right_entry['start'], "right_end": right_entry['end'],
                "pad": pad,
                "n_hits_left": nL, "n_hits_right": nR,
                "L_repBegin": "", "L_repEnd": "", "L_repLeft": "",
                "R_repBegin": "", "R_repEnd": "", "R_repLeft": "",
                "decision": None, "reason": "no_internal_hits_for_subfam"
            })
        return None

    # choose terminal hits (within each side)
    L_term = max(hits_L, key=lambda h: h['g_end'])    # rightmost within left
    R_term = min(hits_R, key=lambda h: h['g_start'])  # leftmost within right

    if strand == '+':
        decision = (L_term['repLeft'] > R_term['repLeft'])
        reason = "repLeft_decrease_+"
    else:
        # On '-' we use the "start-like" monotone coord across the junction.
        # From parse_rmsk_out, repBegin is the monotone field for '-'.
        decision = (R_term['repBegin'] < L_term['repBegin'])
        reason = "repBegin_smaller_on_right_-"

    if debug_rows is not None:
        debug_rows.append({
            "chrom": chrom, "strand": strand, "subfamily": subfam,
            "left_start": left_entry['start'], "left_end": left_entry['end'],
            "right_start": right_entry['start'], "right_end": right_entry['end'],
            "pad": pad,
            "n_hits_left": nL, "n_hits_right": nR,
            "L_repBegin": L_term['repBegin'], "L_repEnd": L_term['repEnd'], "L_repLeft": L_term['repLeft'],
            "R_repBegin": R_term['repBegin'], "R_repEnd": R_term['repEnd'], "R_repLeft": R_term['repLeft'],
            "decision": decision, "reason": reason
        })

    return decision


# ==============================
# Merging (two-stage, same subfamily)
# ==============================

def merge_entries_with_guardrail_same_subfam(groups, allowed_gap,
                                             rmsk_ix_pack=None,
                                             require_rmsk=False, pad=0,
                                             debug_rows=None):
    merged_out = []
    for (chrom, strand, subfamily), entries in tqdm(list(groups.items()),
                                                    desc="Merging groups", unit="group"):
        if not entries:
            continue
        entries.sort(key=lambda x: x['start'])
        merged = [entries[0]]

        for current in entries[1:]:
            last = merged[-1]
            can_merge_basic = (
                current['start'] <= last['end'] + allowed_gap and
                current['chrom'] == last['chrom'] and
                current['strand'] == last['strand'] and
                current['subfamily'] == last['subfamily']
            )
            if not can_merge_basic:
                merged.append(current)
                continue

            if rmsk_ix_pack is None:
                last['end'] = max(last['end'], current['end'])
                continue

            cont = model_continuity_ok_boundary_same_subfam(
                rmsk_ix_pack, last, current, pad=pad, debug_rows=debug_rows
            )
            if cont is True:
                last['end'] = max(last['end'], current['end'])
            elif cont is False:
                merged.append(current)
            else:  # cont is None
                if require_rmsk:
                    merged.append(current)
                else:
                    last['end'] = max(last['end'], current['end'])

        merged_out.extend(merged)
    return merged_out

# ==============================
# Main
# ==============================

def main(args):
    # Load BED
    with open(args.input_bed, 'r') as f:
        entries = [parse_bed_line(line) for line in f if line.strip()]

    # Build RMSK indices if provided
    rmsk_ix_pack = None
    if args.rmsk_out:
        rmsk_records = parse_rmsk_out(args.rmsk_out)
        (idx_cs, starts_cs), (idx_cst, starts_cst) = build_rmsk_index(rmsk_records)
        rmsk_ix_pack = (idx_cs, starts_cs, idx_cst, starts_cst)

    # --- DEBUG collector before any merges ---
    debug_rows = []

    # Stage A
    groups_A = defaultdict(list)
    for e in entries:
        groups_A[(e['chrom'], e['strand'], e['subfamily'])].append(e)

    merged_A = merge_entries_with_guardrail_same_subfam(
        groups_A, allowed_gap=args.short_gap,
        rmsk_ix_pack=None, require_rmsk=False, pad=0
    )

    # Stage B
    groups_B = defaultdict(list)
    for e in merged_A:
        groups_B[(e['chrom'], e['strand'], e['subfamily'])].append(e)

    merged_B = merge_entries_with_guardrail_same_subfam(
        groups_B, allowed_gap=args.long_gap,
        rmsk_ix_pack=rmsk_ix_pack,
        require_rmsk=args.require_rmsk,
        pad=args.pad,
        debug_rows=debug_rows,
    )

    # Sort & write BED
    merged_B.sort(key=lambda x: (x['chrom'], x['start'], x['end'], x['strand'], x['subfamily']))
    with open(args.output_bed, 'w') as f:
        for entry in merged_B:
            f.write(format_bed_line(entry) + '\n')

    # Write TSV if requested
    if args.debug_tsv:
        import csv
        with open(args.debug_tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "chrom","strand","subfamily",
                "left_start","left_end","right_start","right_end","pad",
                "n_hits_left","n_hits_right",
                "L_repBegin","L_repEnd","L_repLeft",
                "R_repBegin","R_repEnd","R_repLeft",
                "decision","reason"
            ], delimiter="\t")
            w.writeheader()
            for row in debug_rows:
                w.writerow(row)



if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Two-stage Telebuilder-like merge of ERV internal BED regions, "
                    "enforcing SAME subfamily and continuity using INTERNAL-only RMSK hits."
    )
    ap.add_argument("input_bed", help="Input BED (6 cols). Name should contain subfamily like 'HERVK13-1-int_...'.")
    ap.add_argument("output_bed", help="Output BED.")
    ap.add_argument("--short-gap", type=int, default=10,
                    help="Short unconditional merge distance [default: 10].")
    ap.add_argument("--long-gap", type=int, default=2500,
                    help="Long conditional merge distance with RMSK continuity [default: 2500].")
    ap.add_argument("--rmsk-out",
                    help="RepeatMasker .out file to enforce model-continuity guardrail (internal-only indices).")
    ap.add_argument("--require-rmsk", action="store_true",
                    help="Disallow long-gap merges lacking INTERNAL RMSK evidence for this subfamily.")
    ap.add_argument("--pad", type=int, default=10,
                    help="Padding (bp) for RMSK queries [default: 10].")
    ap.add_argument("--debug-tsv", help="Write per-boundary continuity decisions to this TSV.")

    args = ap.parse_args()

    main(args)
