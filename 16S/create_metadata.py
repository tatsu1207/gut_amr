#!/usr/bin/env python3
"""Create metadata file from ASV table sample names for the KJM obesity+antibiotics study."""

import re
import sys

# Read sample names from ASV table header
with open("/home/unnot/vault/kjm/miseq/dada2_output/asv_table.tsv") as f:
    header = f.readline().strip().split("\t")
samples = header[2:]  # skip ASV_ID and sequence columns

# Group code mapping
group_map = {
    "CTL": ("CTL", "ND", "None"),
    "HFD": ("HFD", "HFD", "None"),
    "AC": ("AMX_CTL", "ND", "Amoxicillin"),
    "AH": ("AMX_HFD", "HFD", "Amoxicillin"),
    "CC": ("CIP_CTL", "ND", "Ciprofloxacin"),
    "CH": ("CIP_HFD", "HFD", "Ciprofloxacin"),
    "EC": ("ERY_CTL", "ND", "Erythromycin"),
    "EH": ("ERY_HFD", "HFD", "Erythromycin"),
    "AEar": ("AMX_Early", "ND", "Amoxicillin"),
    "CEar": ("CIP_Early", "ND", "Ciprofloxacin"),
    "CSL": ("CTL", "ND", "None"),  # CTL at SD14
}

# Timepoint mapping to experiment day
timepoint_day = {
    "T0h": 0, "T6h": 0.25, "T24h": 1, "T72h": 3,
    "S72h": 6, "SD7": 10, "SD10": 13, "SD14": 17,
}

# Timepoint phase
def get_phase(tp):
    if tp in ("T0h",):
        return "Baseline"
    elif tp in ("T6h", "T24h", "T72h"):
        return "Treatment"
    else:
        return "Recovery"

rows = []
for s in samples:
    experiment = "other"
    group = diet = antibiotic = timepoint = mouse_id = "NA"
    day = "NA"
    phase = "NA"

    # ARG samples - main study
    if s.startswith("ARG_"):
        rest = s[4:]  # remove ARG_ prefix

        # Pattern 1: ARG_{group}{mouse}_{timepoint} (e.g., ARG_AC1_S72h, ARG_CTL5_T0h)
        # Pattern 2: ARG_{timepoint}_{tissue}_{group}{mouse} (e.g., ARG_D0_C_CTL1)
        # Pattern 3: ARG_{timepoint}_{group}{mouse} (e.g., ARG_D0_CTL1, ARG_S72h_CTL1)
        # Pattern 4: ARG_{timepoint}_{tissue}_{group}{mouse} (e.g., ARG_T6h_C_EAR1)

        # Try D0 prefix patterns first
        m = re.match(r"D0_([A-Z])_([A-Z]+)(\d+)$", rest)
        if m:
            tissue, grp, mid = m.groups()
            experiment = "ARG_D0"
            mouse_id = mid
            timepoint = "T0h"
            day = "0"
            phase = "Baseline"
            tissue_map = {"A": "adipose", "C": "cecum", "N": "nasal"}
            if grp == "CTL":
                group, diet, antibiotic = "CTL", "ND", "None"
            elif grp == "HFD":
                group, diet, antibiotic = "HFD", "HFD", "None"
            elif grp == "EAR":
                if tissue == "C":
                    group, diet, antibiotic = "CIP_Early", "ND", "Ciprofloxacin"
                elif tissue == "N":
                    group, diet, antibiotic = "AMX_Early", "ND", "Amoxicillin"
                else:
                    group, diet, antibiotic = "Early_" + tissue, "ND", "None"
            rows.append([s, experiment, group, diet, antibiotic, timepoint, str(day), phase, mouse_id])
            continue

        # D0 without tissue
        m = re.match(r"D0_(CTL|HFD)(\d+)$", rest)
        if m:
            grp, mid = m.groups()
            experiment = "ARG_D0"
            mouse_id = mid
            timepoint = "T0h"
            day = "0"
            phase = "Baseline"
            if grp == "CTL":
                group, diet, antibiotic = "CTL", "ND", "None"
            else:
                group, diet, antibiotic = "HFD", "HFD", "None"
            rows.append([s, experiment, group, diet, antibiotic, timepoint, str(day), phase, mouse_id])
            continue

        # Timepoint prefix with N/C + EAR (C_EAR = CIP-Early, N_EAR = AMX-Early)
        # Must come before generic tissue pattern
        m = re.match(r"(T6h|T24h|T72h|S72h)_([NC])_EAR(\d+)$", rest)
        if m:
            tp, tissue, mid = m.groups()
            experiment = "ARG"
            mouse_id = mid
            timepoint = tp
            day = str(timepoint_day.get(tp, "NA"))
            phase = get_phase(tp)
            if tissue == "C":
                group, diet, antibiotic = "CIP_Early", "ND", "Ciprofloxacin"
            else:
                group, diet, antibiotic = "AMX_Early", "ND", "Amoxicillin"
            rows.append([s, experiment, group, diet, antibiotic, timepoint, str(day), phase, mouse_id])
            continue

        # Timepoint prefix without tissue (e.g., T6h_CTL1, S72h_HFD3)
        m = re.match(r"(T6h|T24h|T72h|S72h)_(CTL|HFD)(\d+)$", rest)
        if m:
            tp, grp, mid = m.groups()
            experiment = "ARG"
            mouse_id = mid
            timepoint = tp
            day = str(timepoint_day.get(tp, "NA"))
            phase = get_phase(tp)
            if grp == "CTL":
                group, diet, antibiotic = "CTL", "ND", "None"
            else:
                group, diet, antibiotic = "HFD", "HFD", "None"
            rows.append([s, experiment, group, diet, antibiotic, timepoint, str(day), phase, mouse_id])
            continue

        # Standard: group+mouse_timepoint (e.g., AC1_S72h, CTL5_T0h, HFD2_SD14)
        m = re.match(r"(CSL|AEar|CEar|CTL|HFD|AC|AH|CC|CH|EC|EH)(\d+)_(T0h|T6h|T24h|T72h|S72h|SD7|SD10|SD14)$", rest)
        if m:
            grp, mid, tp = m.groups()
            experiment = "ARG"
            mouse_id = mid
            timepoint = tp
            day = str(timepoint_day.get(tp, "NA"))
            phase = get_phase(tp)
            if grp in group_map:
                group, diet, antibiotic = group_map[grp]
            rows.append([s, experiment, group, diet, antibiotic, timepoint, str(day), phase, mouse_id])
            continue

        # If nothing matched
        rows.append([s, "ARG_unmatched", "NA", "NA", "NA", "NA", "NA", "NA", "NA"])
        continue

    # Non-ARG samples (B_ob, D0_, HD/ND 2nd experiment)
    rows.append([s, "other", "NA", "NA", "NA", "NA", "NA", "NA", "NA"])

# Write output
outfile = "/home/unnot/vault/kjm/miseq/dada2_output/metadata.tsv"
with open(outfile, "w") as f:
    f.write("sample_id\texperiment\tgroup\tdiet\tantibiotic\ttimepoint\tday\tphase\tmouse_id\n")
    for row in rows:
        f.write("\t".join(row) + "\n")

# Summary
from collections import Counter
exp_counts = Counter(r[1] for r in rows)
grp_counts = Counter(r[2] for r in rows if r[1] == "ARG")
tp_counts = Counter(r[5] for r in rows if r[1] == "ARG")
print(f"Total samples: {len(rows)}")
print(f"\nBy experiment: {dict(exp_counts)}")
print(f"\nARG group counts: {dict(grp_counts)}")
print(f"\nARG timepoint counts: {dict(tp_counts)}")
unmatched = [r[0] for r in rows if "unmatched" in r[1]]
if unmatched:
    print(f"\nUnmatched: {unmatched[:10]}")
