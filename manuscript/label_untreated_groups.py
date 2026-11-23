"""Label the untreated groups explicitly in the published workbook, with a label pandas cannot coerce.

A blank `antibiotic` cell means "no antibiotic", but pandas reads it as NaN and statsmodels then drops
those rows silently, so a reader recomputing the ANOVAs from the sheet fits a three-level factor on 18
samples instead of a four-level factor on 24 and gets different p-values.

The obvious fix - writing "None" - does NOT work: the literal string "None" is in pandas' default
na_values list, so it is read back as NaN and the trap survives. The labels below are not in that list
('NONE' and 'none' are safe; only the exact capitalisation 'None' is not). Verify with:
    python3 -c "from pandas._libs.parsers import STR_NA_VALUES; print('NONE' in STR_NA_VALUES)"

Long-form sheets (antibiotic spelled out) get "No antibiotic"; code sheets (AMX/CIP/ERY) get "NONE",
which is the spelling the already-labelled sheets use. 'UNK' (unknown bin) and 'none (assigned group
not recorded)' mean something else and are left alone. Also refreshes SD15b from
richness_depth_models.py, whose p-values now carry three significant digits.
"""
import shutil
import openpyxl
import pandas as pd

XLSX = "/home/unnot/vault/KJM_info/manuscript/supplementary_tables.xlsx"
TMP = "/tmp/claude-1002/-home-unnot-analyses-KJM/c7eb338e-857b-4770-bd71-9145048c8a2e/scratchpad/st_out.xlsx"
TSV = "/home/unnot/vault/kjm/hiseq/humann/tables/richness_depth_models.tsv"

LONG = {"SD1_samples", "SD2d_pairwise_PERMDISP", "SD2e_alpha_diversity", "SD2f_BC_from_baseline"}
CODE = {"SD12a_butyrate_abundance", "SD14b_total_ARG_by_sample", "SD15_humann_pathways"}
HARMONISE = {"none": "NONE"}          # pandas-safe either way; one spelling across the code sheets
KEEP = {"UNK", "none (assigned group not recorded)"}
UNTREATED = {"CTL", "HFD"}

wb = openpyxl.load_workbook(XLSX)
tf = th = 0
for name in wb.sheetnames:
    ws = wb[name]
    hdr = [c.value for c in ws[1]]
    if "antibiotic" not in hdr:
        continue
    ai = hdr.index("antibiotic") + 1
    gi = hdr.index("group") + 1 if "group" in hdr else None
    label = "No antibiotic" if name in LONG else "NONE"
    filled = harm = 0
    for r in range(2, ws.max_row + 1):
        cell = ws.cell(row=r, column=ai)
        v = cell.value
        if v is None or (isinstance(v, str) and not v.strip()):
            if name not in LONG | CODE:
                continue
            grp = ws.cell(row=r, column=gi).value if gi else None
            if grp is None or str(grp).strip() in UNTREATED:
                cell.value = label
                filled += 1
        elif isinstance(v, str) and v in HARMONISE and v not in KEEP:
            cell.value = HARMONISE[v]
            harm += 1
    if filled or harm:
        print(f"  {name:34s} filled {filled:5d} as {label!r:16s} harmonised {harm:4d}")
    tf += filled
    th += harm
print(f"total: {tf} blanks labelled, {th} relabelled")

b = pd.read_csv(TSV, sep="\t")
ws = wb["SD15b_richness_depth"]
hdr = [c.value for c in ws[1]]
assert hdr == list(b.columns), f"SD15b header changed: {hdr} vs {list(b.columns)}"
old = ws.max_row
for r in range(2, old + 1):
    for c in range(1, len(hdr) + 1):
        ws.cell(row=r, column=c).value = None
for i, row in enumerate(b.itertuples(index=False), start=2):
    for c, v in enumerate(row, start=1):
        ws.cell(row=i, column=c).value = None if pd.isna(v) else v
print(f"SD15b refreshed: {old - 1} -> {len(b)} data rows")

wb.save(TMP)
shutil.copyfile(TMP, XLSX)   # NAS has no atomic replace
print("saved")
