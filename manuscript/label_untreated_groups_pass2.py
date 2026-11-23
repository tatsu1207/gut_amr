"""Second pass: replace the literal string 'None' anywhere in the workbook's antibiotic columns.

'None' is in pandas' default na_values, so a cell holding that exact text is read back as NaN just as a
blank cell is. SD12a_butyrate_abundance stored it as text rather than leaving the cell empty, which is
why the first pass (which only looked for blanks) skipped it.
"""
import shutil
import openpyxl

XLSX = "/home/unnot/vault/KJM_info/manuscript/supplementary_tables.xlsx"
TMP = "/tmp/claude-1002/-home-unnot-analyses-KJM/c7eb338e-857b-4770-bd71-9145048c8a2e/scratchpad/st_out2.xlsx"
LONG = {"SD1_samples", "SD2d_pairwise_PERMDISP", "SD2e_alpha_diversity", "SD2f_BC_from_baseline"}
NA_TEXT = {"None", "NaN", "NULL", "null", "nan", "NA", "N/A"}

wb = openpyxl.load_workbook(XLSX)
total = 0
for name in wb.sheetnames:
    ws = wb[name]
    hdr = [c.value for c in ws[1]]
    if "antibiotic" not in hdr:
        continue
    ai = hdr.index("antibiotic") + 1
    label = "No antibiotic" if name in LONG else "NONE"
    n = 0
    for r in range(2, ws.max_row + 1):
        cell = ws.cell(row=r, column=ai)
        if isinstance(cell.value, str) and cell.value.strip() in NA_TEXT:
            cell.value = label
            n += 1
    if n:
        print(f"  {name:34s} {n:5d} cells -> {label!r}")
    total += n
print(f"total: {total} NA-coercible strings replaced")
wb.save(TMP)
shutil.copyfile(TMP, XLSX)
print("saved")
