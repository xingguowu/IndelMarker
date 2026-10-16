#!/usr/bin/env python3
"""Convert Primer3 Boulder-IO output into the final marker-design workbook."""
import argparse
from pathlib import Path
from html import escape
from zipfile import ZIP_DEFLATED, ZipFile


def fasta(path):
    sequences, name, chunks = {}, None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name:
        sequences[name] = "".join(chunks).upper()
    return sequences


def variants(path, reference, flank):
    records = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, *_ = line.rstrip().split("\t")
            if chrom not in reference or "," in alt:
                continue
            pos = int(pos)
            left = max(0, pos - 1 - flank)
            right = min(len(reference[chrom]), pos - 1 + len(ref) + flank)
            target = pos - 1 - left
            template = reference[chrom][left:right]
            marked = template[:target] + "[" + template[target:target + len(ref)] + "]" + template[target + len(ref):]
            records[f"{chrom}_{pos}"] = [chrom, pos, ref, alt, abs(len(ref) - len(alt)), marked]
    return records


def boulder(path):
    entries, current = {}, {}
    if not Path(path).exists() or Path(path).stat().st_size == 0:
        return entries
    with open(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "=":
                identifier = current.get("SEQUENCE_ID")
                if identifier:
                    entries[identifier] = current
                current = {}
            elif "=" in line:
                key, value = line.split("=", 1)
                current[key] = value
    return entries


def primer_values(entry, index):
    left = entry.get(f"PRIMER_LEFT_{index}_SEQUENCE", "")
    right = entry.get(f"PRIMER_RIGHT_{index}_SEQUENCE", "")
    if not left or not right:
        return [""] * 9
    left_pos, left_len = entry.get(f"PRIMER_LEFT_{index}", ",").split(",")
    right_pos, right_len = entry.get(f"PRIMER_RIGHT_{index}", ",").split(",")
    return [left, left_len, entry.get(f"PRIMER_LEFT_{index}_TM", ""), right, right_len,
            entry.get(f"PRIMER_RIGHT_{index}_TM", ""), left_pos, right_pos,
            entry.get(f"PRIMER_PAIR_{index}_PRODUCT_SIZE", "")]


def write_xlsx(path, rows):
    """Write a compact, standards-compliant one-sheet XLSX with no extra dependency."""
    def column_name(index):
        name = ""
        while index:
            index, remainder = divmod(index - 1, 26)
            name = chr(65 + remainder) + name
        return name
    def cell(col, row, value):
        text = escape(str(value or ""))
        return f'<c r="{col}{row}" t="inlineStr"><is><t>{text}</t></is></c>'
    sheet_rows = []
    for row_number, row in enumerate(rows, 1):
        cells = "".join(cell(column_name(i), row_number, value) for i, value in enumerate(row, 1))
        sheet_rows.append(f'<row r="{row_number}">{cells}</row>')
    sheet = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
             '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
             '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" state="frozen"/></sheetView></sheetViews>'
             f'<sheetData>{"".join(sheet_rows)}</sheetData></worksheet>')
    content_types = '<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>'
    rels = '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>'
    workbook = '<?xml version="1.0" encoding="UTF-8"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Primer markers" sheetId="1" r:id="rId1"/></sheets></workbook>'
    workbook_rels = '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>'
    with ZipFile(path, "w", ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types)
        archive.writestr("_rels/.rels", rels)
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        archive.writestr("xl/worksheets/sheet1.xml", sheet)


parser = argparse.ArgumentParser()
parser.add_argument("--vcf", required=True)
parser.add_argument("--reference", required=True)
parser.add_argument("--primer3-output", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--flank", type=int, default=400)
args = parser.parse_args()

records = variants(args.vcf, fasta(args.reference), args.flank)
results = boulder(args.primer3_output)
headers = ["Chromosome", "Position", "Reference", "Alternative", "InDel length (bp)", "Template sequence (target in brackets)"]
for n in range(1, 4):
    headers += [f"Forward primer {n} (5'-3')", "Length", "Tm (°C)", f"Reverse primer {n} (5'-3')", "Length", "Tm (°C)", "Forward start", "Reverse end", "Product size (bp)"]

rows = [headers]
for identifier, metadata in records.items():
    row = list(metadata)
    entry = results.get(identifier, {})
    for n in range(3):
        row += primer_values(entry, n)
    rows.append(row)
write_xlsx(args.output, rows)
