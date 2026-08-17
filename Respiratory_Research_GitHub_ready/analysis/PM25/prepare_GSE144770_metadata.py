#!/usr/bin/env python3
import os
import csv
import re
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(os.environ.get("PROJECT_ROOT", ".")).resolve() / "03_PM25/SECOND_VALIDATION/GSE144770"
ARCHIVE = ROOT / "data/GSE144770_family.xml.tgz"
OUTPUT = ROOT / "data/GSE144770_sample_metadata.csv"
NS = {"m": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}

with tarfile.open(ARCHIVE) as archive:
    xml_file = archive.extractfile("GSE144770_family.xml")
    tree = ET.parse(xml_file)

rows = []
for sample in tree.getroot().findall("m:Sample", NS):
    title = sample.findtext("m:Title", default="", namespaces=NS).strip()
    accession = sample.findtext("m:Accession", default="", namespaces=NS).strip()
    description = sample.findtext("m:Description", default="", namespaces=NS)
    matrix_id = description.strip().splitlines()[-1].strip()
    characteristics = {
        item.attrib.get("tag", ""): (item.text or "").strip()
        for item in sample.findall(".//m:Characteristics", NS)
    }
    donor_match = re.search(r"Donor\s+(\d+)", title)
    donor = int(donor_match.group(1)) if donor_match else None
    treatment = characteristics.get("treatment", "")
    treatment_class = {
        "OECtrl": "OE_control", "0.75OE": "OE_low", "7.5OE": "OE_moderate",
        "75OE": "OE_high", "NIST": "NIST", "WECtrl": "water_control"
    }.get(treatment, treatment)
    dose_step = {"OECtrl": 0, "0.75OE": 1, "7.5OE": 2, "75OE": 3}.get(treatment, "")
    rows.append({
        "matrix_id": matrix_id, "gsm": accession, "title": title, "donor": donor,
        "treatment_code": treatment, "treatment_class": treatment_class,
        "dose_step": dose_step, "asthma_status": characteristics.get("asthma status", ""),
        "ALI_age": characteristics.get("age", ""), "tissue": characteristics.get("tissue", "")
    })

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} samples to {OUTPUT}")
