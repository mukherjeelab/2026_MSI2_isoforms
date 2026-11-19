"""Record frozen source URLs, download dates, and SHA256 checksums."""
from pathlib import Path
import csv
import hashlib
import datetime

root = Path(__file__).resolve().parents[2]
sources = {
    'data/FredHutchProteomics.xlsx': 'Local source workbook; upstream filtering metadata not supplied',
    'code/Proetomics.Rmd': 'Local source of exploratory isoform-specific classification rules',
    'data/hpa_25.1/subcellular_location.tsv.zip': 'https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip',
    'data/hpa_25.1/subcellular_locations.tsv.zip': 'https://www.proteinatlas.org/download/tsv/subcellular_locations.tsv.zip',
    'data/hpa_25.1/proteinatlas.tsv.zip': 'https://www.proteinatlas.org/download/proteinatlas.tsv.zip',
    'data/hpa_25.1/subcellular_data_page.html': 'https://www.proteinatlas.org/humanproteome/subcellular/data',
    'data/hpa_25.1/subcellular_overview.html': 'https://www.proteinatlas.org/humanproteome/subcellular',
}
with (root / 'output/proteomics_hpa/source_manifest.tsv').open('w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['path', 'source', 'hpa_version', 'file_mtime_utc', 'bytes', 'sha256'])
    for name, url in sources.items():
        p = root / name
        writer.writerow([name, url, '25.1' if 'hpa_25.1' in name else 'NA',
                         datetime.datetime.fromtimestamp(p.stat().st_mtime, datetime.timezone.utc).isoformat(),
                         p.stat().st_size, hashlib.sha256(p.read_bytes()).hexdigest()])
