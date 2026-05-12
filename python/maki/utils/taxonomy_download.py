import urllib.request
import tarfile
from pathlib import Path


def download_taxonomy(source: str, output: str):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)

    if source == "ncbi":
        url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        archive = output / "taxdump.tar.gz"

        urllib.request.urlretrieve(url, archive)

        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(output)

    elif source == "gtdb":
        url = "https://data.gtdb.ecogenomic.org/releases/latest/ar53_taxonomy.tsv"
        out_file = output / "gtdb_taxonomy.tsv"

        urllib.request.urlretrieve(url, out_file)

    else:
        raise ValueError("Supported sources: ncbi, gtdb")
