# Jill E. Moore
# Moore Lab - UMass Chan
# ENCODE4 cCRE Pipeline
# December 2024

import sys
import json
import urllib.request
from pathlib import Path


def load_encode_json(accession, local_base="/data/projects/encode/json/exps"):
    """
    Try local JSON first, then ENCODE API.
    accession can be experiment accession or file accession.
    """
    possible_paths = [
        Path(local_base) / f"{accession}.json",
        Path(local_base) / accession / f"{accession}.json",
    ]

    for path in possible_paths:
        if path.exists():
            with open(path) as f:
                return json.load(f)

    url = f"https://www.encodeproject.org/{accession}/?format=json"
    with urllib.request.urlopen(url) as response:
        return json.loads(response.read())


def extract_qc_metrics(file_accession):
    """
    Extract FRiP and filtered read depth from ENCODE QC metrics.
    """
    frip = "NA"
    read_depth = "NA"

    try:
        file_data = load_encode_json(file_accession)
    except Exception as e:
        print(f"WARNING loading file {file_accession}: {e}", file=sys.stderr)
        return frip, read_depth

    filtered_depth = None
    fallback_depth = None

    for qc in file_data.get("quality_metrics", []):
        if frip == "NA" and "frip" in qc:
            frip = qc["frip"]

        processing_stage = qc.get("processing_stage", "filtered").lower()

        depth = None
        for key in [
            "mapped_reads",
            "usable_fragments",
            "total_usable_reads",
            "uniquely_mapped_reads",
            "total_reads",
        ]:
            if key in qc:
                depth = qc[key]
                break

        if depth is None:
            continue

        if processing_stage != "unfiltered":
            if filtered_depth is None:
                filtered_depth = depth
        elif fallback_depth is None:
            fallback_depth = depth

    if filtered_depth is not None:
        read_depth = filtered_depth
    elif fallback_depth is not None:
        read_depth = fallback_depth

    return frip, read_depth


def rep_score(entry):
    """
    Prefer files made from all replicates: 1,2 > single replicate.
    Works for 'replicate' or 'isogenic_replicate' fields.
    """
    reps = (
        entry.get("biological_replicates")
        or entry.get("technical_replicates")
        or entry.get("isogenic_replicates")
        or entry.get("isogenic_replicate")
        or []
    )

    if isinstance(reps, int):
        reps = [reps]

    if isinstance(reps, str):
        reps = [x.strip() for x in reps.split(",")]

    reps = set(map(str, reps))

    # Best: combined replicate file
    if len(reps) >= 2:
        return 2

    # Single-replicate file
    if len(reps) == 1:
        return 1

    return 0


def get_fold_change_bigwig(data, genome):
    """
    Get released fold-change bigWig accession and URL.
    Prefer combined replicate/isogenic replicate files, e.g. 1,2.
    """
    candidates = []

    for entry in data.get("files", []):
        if entry.get("status") != "released":
            continue

        if entry.get("assembly") != genome:
            continue

        if entry.get("file_format") != "bigWig":
            continue

        output_type = str(entry.get("output_type", "")).lower()

        if (
            "fold change" in output_type
            or "fold-change" in output_type
            or "fold enrichment" in output_type
            or "fold-enrichment" in output_type
        ):
            candidates.append(entry)

    if not candidates:
        return "NA", "NA"

    candidates = sorted(
        candidates,
        key=lambda x: (
            rep_score(x),
            x.get("preferred_default", False),
            x.get("file_size", 0) or 0,
        ),
        reverse=True,
    )

    chosen = candidates[0]

    acc = chosen.get("accession", "NA")
    href = chosen.get("href", "NA")

    url = "NA"
    if href != "NA":
        url = "https://www.encodeproject.org" + href

    return acc, url

def Retrieve_Biosample_Summary(dataset, genome):
    try:
        data = load_encode_json(dataset)
    except Exception as e:
        print(f"ERROR loading experiment {dataset}: {e}", file=sys.stderr)
        return

    name = data.get("biosample_summary", "NA")
    lab = data.get("lab", {}).get("title", "NA")

    tags = data.get("internal_tags", [])
    tag = "ENTEx" if "ENTEx" in tags else "NA"

    donor = "NA"
    bio = "NA"
    typ = data.get("biosample_ontology", {}).get("classification", "NA")

    try:
        replicates = data.get("replicates", [])
        if replicates:
            lib = replicates[0].get("library", {})
            biosample = lib.get("biosample", {})
            donor = biosample.get("donor", {}).get("accession", "NA")
            bio = biosample.get("accession", "NA")
    except Exception:
        pass

    bed = "NA"
    frip = "NA"
    read_depth = "NA"

    fold_change_bw, fold_change_bw_url = get_fold_change_bigwig(data, genome)

    # ------------------------------------------------------------
    # First pass:
    # Preferred released narrowPeak
    # ------------------------------------------------------------
    for entry in data.get("files", []):
        if (
            entry.get("file_type") == "bed narrowPeak"
            and entry.get("status") == "released"
            and entry.get("assembly") == genome
            and entry.get("preferred_default", False) is True
        ):
            bed = entry.get("accession", "NA")

            if bed != "NA":
                bed_frip, bed_depth = extract_qc_metrics(bed)

                if bed_frip != "NA":
                    frip = bed_frip

                if bed_depth != "NA":
                    read_depth = bed_depth

            break

    # ------------------------------------------------------------
    # Fallback:
    # Any released narrowPeak
    # ------------------------------------------------------------
    if bed == "NA":
        for entry in data.get("files", []):
            if (
                entry.get("file_type") == "bed narrowPeak"
                and entry.get("status") == "released"
                and entry.get("assembly") == genome
            ):
                bed = entry.get("accession", "NA")

                if bed != "NA":
                    bed_frip, bed_depth = extract_qc_metrics(bed)

                    if bed_frip != "NA":
                        frip = bed_frip

                    if bed_depth != "NA":
                        read_depth = bed_depth

                break

    # ------------------------------------------------------------
    # Search all released files for best depth
    # ------------------------------------------------------------
    if read_depth == "NA":
        for entry in data.get("files", []):
            if entry.get("status") != "released":
                continue

            file_acc = entry.get("accession")
            if not file_acc:
                continue

            _, file_depth = extract_qc_metrics(file_acc)

            if file_depth != "NA":
                read_depth = file_depth
                break

    # ------------------------------------------------------------
    # Search all released files for FRiP
    # ------------------------------------------------------------
    if frip == "NA":
        for entry in data.get("files", []):
            if entry.get("status") != "released":
                continue

            file_acc = entry.get("accession")
            if not file_acc:
                continue

            file_frip, _ = extract_qc_metrics(file_acc)

            if file_frip != "NA":
                frip = file_frip
                break

    print(
        dataset + "\t" +
        bed + "\t" +
        fold_change_bw + "\t" +
        fold_change_bw_url + "\t" +
        name + "\t" +
        "ATAC-seq" + "\t" +
        str(frip) + "\t" +
        str(read_depth) + "\t" +
        lab + "\t" +
        typ + "\t" +
        tag + "\t" +
        donor + "\t" +
        bio
    )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

genome_input = "hg38"

#if len(sys.argv) > 1:
  #  genome_input = sys.argv[1]

if genome_input == "hg38":
    species = "Homo+sapiens"
    genome = "GRCh38"
elif genome_input == "mm10":
    species = "Mus+musculus"
    genome = "mm10"
else:
    raise ValueError("genome must be hg38 or mm10")


print(
    "experiment_accession\t"
    "peak_file_accession\t"
    "fold_change_bigwig_accession\t"
    "fold_change_bigwig_url\t"
    "biosample_summary\t"
    "assay\t"
    "frip\t"
    "read_depth\t"
    "lab\t"
    "biosample_type\t"
    "tag\t"
    "donor\t"
    "biosample"
)


url = (
    "https://www.encodeproject.org/search/?type=Experiment"
    "&status=released"
    "&perturbed=false"
    "&assay_title=ATAC-seq"
    "&replicates.library.biosample.donor.organism.scientific_name=" + species +
    "&format=json"
    "&limit=all"
)

with urllib.request.urlopen(url) as response:
    search_data = json.loads(response.read())

for entry in search_data.get("@graph", []):
    accession = entry.get("accession")

    if accession:
        Retrieve_Biosample_Summary(accession, genome)