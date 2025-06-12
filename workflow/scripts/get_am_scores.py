from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import requests
import sys
import pandas as pd
from io import StringIO
from tqdm import tqdm


def get_uniprot_id(
    gene,
    organism_id=9606,  # default to human
):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene} AND organism_id:{organism_id} AND reviewed:true",
        "fields": "accession",
        "format": "tsv",
    }
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    lines = resp.text.strip().split("\n")
    return [line for line in lines[1:]][0]


def fetch_am_score(gene: str, variant: str) -> float:
    uniprot_id = get_uniprot_id(gene)
    url = f"https://alphafold.ebi.ac.uk/api/annotations/{uniprot_id}?annotation_type=MUTAGEN&key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94"
    am_resp = requests.get(url)
    am_resp.raise_for_status()

    path_scores_csv_url = am_resp.json()["annotation"][0]["source_url"]
    if not path_scores_csv_url:
        raise ValueError(f"No AM score found for {uniprot_id} with variant {variant}")
    path_scores_resp = requests.get(path_scores_csv_url)
    path_scores_resp.raise_for_status()

    scores_df = pd.read_csv(StringIO(path_scores_resp.text))
    score = scores_df.loc[scores_df["protein_variant"] == variant, "am_pathogenicity"]
    return float(score.iloc[0])


# Use Snakemake parameters directly
perturbations_list = snakemake.input["perturbations"]
out_path = snakemake.output["am_scores"]

# Read the sample_to_variant file with pandas
sample_df = pd.read_csv(perturbations_list, sep="\t", header=None)

results = []
cache = {}
for idx, row in tqdm(
    sample_df.iterrows(), total=sample_df.shape[0], desc="Downloading AM scores"
):
    sample = row[0]
    gene = row[1]
    variant = row[2]
    cache_key = (gene, variant)
    if cache_key in cache:
        score = cache[cache_key]
    else:
        score = fetch_am_score(gene, variant)
        cache[cache_key] = score
    results.append((sample, gene, variant, score))

results_df = pd.DataFrame(results, columns=["sample", "gene", "variant", "score"])
results_df.to_csv(out_path, sep="\t", index=False, header=False)
print(f"Downloaded AM scores for {len(results)} sample-variant pairs.")
