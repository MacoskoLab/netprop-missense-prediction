rule convert_network_to_h5ad:
    input:
        tsv="results/{run}/predicted_perturbed_network.tsv",
    output:
        h5ad="results/{run}/predicted_perturbed_network.h5ad",
    params:
        dataset_id=lambda wildcards: "GSE158067",
        method_id=lambda wildcards: "perturb_algo",
    conda:
        "../envs/convert_adata.yml"
    shell:
        """
        python -c '
import pandas as pd, anndata as ad
df = pd.read_csv("{input.tsv}", sep="\t", names=["regulatoryGene","targetGene","weight"])
adata = ad.AnnData()
adata.uns["dataset_id"] = "{params.dataset_id}"
adata.uns["method_id"]  = "{params.method_id}"
adata.uns["prediction"] = df[["regulatoryGene","targetGene","weight"]].values.tolist()
adata.write_h5ad("{output.h5ad}")
'
        """


rule convert_sc_to_h5ad:
    input:
        expr="resources/GSE158067/GSE158067_gene_exp_mtx_filtered.txt",
        meta="resources/GSE158067/GSE158067_cell_to_perturbation.tsv",
    output:
        h5ad="results/{run}/GSE158067_sc.h5ad",
    params:
        dataset_id=lambda wildcards: "GSE158067",
    conda:
        "../envs/convert_adata.yml"
    shell:
        """
        python -c '
import pandas as pd, anndata as ad
expr = pd.read_csv("{input.expr}", sep="\t", index_col=0).T
meta = pd.read_csv("{input.meta}", sep="\t", index_col=0)
meta = meta.loc[expr.index]
adata = ad.AnnData(X=expr.values, obs=meta.rename(columns={{"Perturbation":"perturbation"}}))
adata.uns["dataset_id"] = "{params.dataset_id}"
adata.layers["X_norm"] = adata.X
adata.write_h5ad("{output.h5ad}")
'
        """


rule evaluate_perturb_with_geneRNIB:
    input:
        pred_h5ad="results/{run}/predicted_perturbed_network.h5ad",
        eval_sc="results/{run}/GSE158067_sc.h5ad",
        tf_all="resources/task_grn_inference/resources_test/grn_benchmark/prior/tf_all.csv",
    output:
        scores="results/{run}/geneRNIB_scores.yaml",
    params:
        temp="results/{run}/tmp",
    shell:
        """
        mkdir -p {params.temp} &&
        resources/task_grn_inference/scripts/single_grn_evaluation.sh \
          {input.pred_h5ad} op \
          --evaluation_data_sc {input.eval_sc} \
          --tf_all {input.tf_all} \
          --temp_dir {params.temp} \
          --score {output.scores}
        """
