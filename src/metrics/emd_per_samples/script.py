import anndata as ad
import cytonormpy as cnp

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
    "input_integrated": "resources_test/task_cyto_batch_integration/starter_file/integrated.h5ad",
    "output": "output.h5ad",
    "samples_to_compare": "Tube1_Batch1_WT,Tube1_Batch2_WT",
    "layer": "integrated",
}
meta = {"name": "emd_per_samples"}
## VIASH END

print("Reading input files", flush=True)

adata = ad.read_h5ad(par["input_integrated"])

samples_to_compare = [x.strip() for x in par["samples_to_compare"].split(",")]

layer = par["layer"]

markers_to_assess = adata.var[adata.var["to_correct"]].index.to_numpy()

print("Compute metrics", flush=True)

# have to change the "sample" column to file_name for emd_comparison_from_anndata to work.
# Otherwise the _calculate_emd_per_frame used in cytonormpy will error because they
# harcoded the column file_name and use it in assert.
# See line 176 of https://github.com/TarikExner/CytoNormPy/blob/main/cytonormpy/_evaluation/_emd_utils.py#L173
adata.obs["file_name"] = adata.obs["sample"]

df = cnp.emd_from_anndata(
    adata=adata,
    file_list=samples_to_compare,
    channels=markers_to_assess,
    layer=layer,
    sample_identifier_column="file_name",
)

uns_metric_ids = [f"EMD_per_samples_{x}" for x in df.columns]
uns_metric_values = df.loc["all_cells"].to_numpy()
uns_method_id = adata.uns["method_id"] if "method_id" in adata.uns else "unintegrated"


print("Write output AnnData to file", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "method_id": uns_method_id,
        "sample_ids": samples_to_compare,
        "metric_ids": uns_metric_ids,
        "metric_values": uns_metric_values,
    }
)
output.write_h5ad(par["output"], compression="gzip")