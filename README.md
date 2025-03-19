# Cyto Batch Integration


<!--
This file is automatically generated from the tasks's api/*.yaml files.
Do not edit this file directly.
-->

Benchmarking of batch integration algorithms for cytometry data.

Repository:
[openproblems-bio/task_cyto_batch_integration](https://github.com/openproblems-bio/task_cyto_batch_integration)

## Description

Cytometry is a non-sequencing single cell profiling technique commonly
used in clinical studies. It is very sensitive to batch effects, which
can lead to biases in the interpretation of the result. Batch
integration algorithms are often used to mitigate this effect.

In this project, we are building a pipeline for reproducible and
continuous benchmarking of batch integration algorithms for cytometry
data. As input, methods require cleaned and normalised (using arc-sinh
or logicle transformation) data with multiple batches, cell type labels,
and biological subjects, with paired samples from a subject profiled
across multiple batches. The batch integrated output must be an
integrated marker by cell matrix stored in Anndata format. All markers
in the input data must be returned, regardless of whether they were
integrated or not. This output is then evaluated using metrics that
assess how well the batch effects were removed and how much biological
signals were preserved.

## Authors & contributors

| name               | roles              |
|:-------------------|:-------------------|
| Luca Leomazzi      | author, maintainer |
| Givanna Putri      | author, maintainer |
| Robrecht Cannoodt  | author             |
| Katrien Quintelier | contributor        |
| Sofie Van Gassen   | contributor        |

## API

``` mermaid
flowchart TB
  file_common_dataset("<a href='https://github.com/openproblems-bio/task_cyto_batch_integration#file-format-common-dataset'>Common Dataset</a>")
  comp_data_processor[/"<a href='https://github.com/openproblems-bio/task_cyto_batch_integration#component-type-data-processor'>Data processor</a>"/]
  file_unintegrated_censored("<a href='https://github.com/openproblems-bio/task_cyto_batch_integration#file-format-unintegrated-censored'>Unintegrated Censored</a>")
  file_unintegrated("<a href='https://github.com/openproblems-bio/task_cyto_batch_integration#file-format-unintegrated'>Unintegrated</a>")
  file_validation("<a href='https://github.com/openproblems-bio/task_cyto_batch_integration#file-format-validation'>Validation</a>")
  comp_method[/"<a href='https://github.com/openproblems-bio/task_cyto_batch_integration#component-type-method'>Method</a>"/]
  comp_control_method[/"<a href='https://github.com/openproblems-bio/task_cyto_batch_integration#component-type-control-method'>Control Method</a>"/]
  comp_metric[/"<a href='https://github.com/openproblems-bio/task_cyto_batch_integration#component-type-metric'>Metric</a>"/]
  file_integrated("<a href='https://github.com/openproblems-bio/task_cyto_batch_integration#file-format-integrated'>Integrated</a>")
  file_score("<a href='https://github.com/openproblems-bio/task_cyto_batch_integration#file-format-score'>Score</a>")
  file_common_dataset---comp_data_processor
  comp_data_processor-->file_unintegrated_censored
  comp_data_processor-->file_unintegrated
  comp_data_processor-->file_validation
  file_unintegrated_censored---comp_method
  file_unintegrated---comp_control_method
  file_unintegrated---comp_metric
  file_validation---comp_control_method
  file_validation---comp_metric
  comp_method-->file_integrated
  comp_control_method-->file_integrated
  comp_metric-->file_score
  file_integrated---comp_metric
```

## File format: Common Dataset

A subset of the common dataset.

Example file:
`resources_test/task_cyto_batch_integration/cyto_spleen_subset/common_dataset.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'cell_type', 'batch', 'sample', 'donor', 'group', 'is_control', 'is_validation'
     var: 'numeric_id', 'channel', 'marker', 'marker_type', 'to_correct'
     layers: 'preprocessed'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `obs["cell_type"]` | `string` | Cell type information. |
| `obs["batch"]` | `string` | Batch information. |
| `obs["sample"]` | `string` | Sample ID. |
| `obs["donor"]` | `string` | Donor ID. |
| `obs["group"]` | `string` | Biological group of the donor. |
| `obs["is_control"]` | `integer` | Whether the sample the cell came from can be used as a control for batch effect correction. 0: cannot be used as a control. \>= 1: can be used as a control. For cells with \>= 1: cells with the same value come from the same donor. Different values indicate different donors. |
| `obs["is_validation"]` | `boolean` | Whether the cell will be used as validation data or not. If FALSE, then the cell will only be included in unintegrated and unintegrated_censored. If TRUE, then the cell will only be included in validation. |
| `var["numeric_id"]` | `integer` | Numeric ID associated with each marker. |
| `var["channel"]` | `string` | The channel / detector of the instrument. |
| `var["marker"]` | `string` | (*Optional*) The marker name associated with the channel. |
| `var["marker_type"]` | `string` | Whether the marker is a functional or lineage marker. |
| `var["to_correct"]` | `boolean` | Whether the marker will be batch corrected. |
| `layers["preprocessed"]` | `double` | preprocessed data, e.g. already compensated, transformed and debris/doublets removed. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["dataset_name"]` | `string` | Nicely formatted name. |
| `uns["dataset_url"]` | `string` | (*Optional*) Link to the original source of the dataset. |
| `uns["dataset_reference"]` | `string` | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]` | `string` | Short description of the dataset. |
| `uns["dataset_description"]` | `string` | Long description of the dataset. |
| `uns["dataset_organism"]` | `string` | (*Optional*) The organism of the sample in the dataset. |

</div>

## Component type: Data processor

A data processor.

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input` | `file` | A subset of the common dataset. |
| `--output_unintegrated_censored` | `file` | (*Output*) An unintegrated dataset with certain columns (cells metadata), such as the donor information, hidden. These columns are intentionally hidden to prevent bias. The batch correction algorithm should not have to rely on these information to properly integrate different batches. This dataset is used as the input for the batch correction algorithm. The cells therein are identical to those in the unintegrated dataset. |
| `--output_unintegrated` | `file` | (*Output*) The complete unintegrated dataset, including all cells’ metadata (columns) from the unintegrated_censored dataset. The cells in this dataset are the same to those in the unintegrated_censored dataset. |
| `--output_validation` | `file` | (*Output*) Hold-out dataset for validation. |

</div>

## File format: Unintegrated Censored

An unintegrated dataset with certain columns (cells metadata), such as
the donor information, hidden. These columns are intentionally hidden to
prevent bias. The batch correction algorithm should not have to rely on
these information to properly integrate different batches. This dataset
is used as the input for the batch correction algorithm. The cells
therein are identical to those in the unintegrated dataset.

Example file:
`resources_test/task_cyto_batch_integration/cyto_spleen_subset/unintegrated_censored.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'sample', 'is_control'
     var: 'numeric_id', 'channel', 'marker', 'marker_type', 'to_correct'
     layers: 'preprocessed'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `obs["batch"]` | `string` | Batch information. |
| `obs["sample"]` | `string` | Sample ID. |
| `obs["is_control"]` | `integer` | Whether the sample the cell came from can be used as a control for batch effect correction. 0: cannot be used as a control. \>= 1: can be used as a control. For cells with \>= 1: cells with the same value come from the same donor. Different values indicate different donors. |
| `var["numeric_id"]` | `integer` | Numeric ID associated with each marker. |
| `var["channel"]` | `string` | The channel / detector of the instrument. |
| `var["marker"]` | `string` | (*Optional*) The marker name associated with the channel. |
| `var["marker_type"]` | `string` | Whether the marker is a functional or lineage marker. |
| `var["to_correct"]` | `boolean` | Whether the marker will be batch corrected. |
| `layers["preprocessed"]` | `double` | preprocessed data, e.g. already compensated, transformed and debris/doublets removed. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["dataset_name"]` | `string` | Nicely formatted name. |
| `uns["dataset_url"]` | `string` | (*Optional*) Link to the original source of the dataset. |
| `uns["dataset_reference"]` | `string` | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]` | `string` | Short description of the dataset. |
| `uns["dataset_description"]` | `string` | Long description of the dataset. |
| `uns["dataset_organism"]` | `string` | (*Optional*) The organism of the sample in the dataset. |

</div>

## File format: Unintegrated

The complete unintegrated dataset, including all cells’ metadata
(columns) from the unintegrated_censored dataset. The cells in this
dataset are the same to those in the unintegrated_censored dataset.

Example file:
`resources_test/task_cyto_batch_integration/cyto_spleen_subset/unintegrated.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'cell_type', 'batch', 'sample', 'donor', 'group', 'is_control'
     var: 'numeric_id', 'channel', 'marker', 'marker_type', 'to_correct'
     layers: 'preprocessed'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `obs["cell_type"]` | `string` | Cell type information. |
| `obs["batch"]` | `string` | Batch information. |
| `obs["sample"]` | `string` | Sample ID. |
| `obs["donor"]` | `string` | Donor ID. |
| `obs["group"]` | `string` | Biological group of the donor. |
| `obs["is_control"]` | `integer` | Whether the sample the cell came from can be used as a control for batch effect correction. 0: cannot be used as a control. \>= 1: can be used as a control. For cells with \>= 1: cells with the same value come from the same donor. Different values indicate different donors. |
| `var["numeric_id"]` | `integer` | Numeric ID associated with each marker. |
| `var["channel"]` | `string` | The channel / detector of the instrument. |
| `var["marker"]` | `string` | (*Optional*) The marker name associated with the channel. |
| `var["marker_type"]` | `string` | Whether the marker is a functional or lineage marker. |
| `var["to_correct"]` | `boolean` | Whether the marker will be batch corrected. |
| `layers["preprocessed"]` | `double` | preprocessed data, e.g. already compensated, transformed and debris/doublets removed. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["dataset_name"]` | `string` | Nicely formatted name. |
| `uns["dataset_url"]` | `string` | (*Optional*) Link to the original source of the dataset. |
| `uns["dataset_reference"]` | `string` | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]` | `string` | Short description of the dataset. |
| `uns["dataset_description"]` | `string` | Long description of the dataset. |
| `uns["dataset_organism"]` | `string` | (*Optional*) The organism of the sample in the dataset. |

</div>

## File format: Validation

Hold-out dataset for validation.

Example file:
`resources_test/task_cyto_batch_integration/cyto_spleen_subset/validation.h5ad`

Description:

Dataset containing cells from samples that were held out for evaluating
batch integration output. The cells that are in this dataset belong to
samples which are not included in the unintegrated or
unintegrated_censored datasets. For example, if samples from donor A are
present in batch 1 and 2, the sample from batch 1 may be used as input
for the batch correction algorithm (and thus present in unintegrated and
unintegrated_censored datasets). The sample from batch 2, may not be
included as an input for the batch correction algorithm, but is needed
to validate whether whether the algorithm managed to correct the batch
effect in batch 2 towards batch 1. This sample will then be included in
this dataset (but not in unintegrated and unintegrated_censored
datasets).

Format:

<div class="small">

    AnnData object
     obs: 'cell_type', 'batch', 'sample', 'donor', 'group', 'is_control'
     var: 'numeric_id', 'channel', 'marker', 'marker_type', 'to_correct'
     layers: 'preprocessed'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `obs["cell_type"]` | `string` | Cell type information. |
| `obs["batch"]` | `string` | Batch information. |
| `obs["sample"]` | `string` | Sample ID. |
| `obs["donor"]` | `string` | Donor ID. |
| `obs["group"]` | `string` | Biological group of the donor. |
| `obs["is_control"]` | `integer` | Whether the sample the cell came from can be used as a control for batch effect correction. 0: cannot be used as a control. \>= 1: can be used as a control. For cells with \>= 1: cells with the same value come from the same donor. Different values indicate different donors. |
| `var["numeric_id"]` | `integer` | Numeric ID associated with each marker. |
| `var["channel"]` | `string` | The channel / detector of the instrument. |
| `var["marker"]` | `string` | (*Optional*) The marker name associated with the channel. |
| `var["marker_type"]` | `string` | Whether the marker is a functional or lineage marker. |
| `var["to_correct"]` | `boolean` | Whether the marker will be batch corrected. |
| `layers["preprocessed"]` | `double` | preprocessed data, e.g. already compensated, transformed and debris/doublets removed. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["dataset_name"]` | `string` | Nicely formatted name. |
| `uns["dataset_url"]` | `string` | (*Optional*) Link to the original source of the dataset. |
| `uns["dataset_reference"]` | `string` | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]` | `string` | Short description of the dataset. |
| `uns["dataset_description"]` | `string` | Long description of the dataset. |
| `uns["dataset_organism"]` | `string` | (*Optional*) The organism of the sample in the dataset. |

</div>

## Component type: Method

A method for integrating batch effects in cytometry data.

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input` | `file` | An unintegrated dataset with certain columns (cells metadata), such as the donor information, hidden. These columns are intentionally hidden to prevent bias. The batch correction algorithm should not have to rely on these information to properly integrate different batches. This dataset is used as the input for the batch correction algorithm. The cells therein are identical to those in the unintegrated dataset. |
| `--output` | `file` | (*Output*) Integrated dataset which batch effect was corrected by an algorithm. |

</div>

## Component type: Control Method

Quality control methods for verifying the pipeline.

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input_unintegrated` | `file` | The complete unintegrated dataset, including all cells’ metadata (columns) from the unintegrated_censored dataset. The cells in this dataset are the same to those in the unintegrated_censored dataset. |
| `--input_validation` | `file` | Hold-out dataset for validation. |
| `--output` | `file` | (*Output*) Integrated dataset which batch effect was corrected by an algorithm. |

</div>

## Component type: Metric

A task template metric.

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input_validation` | `file` | Hold-out dataset for validation. |
| `--input_unintegrated` | `file` | The complete unintegrated dataset, including all cells’ metadata (columns) from the unintegrated_censored dataset. The cells in this dataset are the same to those in the unintegrated_censored dataset. |
| `--input_integrated` | `file` | Integrated dataset which batch effect was corrected by an algorithm. |
| `--output` | `file` | (*Output*) File indicating the score of a metric. |

</div>

## File format: Integrated

Integrated dataset which batch effect was corrected by an algorithm

Example file:
`resources_test/task_cyto_batch_integration/cyto_spleen_subset/integrated.h5ad`

Format:

<div class="small">

    AnnData object
     layers: 'integrated'
     uns: 'dataset_id', 'method_id', 'parameters'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `layers["integrated"]` | `double` | The integrated data as returned by a batch correction method. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["method_id"]` | `string` | A unique identifier for the method. |
| `uns["parameters"]` | `object` | (*Optional*) The parameters used for the integration. |

</div>

## File format: Score

File indicating the score of a metric.

Example file:
`resources_test/task_cyto_batch_integration/cyto_spleen_subset/score.h5ad`

Format:

<div class="small">

    AnnData object
     uns: 'dataset_id', 'method_id', 'metric_ids', 'metric_values'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["method_id"]` | `string` | A unique identifier for the batch correction method. |
| `uns["metric_ids"]` | `string` | One or more unique metric identifiers. |
| `uns["metric_values"]` | `double` | The metric values obtained. Must be of same length as ‘metric_ids’. |

</div>

