# Cyto Batch Integration


<!--
This file is automatically generated from the tasks's api/*.yaml files.
Do not edit this file directly.
-->

A one sentence summary of purpose and methodology. Used for creating an
overview tables.

Repository:
[openproblems-bio/task_cyto_batch_integration](https://github.com/openproblems-bio/task_cyto_batch_integration)

## Description

Provide a clear and concise description of your task, detailing the
specific problem it aims to solve. Outline the input data types, the
expected output, and any assumptions or constraints. Be sure to explain
any terminology or concepts that are essential for understanding the
task.

Explain the motivation behind your proposed task. Describe the
biological or computational problem you aim to address and why it’s
important. Discuss the current state of research in this area and any
gaps or challenges that your task could help address. This section
should convince readers of the significance and relevance of your task.

## Authors & contributors

| name     | roles              |
|:---------|:-------------------|
| John Doe | author, maintainer |

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
`resources_test/task_cyto_batch_integration/starter_file/common_dataset.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'cell_type', 'batch', 'sample', 'donor', 'group'
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

| Name                             | Type   | Description                      |
|:---------------------------------|:-------|:---------------------------------|
| `--input`                        | `file` | A subset of the common dataset.  |
| `--output_unintegrated_censored` | `file` | (*Output*) Unintegrated dataset. |
| `--output_unintegrated`          | `file` | (*Output*) Unintegrated dataset. |
| `--output_validation`            | `file` | (*Output*) Validation dataset.   |

</div>

## File format: Unintegrated Censored

Unintegrated dataset

Example file:
`resources_test/task_cyto_batch_integration/cxg_mouse_pancreas_atlas/train.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'sample', 'donor'
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
| `obs["donor"]` | `string` | (*Optional*) Donor ID. |
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

Unintegrated dataset

Example file:
`resources_test/task_cyto_batch_integration/cxg_mouse_pancreas_atlas/train.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'cell_type', 'batch', 'sample', 'donor', 'group'
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

Validation dataset

Example file:
`resources_test/task_cyto_batch_integration/cxg_mouse_pancreas_atlas/solution.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'cell_type', 'batch', 'sample', 'donor', 'group'
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

A method.

Arguments:

<div class="small">

| Name       | Type   | Description                    |
|:-----------|:-------|:-------------------------------|
| `--input`  | `file` | Unintegrated dataset.          |
| `--output` | `file` | (*Output*) Integrated dataset. |

</div>

## Component type: Control Method

Quality control methods for verifying the pipeline.

Arguments:

<div class="small">

| Name                   | Type   | Description                    |
|:-----------------------|:-------|:-------------------------------|
| `--input_unintegrated` | `file` | Unintegrated dataset.          |
| `--input_validation`   | `file` | Validation dataset.            |
| `--output`             | `file` | (*Output*) Integrated dataset. |

</div>

## Component type: Metric

A task template metric.

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input_validation` | `file` | Validation dataset. |
| `--input_unintegrated` | `file` | Unintegrated dataset. |
| `--input_integrated` | `file` | Integrated dataset. |
| `--output` | `file` | (*Output*) File indicating the score of a metric. |

</div>

## File format: Integrated

Integrated dataset

Example file:
`resources_test/task_cyto_batch_integration/cxg_mouse_pancreas_atlas/prediction.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'label_pred'
     uns: 'dataset_id', 'normalization_id', 'method_id'

</div>

Data structure:

<div class="small">

| Slot                      | Type     | Description                          |
|:--------------------------|:---------|:-------------------------------------|
| `obs["label_pred"]`       | `string` | Predicted labels for the test cells. |
| `uns["dataset_id"]`       | `string` | A unique identifier for the dataset. |
| `uns["normalization_id"]` | `string` | Which normalization was used.        |
| `uns["method_id"]`        | `string` | A unique identifier for the method.  |

</div>

## File format: Score

File indicating the score of a metric.

Example file:
`resources_test/task_cyto_batch_integration/cxg_mouse_pancreas_atlas/score.h5ad`

Format:

<div class="small">

    AnnData object
     uns: 'dataset_id', 'normalization_id', 'method_id', 'metric_ids', 'metric_values'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["normalization_id"]` | `string` | Which normalization was used. |
| `uns["method_id"]` | `string` | A unique identifier for the method. |
| `uns["metric_ids"]` | `string` | One or more unique metric identifiers. |
| `uns["metric_values"]` | `double` | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

</div>

