workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

// construct list of methods and control methods
methods = [
  shuffle_integration,
  shuffle_integration_by_batch,
  shuffle_integration_by_cell_type,
  harmonypy,
  limma_remove_batch_effect,
  no_integration,
  perfect_integration_horizontal,
  perfect_integration_vertical,
  combat,
  cycombine_nocontrols,
  gaussnorm,
  cytonorm_controls,
  mnn
]

// construct list of metrics
metrics = [
  emd,
  n_inconsistent_peaks,
  average_batch_r2,
  flowsom_mapping_similarity,
  cms
]

workflow run_wf {
  take:
  input_ch

  main:

  // input_ch is a channel containing:
  // [id, state]
  // where id is a unique string
  // and state is a dictionary
  // in this case, state:
  // [
  //   input_unintegrated_censored: ...,
  //   input_unintegrated: ...,
  // ]

  /****************************
   * EXTRACT DATASET METADATA *
   ****************************/
  dataset_ch = input_ch
    // store join id
    | map{ id, state -> 
      [id, state + ["_meta": [join_id: id]]]
    }

    // // for your information:
    // | harmonypy.run(
    //   fromState: [input: "input_unintegrated_censored"],
    //   args: [alpha: 10],
    //   toState: [output_method: "output"]
    // )
    // 
    // | harmonypy.run(
    //   fromState: { id, state ->
    //     [
    //       input: state.input_unintegrated_censored,
    //       alpha: 10
    //     ]
    //   },
    //   toState: { id, output, state ->
    //     state + [
    //       output_method: output.output
    //     ]
    //   }
    // )

    // extract the dataset metadata
    | extract_uns_metadata.run(
      fromState: [input: "input_unintegrated"],
      toState: { id, output, state ->
        state + [
          dataset_uns: readYaml(output.output).uns
        ]
      }
    )

  /***************************
   * RUN METHODS AND METRICS *
   ***************************/
  score_ch = dataset_ch

    // run all methods
    | runEach(
      components: methods,

      // use the 'filter' argument to only run a method on the normalisation the component is asking for
      filter: { id, state, comp ->
        !state.method_ids || state.method_ids.contains(comp.config.name)
      },

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, comp ->
        id + "." + comp.config.name
      },

      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: { id, state, comp ->
        if (comp.config.info.type == "control_method") {
          [
            input_unintegrated: state.input_unintegrated,
            input_validation: state.input_validation
          ]
        } else {
          [
            input: state.input_unintegrated_censored
          ]
        }
      },

      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          method_id: comp.config.name,
          method_output: output.output
        ]
      }
    )

    // run all metrics
    | runEach(
      components: metrics,
      id: { id, state, comp ->
        id + "." + comp.config.name
      },
      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
        input_validation: "input_validation", 
        input_unintegrated: "input_unintegrated",
        input_integrated: "method_output",
      ],
      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          metric_id: comp.config.name,
          metric_output: output.output
        ]
      }
    )

    // extract the scores
    | extract_uns_metadata.run(
      key: "extract_scores",
      fromState: [input: "metric_output"],
      toState: { id, output, state ->
        state + [
          score_uns: readYaml(output.output).uns
        ]
      }
    )

    | joinStates { ids, states ->
      // store the scores in a file
      def score_uns = states.collect{it.score_uns}
      def score_uns_yaml_blob = toYamlBlob(score_uns)
      def score_uns_file = tempFile("score_uns.yaml")
      score_uns_file.write(score_uns_yaml_blob)

      ["output", [output_scores: score_uns_file]]
    }

  /******************************
   * GENERATE OUTPUT YAML FILES *
   ******************************/
  // TODO: can we store everything below in a separate helper function?

  // extract the dataset metadata
  meta_ch = dataset_ch
    | joinStates { ids, states ->
      // store the dataset metadata in a file
      def dataset_uns = states.collect{it.dataset_uns}
      def dataset_uns_yaml_blob = toYamlBlob(dataset_uns)
      def dataset_uns_file = tempFile("dataset_uns.yaml")
      dataset_uns_file.write(dataset_uns_yaml_blob)

      // store the method configs in a file
      def method_configs = methods.collect{it.config}
      def method_configs_yaml_blob = toYamlBlob(method_configs)
      def method_configs_file = tempFile("method_configs.yaml")
      method_configs_file.write(method_configs_yaml_blob)

      // store the metric configs in a file
      def metric_configs = metrics.collect{it.config}
      def metric_configs_yaml_blob = toYamlBlob(metric_configs)
      def metric_configs_file = tempFile("metric_configs.yaml")
      metric_configs_file.write(metric_configs_yaml_blob)

      // store the task info in a file
      def viash_file = meta.resources_dir.resolve("_viash.yaml")

      // create output state
      def new_state = [
        output_dataset_info: dataset_uns_file,
        output_method_configs: method_configs_file,
        output_metric_configs: metric_configs_file,
        output_task_info: viash_file,
        _meta: states[0]._meta
      ]

      ["output", new_state]
    }

  // merge all of the output data
  output_ch = score_ch
    | mix(meta_ch)
    | joinStates{ ids, states ->
      def mergedStates = states.inject([:]) { acc, m -> acc + m }
      [ids[0], mergedStates]
    }

  emit:
  output_ch
}
