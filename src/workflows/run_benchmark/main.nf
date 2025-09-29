include { checkItemAllowed } from "${meta.resources_dir}/helper.nf"

workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

// construct list of methods and control methods
methods = [
  // shuffle_integration,
  shuffle_integration_by_batch,
  shuffle_integration_by_cell_type,
  harmonypy,
  limma_remove_batch_effect,
  no_integration,
  perfect_integration,
  combat,
  cycombine_no_controls_to_mid,
  cycombine_no_controls_to_goal,
  cycombine_one_control_to_mid,
  cycombine_one_control_to_goal,
  cycombine_all_controls_to_mid,
  cycombine_all_controls_to_goal,
  gaussnorm,
  batchadjust_one_control,
  batchadjust_all_controls,
  cytonorm_no_controls_to_mid,
  cytonorm_all_controls_to_mid,
  cytonorm_one_control_to_mid,
  cytonorm_no_controls_to_goal,
  cytonorm_all_controls_to_goal,
  cytonorm_one_control_to_goal,
  rpca_to_goal,
  rpca_to_mid,
  cytovi,
  mnnpy
]

// construct list of metrics
metrics = [
  emd,
  n_inconsistent_peaks,
  average_batch_r2,
  flowsom_mapping_similarity,
  lisi,
  bras
]

workflow run_wf {
  take:
  input_ch

  main:

  /****************************
   * EXTRACT DATASET METADATA *
   ****************************/
  dataset_ch = input_ch
    // store join id
    | map{ id, state -> 
      [id, state + ["_meta": [join_id: id]]]
    }

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
  method_outputs_ch = dataset_ch

    // run methods on censored split1
    | runEach(
      components: methods,

      // run only non-control methods & filter by method_ids
      filter: { id, state, comp ->
        def method_check = checkItemAllowed(
          comp.config.name,
          state.methods_include,
          state.methods_exclude,
          "methods_include",
          "methods_exclude"
        )
        def method_filter = comp.config.info.type == "method"
        method_check && method_filter
      },

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, comp ->
        id + "." + comp.config.name
      },

      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [ input: "input_censored_split1" ],

      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          method_id: comp.config.name,
          integrated_split1: output.output
        ]
      }
    )

    // run methods on censored split2
    | runEach(
      components: methods,

      // use the 'filter' argument to only run a method on the normalisation the component is asking for
      filter: { id, state, comp ->
        state.method_id == comp.config.name
      },

      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [ input: "input_censored_split2" ],

      // use 'toState' to publish that component's outputs to the overall state
      toState: [ integrated_split2: "output" ]
    )

  control_method_outputs_ch = dataset_ch

    // run control methods on unintegrated data
    | runEach(
      components: methods,

      // run only control methods & filter by method_ids
      filter: { id, state, comp ->
        def method_check = checkItemAllowed(
          comp.config.name,
          state.methods_include,
          state.methods_exclude,
          "methods_include",
          "methods_exclude"
        )
        def method_filter = comp.config.info.type == "control_method"
        method_check && method_filter
      },

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, comp ->
        id + "." + comp.config.name
      },

      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [ input_unintegrated: "input_unintegrated" ],

      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          method_id: comp.config.name,
          integrated_split1: output.output_integrated_split1,
          integrated_split2: output.output_integrated_split2
        ]
      }
    )


  score_ch = method_outputs_ch
    | mix(control_method_outputs_ch)

    // run all metrics
    | runEach(
      components: metrics,
      id: { id, state, comp ->
        id + "." + comp.config.name
      },
      filter: { id, state, comp ->
        // filter by metric_ids
        def metric_check = checkItemAllowed(
          comp.config.name,
          state.metrics_include,
          state.metrics_exclude,
          "metrics_include",
          "metrics_exclude"
        )
        // filter by method_id
        metric_check
      },
      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
        input_unintegrated: "input_unintegrated",
        input_integrated_split1: "integrated_split1", 
        input_integrated_split2: "integrated_split2"
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
