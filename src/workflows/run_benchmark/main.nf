include { checkItemAllowed } from "${meta.resources_dir}/helper.nf"
include { paramsetsFromVariants; expandParamsets; methodMatchesParamset; checkMethodAllowed } from "${meta.resources_dir}/paramset_helper.nf"

workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

// construct list of methods and control methods
methods = [
  shuffle_integration_globally,
  no_integration,
  perfect_integration,
  combat,
  cycombine,
  gaussnorm,
  batchadjust,
  cytonorm,
  harmonypy,
  limma_remove_batch_effect,
  rpca,
]

// construct list of metrics
metrics = [
  emd,
  ratio_consistent_peaks,
  average_batch_r2,
  flowsom_mapping_similarity,
  lisi,
  functional_marker_preservation
]

// serialise data to a yaml file in the temp dir
def writeYamlFile(data, String filename) {
  def file = tempFile(filename)
  file.write(toYamlBlob(data))
  file
}

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

    // expand the channel so parameterised methods run once per paramset.
    // the paramsets are read from the --paramsets file if provided, and
    // default to the method components' info.variants.
    | flatMap { id, state ->
      def method_paramsets = state.paramsets
        ? readYaml(state.paramsets)
        : paramsetsFromVariants(methods)
      expandParamsets(id, state, method_paramsets)
    }

    // run methods on censored split1
    | runEach(
      components: methods,

      // run only non-control methods, match paramset-tagged states to their
      // method & filter by method_ids
      filter: { id, state, comp ->
        comp.config.info.type == "method" &&
          methodMatchesParamset(state, comp.config.name) &&
          checkMethodAllowed(comp.config.name, state.paramset_name, state.methods_include, state.methods_exclude)
      },

      // define a new 'id' by appending the method name and paramset name to the dataset id
      id: { id, state, comp ->
        id + "." + comp.config.name + (state.paramset_name ? "." + state.paramset_name : "")
      },

      // pass the paramset arguments (if any) along with the input
      fromState: { id, state, comp ->
        [ input: state.input_censored_split1 ] + (state.paramset ?: [:])
      },

      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          method_id: comp.config.name,
          integrated_split1: output.output
        ]
      }
    )

    // run the same method with the same paramset on censored split2
    | runEach(
      components: methods,

      filter: { id, state, comp ->
        state.method_id == comp.config.name
      },

      fromState: { id, state, comp ->
        [ input: state.input_censored_split2 ] + (state.paramset ?: [:])
      },

      toState: [ integrated_split2: "output" ]
    )

  control_method_outputs_ch = dataset_ch

    // run control methods on unintegrated data
    | runEach(
      components: methods,

      // run only control methods & filter by method_ids
      filter: { id, state, comp ->
        comp.config.info.type == "control_method" &&
          checkItemAllowed(comp.config.name, state.methods_include, state.methods_exclude, "methods_include", "methods_exclude")
      },

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, comp ->
        id + "." + comp.config.name
      },

      fromState: [ input_unintegrated: "input_unintegrated" ],

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
      // filter by metric_ids
      filter: { id, state, comp ->
        checkItemAllowed(comp.config.name, state.metrics_include, state.metrics_exclude, "metrics_include", "metrics_exclude")
      },
      fromState: [
        input_unintegrated: "input_unintegrated",
        input_integrated_split1: "integrated_split1",
        input_integrated_split2: "integrated_split2"
      ],
      toState: { id, output, state, comp ->
        state + [
          metric_id: comp.config.name,
          metric_output: output.output
        ]
      }
    )

    // extract the scores, tagged with the paramset used for the method
    // (null for control methods and methods without paramsets)
    | extract_uns_metadata.run(
      key: "extract_scores",
      fromState: [input: "metric_output"],
      toState: { id, output, state ->
        def uns = readYaml(output.output).uns
        uns.paramset_name = state.paramset_name
        uns.paramset = state.paramset
        state + [
          score_uns: uns
        ]
      }
    )

    // store the scores in a file
    | joinStates { ids, states ->
      ["output", [output_scores: writeYamlFile(states.collect{it.score_uns}, "score_uns.yaml")]]
    }

  /******************************
   * GENERATE OUTPUT YAML FILES *
   ******************************/
  meta_ch = dataset_ch
    | joinStates { ids, states ->
      // gather the task info and annotate it with the commit and timestamp
      def task_info = readYaml(meta.resources_dir.resolve("_viash.yaml"))
      if (workflow.commitId) {
        task_info.commit = workflow.commitId
      }
      task_info.timestamp = workflow.start.toInstant().truncatedTo(java.time.temporal.ChronoUnit.SECONDS).toString()

      // store the dataset metadata, component configs and task info in files
      def new_state = [
        output_dataset_info: writeYamlFile(states.collect{it.dataset_uns}, "dataset_uns.yaml"),
        output_method_configs: writeYamlFile(methods.collect{it.config}, "method_configs.yaml"),
        output_metric_configs: writeYamlFile(metrics.collect{it.config}, "metric_configs.yaml"),
        output_task_info: writeYamlFile(task_info, "task_info.yaml"),
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
