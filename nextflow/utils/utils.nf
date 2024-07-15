#!/usr/bin/env nextflow

def print_params(description, separator="-"*description.length(), indent=4,
                 params=params, sort=false, nullable=true, skip=['help']) {
  /*
    Print script parameters

    @param description Short pipeline description
    @param separator   Separator used between description and params
                       (default: '-' based on description length)
    @param indent      Number of spaces used to indent the message (default: 4)
    @param params      Map of params (default: script params)
    @param sort        Boolean to sort params alphabetically (default: false)
    @param nullable    Boolean to allow null params; if false, raises an error
                       if any param is null (default: true)
    @param skip        List of params not to print; if null, all params are
                       printed (default: ['help'])
  */
  indent = " " * indent
  log.info "\n${indent}${description}\n${indent}${separator}"

  // ignore specific parameters
  p = params.findAll { it.key !in skip }

  // sort parameters (by default, they are ordered as defined in the script)
  if (sort) { p = p.sort() }

  // get max length of keys
  max = p.collect { it.key.length() }.max()
  
  for (i in p) {
    // print parameter
    log.info "${indent}${i.key.padRight(max)} : ${i.value}"

    if (nullable && i.value == null) {
      // raise error if param is null
      exit 1, "ERROR: parameter --${i.key} not defined"
    }
  }
  log.info ""
}

def print_summary() {
  /*
    Print workflow summary when it finishes
  */
  workflow.onComplete {
    println ( workflow.success ? """
    Workflow summary
    ----------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    Exit status : ${workflow.exitStatus}
    """ : """
    Failed      : ${workflow.errorReport}
    Exit status : ${workflow.exitStatus}
    """
    )
  }
}
