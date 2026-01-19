
process DIST {
  tag("dist")
  publishDir "${params.outdir ?: 'results'}/dist", mode: 'copy'
  container params.container ?: null
  conda     params.conda_env ?: null
  cpus params.threads ?: 1

  input:
    path databaseDir

  output:
    path "shared-kmers.txt.gz", emit: shared_kmers

  script:
    def threadsArg = params.threads ? "--threads ${task.cpus}" : ""
    def tmpArg     = params.temp_dir ? "--temp-dir ${params.temp_dir}" : ""
    """
    ${params.tool ?: 'maki'} dist \
      ${threadsArg} ${tmpArg} \
      --db ${databaseDir} \
      --out shared-kmers.txt.gz
    """
}
