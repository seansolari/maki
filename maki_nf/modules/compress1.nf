
process COMPRESS1 {
  tag("compress1")
  publishDir "${params.outdir ?: 'results'}/compress1", mode: 'copy'
  container params.container ?: null
  conda     params.conda_env ?: null
  cpus params.threads ?: 1

  input:
    path databaseDir

  output:
    path "cooc.kmers.txt.gz", emit: coocc_kmer_count

  script:
    def threadsArg = params.threads ? "--threads ${task.cpus}" : ""
    """
    ${params.tool ?: 'maki'} compress1 \
      ${threadsArg} \
      --db ${databaseDir} \
      --out cooc.kmers.txt.gz
    """
}
