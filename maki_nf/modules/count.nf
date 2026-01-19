
process COUNT_INDEX {
  publishDir "${params.outdir ?: 'results'}/count", mode: 'copy'
  input:
    path databaseDir
    val whole
  output:
    path(params.count_out ?: 'unique-kmers.txt.gz'), emit: unique_kmers
  script:
    def threadsArg = params.threads ? "--threads ${params.threads}" : ""
    def logArg     = params.log ? "--log ${params.log}" : ""
    def wholeArg   = whole ? "-W" : ""
    def outName    = params.count_out ?: 'unique-kmers.txt.gz'
    """
    ${params.tool ?: 'maki'} count \
      ${threadsArg} ${logArg} \
      --db ${databaseDir} \
      ${wholeArg} \
      --out ${outName}
    """
}

workflow COUNT {
  take:
    database_ch
    whole_ch
  main:
    COUNT_INDEX(database_ch, whole_ch)
  emit:
    unique_kmers = COUNT_INDEX.out.unique_kmers
}
