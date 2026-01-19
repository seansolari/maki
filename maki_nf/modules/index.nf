
process CONSTRUCT_INDEX {
  publishDir "${params.outdir ?: 'results'}/index", mode: 'copy'
  container params.container ?: null
  conda     params.conda_env ?: null

  input:
    path manifest

  output:
    path(params.index_path ?: 'database'), emit: database

  script:
    def threadsArg = params.threads ? "--threads ${params.threads}" : ""
    def ramArg     = params.ram_limit ? "--ram-gb ${params.ram_limit}" : ""
    def bulkArg    = params.bulk ? "-b" : ""
    def suffArg    = params.suffix_size ? "--suffix-size ${params.suffix_size}" : ""
    def logArg     = params.log ? "--log ${params.log}" : ""
    def forceArg   = params.force ? "-F" : ""
    def outName    = params.index_path ?: 'database'
    def kmer       = params.kmer_size
    if (!kmer) exit 1, 'params.kmer_size is required'
    """
    ${params.tool ?: 'maki'} index \
      ${threadsArg} ${ramArg} ${bulkArg} ${suffArg} ${logArg} ${forceArg} \
      --qf ${manifest} \
      --out ${outName} \
      --kmer ${kmer}
    """
}

workflow INDEX {
  take:
    manifest_ch
  main:
    CONSTRUCT_INDEX(manifest_ch)
  emit:
    database = CONSTRUCT_INDEX.out.database
}
