
process RANDIST {
  tag("randist")
  publishDir "${params.outdir ?: 'results'}/randist", mode: 'copy'
  container params.container ?: null
  conda     params.conda_env ?: null

  input:
    path databaseDir
    path treeFile
    path clusterFile
    path countFile

  output:
    path "feature-table.txt.gz", emit: feature_table

  script:
    def threadsArg = params.threads ? "--threads ${params.threads}" : ""
    def alphaArg = params.alpha ? "--alpha ${params.alpha}" : ""
    """
    ${params.tool ?: 'maki'} randist \
      --db ${databaseDir} \
      --tree ${treeFile} \
      --cluster ${clusterFile} \
      --count ${countFile} \
      ${threadsArg} ${alphaArg} \
      --out feature-table.txt.gz
    """
}
