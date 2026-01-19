
process LCA1 {
  tag("lca1")
  publishDir "${params.outdir ?: 'results'}/lca1", mode: 'copy'
  container params.container ?: null
  conda     params.conda_env ?: null

  input:
    path databaseDir
    path treeFile
    path coocFile

  output:
    path "cooc.kmers.lca1.txt.gz", emit: coocc_summary

  script:
    """
    ${params.tool ?: 'maki'} lca1 \
      --db ${databaseDir} \
      --tree ${treeFile} \
      --count ${coocFile} \
      --out cooc.kmers.lca1.txt.gz
    """
}
