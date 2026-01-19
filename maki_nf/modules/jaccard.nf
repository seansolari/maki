
process CALC_JACCARD {
  publishDir "${params.outdir ?: 'results'}/jaccard", mode: 'copy'
  input:
    path databaseDir
    path treeFile
    path sharedFile
    path countFile
  output:
    path(params.jaccard ?: 'jaccard-summary.txt.gz'), emit: jaccard_summary
  script:
    def takeArg  = params.jaccard_sample ? "--take ${params.jaccard_sample}" : ""
    def alphaArg = (params.alpha instanceof Integer || params.alpha instanceof Float) ? "--alpha ${params.alpha}" : ""
    def outName  = params.jaccard ?: 'jaccard-summary.txt.gz'
    """
    ${params.tool ?: 'maki'} jaccard \
      --db ${databaseDir} \
      --tree ${treeFile} \
      --dist ${sharedFile} \
      --count ${countFile} \
      ${takeArg} ${alphaArg} \
      --out ${outName}
    """
}

workflow JACCARD {
  take:
    db_ch
    tree_ch
    shared_ch
    count_ch
  main:
    CALC_JACCARD(db_ch, tree_ch, shared_ch, count_ch)
  emit:
    jaccard_summary = CALC_JACCARD.out.jaccard_summary
}
