
nextflow.enable.dsl=2

// Include modules
include { INDEX    } from './modules/index.nf'
include { COUNT    } from './modules/count.nf'
include { DIST     } from './modules/dist.nf'
include { JACCARD  } from './modules/jaccard.nf'
include { COMPRESS1 } from './modules/compress1.nf'
include { LCA1     } from './modules/lca1.nf'
include { RANDIST  } from './modules/randist.nf'

// --------------------
// Workflows (entry points)
// --------------------

// Pipeline: jaccard (index -> count --whole -> dist -> jaccard)
workflow jaccard {
    manifest_ch = channel.fromPath(params.genomes)
    phylogeny_ch = channel.fromPath(params.phylogeny)
    db_ch = INDEX(manifest_ch).database
    uniq_ch = COUNT(db_ch, true).unique_kmers
    shared_ch = DIST(db_ch).shared_kmers
    JACCARD(db_ch, phylogeny_ch, shared_ch, uniq_ch)
}

// Pipeline: random_gene_distances (index -> count --features -> randist)
workflow random_gene_distances {
    manifest_ch = channel.fromPath(params.genomes)
    phylogeny_ch = channel.fromPath(params.phylogeny)
    cluster_ch = channel.fromPath(params.cluster_file)
    db_ch = INDEX(manifest_ch).database
    gene_kmers_ch = COUNT(db_ch, false).unique_kmers
    RANDIST(db_ch, phylogeny_ch, cluster_ch, gene_kmers_ch)
}

// Pipeline: kmer_prevalence (index -> compress1 -> lca1)
workflow kmer_prevalence {
    manifest_ch = channel.fromPath(params.genomes)
    phylogeny_ch = channel.fromPath(params.phylogeny)
    db_ch = INDEX(manifest_ch).database
    cooc_ch = COMPRESS1(db_ch).coocc_kmer_count
    LCA1(db_ch, phylogeny_ch, cooc_ch)
}

// Default entry
workflow {
  if (!params.entry) {
    log.info "Choose an entry via -entry [jaccard|random_gene_distances|kmer_prevalence]"
  }
}
