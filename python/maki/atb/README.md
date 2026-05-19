# AllTheBacteria Interface

#### Initialise a sql database

This will automatically download an AllTheBacteria SQLite database from OSF. Supply a specific database version with the `--remote-src` argument.

```{console}
maki.atb init atb.sqlite
```

It also automatically creates file lists for each assembly and Bakta annotation to enable subsequent automated download.

***Warning*** - This file can be large (e.g. ~30Gb).

#### Count samples

Count samples (i.e. assemblies, annotated assemblies) within the SQLite database.

```{console}
maki.atb count --taxon "Bacillus subtilis"
maki.atb count --high-quality --max-contamination 0.1
```

#### List samples

```{console}
# Basic
maki.atb head

# High-quality assemblies
maki.atb head --taxon "Campylobacter jejuni" --high-quality --has-annotation

# Limit output
maki.atb head --taxon "Escherichia coli" --limit 100
```

#### Plan and execute downloads

```{console}
maki.atb download-plan --taxon "Escherichia coli"

# Strict filtering
maki.atb download-plan \
    --taxon "Bacillus subtilis" \
    --high-quality \
    --has-assembly \
    --has-annotation \
    --outfile b_subtilis.plan.csv.gz
```

With this manifest, the download can be executed automatically. As we have requested both assemblies and annotations, the output directory `b_subtilis` will contain both subdirectories `b_subtilis/fa` and `b_subtilis/gff` with assemblies and Bakta annotation files inflated in GFF format, respectively. A sample manifest of all downloaded files, and information in the download manifest, will be written to `b_subtilis/manifest.csv.gz`.

```{console}
maki.atb download b_subtilis.plan.csv.gz --output-dir b_subtilis --concurrency 32
```

These files are sufficient to create a `maki` database. You can of course further filter rows from the CSV download manifest before running the `download` command if you want to do more complex filtering than implemented in this CLI.