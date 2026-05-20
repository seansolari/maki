# Database implementation plan

A database is opened in read mode. If one doesn't exist, it is created. This means clusters are kept in archive format.

The database is converted to edit mode. This means clusters are exported from the archive.

In edit mode, the database can be modified:

  1. A `Database.plan_update` function accepts a manifest that outlines genome accessions and taxonomy data, and returns
    a `JobSet` object, which has organised genome accessions into taxonomic clusters.
      - Each cluster must either be new, or correpond to a cluster created in `updatable` mode.
  2. To each `Job` within the `JobSet`, the user must provide a method to retrieve genome data for that set of
    accessions.
      - A simple method will be implemented for the CLI, which assumes files exist locally.
      - In the `maki.atb` submodule, this enables data to be downloaded on-demand (per cluster), rather than
        all at once.
  3. The user submits the `JobSet` to `Database.fulfil`, which is free to execute jobs individually in parallel. When `fulfil`
    completes, each cluster will have been updated.
  4. The completed `JobSetResult` is ingested by the main database to update accessions and taxonomic data.

The database is then coverted to read mode. This means clusters are archived into a `.tar.xz` file.

When classifying, each cluster is extracted from the archive and run against all samples. This means that only one
cluster at a time needs to be decompressed, and makes moving/archiving databases easier. The database should be able
to list the clusters is currently has in the archive.
