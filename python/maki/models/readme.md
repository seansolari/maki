# Database design

A database consists of a manifest of sequence accessions and taxonomy IDs for each accession (TSV file), some
metadata stored in a JSON file, and a cluster archive. These clusters are collections of sequences grouped at
a specific taxonomic rank. Each cluster comprises a collection of files and metadata, e.g. sequence data,
annotations, taxonomic data, binaries etc.

A database can opened in read mode. If one doesn't exist, it is created. In read mode, clusters are kept in
archive format, i.e. stored in a .tar.xz file. Specific clusters can be extracted and used for analysis as
required. This approach means that not all clusters are extracted at once, which is important because there
are a lot of clusters and they can be quite large in size.

The database can be converted to edit mode. This means all clusters are exported from the archive and kept on
disk in raw format.

In edit mode, the database (and its clusters) can be modified:

  1. The updating process is abstracted through a `DatabasePackage` abstract class. This class provides:
    - A list of sequence accessions and their taxonomy ID which are to be inserted into the database.
    - A member function `retrieve_data(accessions) -> SequencePackage` which returns an object encompassing data
      (e.g. sequence data, annotation data).
      - One basic implementation of the `DatabasePackage -> SequencePackage` interface could use a csv to 
        maintain lists of accessions, taxids, and locations on disk of all relevent data.
      - However, this interface allows more complex workflows. Chiefly, the `retrieve_data` method can be called
        to e.g. download specific data that is needed, so that large databases can be constructed piecewise,
        without all data needing to be on disk at once.
  2. The user calls `Database.update` supplying a `DatabasePackage`, the clusters are updated, and the database
    manifest is updated to include these sequences.

The database is then coverted to read mode. This means updated clusters are archived into a `.tar.xz` file. In read
mode, the database should be able to list the clusters is currently has in the archive, as keys for clusters that
can be extracted.



## To do

  - Implement DatabasePackage for ATB
