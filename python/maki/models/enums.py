
from enum import Enum


class Ranks(str, Enum):
    Species = "species"
    Genus = "genus"
    Family = "family"
    Order = "order"
    Class = "class"
    Phylum = "phylum"
    SuperKingdom = "superkingdom"
    
    
class UpdateMode(str, Enum):
    fixed = "fixed"
    updateable = "updateable"
    

class TaxonomySource(str, Enum):
    gtdb = "gtdb"
    ncbi = "ncbi"
