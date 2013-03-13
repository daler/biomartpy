"""
Simple interface for getting a lookup table of gene ID -> other attributes via
biomaRt
"""
import os
from rpy2 import robjects
from rpy2.robjects import r
import pandas

r.library('biomaRt')


def rpy2_to_pandas(rdf, index_col=None):
    """
    Convert rpy2 dataframe representation to a pandas.DataFrame by dumping to
    an intermediate text file.

    Native support probably not coming anytime soon -- see
    https://github.com/pydata/pandas/issues/1448
    """
    intermediate = 'rdf.txt'
    rdf.to_csvfile(path=intermediate, sep='\t', row_names=False)
    df = pandas.read_table(intermediate, index_col=index_col)
    os.unlink(intermediate)
    return df


def list_marts():
    """
    List available marts

    >>> list_marts().ix[0:5]
                   biomart                            version
    0              ensembl       ENSEMBL GENES 69 (SANGER UK)
    1                  snp   ENSEMBL VARIATION 69 (SANGER UK)
    2  functional_genomics  ENSEMBL REGULATION 69 (SANGER UK)
    3                 vega               VEGA 49  (SANGER UK)
    4     bacteria_mart_16       ENSEMBL BACTERIA 16 (EBI UK)
    5        fungi_mart_16          ENSEMBL FUNGI 16 (EBI UK)

    """

    return rpy2_to_pandas(r.listMarts())


def list_datasets(mart_name):
    """
    Returns a pandas.DataFrame listing datasets in mart name

    >>> list_datasets('ensembl').ix[:3]
                       dataset                              description      version
    0   oanatinus_gene_ensembl   Ornithorhynchus anatinus genes (OANA5)        OANA5
    1    tguttata_gene_ensembl  Taeniopygia guttata genes (taeGut3.2.4)  taeGut3.2.4
    2  cporcellus_gene_ensembl          Cavia porcellus genes (cavPor3)      cavPor3
    3  gaculeatus_gene_ensembl   Gasterosteus aculeatus genes (BROADS1)      BROADS1
    """
    return rpy2_to_pandas(r.listDatasets(mart=r.useMart(mart_name)))


def list_attributes(mart_name, dataset):
    """
    Returns a pandas.DataFrame listing attributes for mart name and dataset

    >>> mart_name = 'ensembl'
    >>> dataset = 'dmelanogaster_gene_ensembl'
    >>> list_attributes(mart_name, dataset)[:3]
                        name            description
    0        ensembl_gene_id        Ensembl Gene ID
    1  ensembl_transcript_id  Ensembl Transcript ID
    2     ensembl_peptide_id     Ensembl Protein ID

    """
    dataset = r.useDataset(dataset, mart=r.useMart(mart_name))
    return rpy2_to_pandas(r.listAttributes(dataset))


def list_filters(mart_name, dataset):
    """
    List filters for mart name and dataset

    >>> mart_name = 'ensembl'
    >>> dataset = 'dmelanogaster_gene_ensembl'
    >>> list_filters(mart_name, dataset)[:3]
                  name      description
    0  chromosome_name  Chromosome name
    1            start  Gene Start (bp)
    2              end    Gene End (bp)

    """
    dataset = r.useDataset(dataset, mart=r.useMart(mart_name))
    return rpy2_to_pandas(r.listFilters(dataset))


def make_lookup(mart_name, dataset, attributes):
    """
    Given a mart name, dataset name, and a list of attributes, return
    a pandas.DataFrame indexed by the first attribute in the list provided.

    >>> attributes = ['flybase_gene_id', 'flybasename_gene', 'description']
    >>> df = make_lookup('ensembl', dataset='dmelanogaster_gene_ensembl',
    ... attributes=attributes)

    >>> df.ix['FBgn0031209']
    flybasename_gene                                                Ir21a
    description         Ionotropic receptor 21a [Source:FlyBase gene n...
    Name: FBgn0031209

    >>> df.ix['FBgn0031209'][attributes[1]]
    'Ir21a'
    """
    mart = r.useDataset(dataset, mart=r.useMart(mart_name))
    attributes = robjects.StrVector(attributes)
    results = r.getBM(attributes=attributes, mart=mart)
    return rpy2_to_pandas(results, index_col=0)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
