"""
Simple interface for getting a lookup table of gene ID -> other attributes via
biomaRt
"""
import os
import rpy2
from rpy2 import robjects
from rpy2.robjects import r
import pandas
import numpy as np

r.library('biomaRt')


def rpy2_to_pandas(rdf, index_col=None):
    """
    Convert rpy2 dataframe representation to a pandas.DataFrame by dumping to
    an intermediate text file.

    Native support probably not coming anytime soon -- see
    https://github.com/pydata/pandas/issues/1448
    """
    intermediate = 'rdf.txt'
    rdf = r['as.data.frame'](rdf)
    rdf.to_csvfile(path=intermediate, sep='\t', row_names=False)
    df = pandas.read_table(intermediate, index_col=index_col)
    os.unlink(intermediate)
    return df


def list_marts():
    """
    List available marts

    >>> list_marts().ix[0:5]
                   biomart                              version
    0              ensembl         ENSEMBL GENES 75 (SANGER UK)
    1                  snp     ENSEMBL VARIATION 75 (SANGER UK)
    2  functional_genomics    ENSEMBL REGULATION 75 (SANGER UK)
    3                 vega                 VEGA 53  (SANGER UK)
    4        fungi_mart_21            ENSEMBL FUNGI 21 (EBI UK)
    5  fungi_variations_21  ENSEMBL FUNGI VARIATION 21 (EBI UK)
    <BLANKLINE>
    [6 rows x 2 columns]
    """

    return rpy2_to_pandas(r.listMarts())


def list_datasets(mart_name, verbose=False):
    """
    Returns a pandas.DataFrame listing datasets in mart name

    >>> list_datasets('ensembl').ix[6:7]
                      dataset                          description  \\
    6  csavignyi_gene_ensembl       Ciona savignyi genes (CSAV2.0)   
    7     fcatus_gene_ensembl  Felis catus genes (Felis_catus_6.2)   
    <BLANKLINE>
               version  
    6          CSAV2.0  
    7  Felis_catus_6.2  
    <BLANKLINE>
    [2 rows x 3 columns]
    """
    return rpy2_to_pandas(
        r.listDatasets(mart=r.useMart(mart_name), verbose=verbose))


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
    <BLANKLINE>
    [3 rows x 2 columns]
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
    <BLANKLINE>
    [3 rows x 2 columns]

    """
    dataset = r.useDataset(dataset, mart=r.useMart(mart_name))
    return rpy2_to_pandas(r.listFilters(dataset))


def make_lookup(mart_name, dataset, attributes, filters=None, values=None,
                unique_rows=True):
    """
    Given a mart name, dataset name, and a list of attributes, return
    a pandas.DataFrame indexed by the first attribute in the list provided.

    In R, filters is a character vector, and values is either a single
    character vector (if only one filter provided) or a list of character
    vectors.

    This function allows `filters` to be a dictionary where keys are filters
    and values are...values.


    >>> mart_name = 'ensembl'
    >>> dataset = 'dmelanogaster_gene_ensembl'
    >>> filters = ['flybase_gene_id', 'chromosome_name']
    >>> attributes = ['flybase_gene_id', 'flybasename_gene', 'chromosome_name']
    >>> values = [['FBgn0031208', 'FBgn0002121', 'FBgn0031209', 'FBgn0051973'], ['2L']]
    >>> df = make_lookup(
    ... mart_name=mart_name,
    ... dataset=dataset,
    ... attributes=attributes,
    ... filters=filters,
    ... values=values)

    Alternatively, make a dictionary of filters: values, in which case you
    don't need to provide `values` separately:

    >>> filters = {
    ... 'flybase_gene_id': ['FBgn0031208', 'FBgn0002121', 'FBgn0031209', 'FBgn0051973'],
    ... 'chromosome_name': ['2L']}

    >>> df2 = make_lookup(
    ... mart_name=mart_name,
    ... dataset=dataset,
    ... attributes=attributes,
    ... filters=filters)

    Confirm that both methods yield identical results:

    >>> assert np.all(df.values == df2.values)

    Check results:

    >>> df.head()
                    flybasename_gene chromosome_name
    flybase_gene_id                                 
    FBgn0002121               l(2)gl              2L
    FBgn0031208              CG11023              2L
    FBgn0031209                Ir21a              2L
    FBgn0051973                 Cda5              2L
    <BLANKLINE>
    [4 rows x 2 columns]


    Indexing by gene ID (or whatever was the first attribute provided):

    >>> df.ix['FBgn0031209']
    flybasename_gene    Ir21a
    chromosome_name        2L
    Name: FBgn0031209, dtype: object



    Extracting data:

    >>> df.ix['FBgn0031209']['flybasename_gene']
    'Ir21a'

    Or get all names:

    >>> df['flybasename_gene']
    flybase_gene_id
    FBgn0002121         l(2)gl
    FBgn0031208        CG11023
    FBgn0031209          Ir21a
    FBgn0051973           Cda5
    Name: flybasename_gene, dtype: object

    """
    mart = r.useDataset(dataset, mart=r.useMart(mart_name))
    attributes = robjects.StrVector(attributes)

    kwargs = dict(
        attributes=attributes,
        uniqueRows=unique_rows,
        mart=mart
    )

    def _filter_and_values_to_RList(d):
        """`d` is a dictionary of filters: values.  Returns a StrVector and
        a ListVector of StrVectors"""
        # Could use ListVector directly with the dict, but want to guarantee
        # positional order of filters and values
        f = robjects.StrVector(d.keys())
        v = robjects.ListVector(
            rpy2.rlike.container.TaggedList(
                d.values(),
                tags=d.keys()
            )
        )
        return f, v

    if isinstance(filters, dict):
        if values is not None:
            raise ValueError("`values` are already specified in the "
                             "`filters` dictionary")
        filter_value_dict = filters
        _filters, _values = _filter_and_values_to_RList(filter_value_dict)
        kwargs['filters'] = _filters
        kwargs['values'] = _values

    elif filters is None:
        if values is not None:
            raise ValueError("`filters` must be specified if `values` "
                             "is specified; alternatively use a dictionary "
                             " for `filters`")

    elif filters and values:
        # values needs to be a list of lists; convert it to one if it's not
        # already
        if not isinstance(values[0], (list, tuple)):
            values = [values]

        # If we got here, then assume filters is a list or tuple
        if len(filters) != len(values):
            raise ValueError('Length of `filters` and `values` must match')

        filter_value_dict = dict(zip(filters, values))
        _filters, _values = _filter_and_values_to_RList(filter_value_dict)
        kwargs['filters'] = _filters
        kwargs['values'] = _values

    else:
        raise ValueError('unhandled case')

    results = r.getBM(**kwargs)
    return rpy2_to_pandas(results, index_col=0)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
