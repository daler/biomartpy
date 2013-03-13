biomartpy
=========
Simple interface to access BioMart from Python (Python -> rpy2 -> R's biomaRt),
originally written to get a lookup table of gene IDs -> various attributes for
downstream work...

Install::

    $ pip install biomartpy

Choose a mart (use `list_marts()` to decide)::

    >>> mart_name = 'ensembl'

Choose a dataset (use `list_datasets(mart_name)` to decide)::

    >>> dataset = 'dmelanogaster_gene_ensembl'

Choose some attributes (use `list_attributes(mart_name, dataset)` to decide)::

    >>> attributes = ['flybase_gene_id', 'flybasename_gene', 'description']

Get a pandas.DataFrame as a lookup table, indexed by the first attribute in the
provided list::

    >>> df = make_lookup(mart_name, dataset, attributes=attributes)

`.ix` to extract rows::

    >>> df.ix['FBgn0031209']
    flybasename_gene                                                Ir21a
    description         Ionotropic receptor 21a [Source:FlyBase gene n...
    Name: FBgn0031209

    >>> df.ix['FBgn0031209']['flybasename_gene']
    'Ir21a'
