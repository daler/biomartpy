biomartpy
=========
Simple interface to access BioMart from Python (Python -> rpy2 -> R's biomaRt
-> ``pandas.DataFrame``), originally written to get a lookup table of gene IDs
-> various attributes for downstream work...

Install from PyPI::

    $ pip install biomartpy


Or from github::

    $ git clone git@github.com:daler/biomartpy.git
    $ cd biomartpy
    $ python setup.py develop

Choose a mart (use ``list_marts()`` to decide)::

    >>> mart_name = 'ensembl'

Choose a dataset (use ``list_datasets(mart_name)`` to decide)::

    >>> dataset = 'dmelanogaster_gene_ensembl'

Choose some attributes (use ``list_attributes(mart_name, dataset)`` to decide)::

    >>> attributes = ['flybase_gene_id', 'flybasename_gene', 'description']

Get a ``pandas.DataFrame`` as a lookup table, indexed by the first attribute in
the provided list::

    >>> df = make_lookup(mart_name, dataset, attributes=attributes)

``.ix`` to extract rows::

    >>> df.ix['FBgn0031209']
    flybasename_gene                                                Ir21a
    description         Ionotropic receptor 21a [Source:FlyBase gene n...
    Name: FBgn0031209

    >>> df.ix['FBgn0031209']['flybasename_gene']
    'Ir21a'

When providing filters and values, you can either provide them in the way
R expects (filters is a list, values is a list-of-lists with one list for each
filter) or as a more convenient dictionary (here, only geting these IDs, and
only for chromosome 2L)::

    >>> filters = {
    ... 'flybase_gene_id': ['FBgn0031208', 'FBgn0002121', 'FBgn0031209', 'FBgn0051973'],
    ... 'chromosome_name': ['2L']}

Set up attributes (here, including ``chromosome_name`` to make sure results are
correct, but attributes and filters don't have to necessarily match)::

    >>> attributes = ['flybase_gene_id', 'flybasename_gene', 'chromosome_name']

Get data::

    >>> df = make_lookup(
    ... mart_name=mart_name,
    ... dataset=dataset,
    ... attributes=attributes,
    ... filters=filters)

Check results::

    >>> df
                    flybasename_gene chromosome_name
    flybase_gene_id                                 
    FBgn0002121               l(2)gl              2L
    FBgn0031208              CG11023              2L
    FBgn0031209                Ir21a              2L
    FBgn0051973                 Cda5              2L

