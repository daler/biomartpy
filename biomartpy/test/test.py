from biomartpy import make_lookup

def test_empty_filters():
    # smoke test for issue #2 -- this previously segfaulted.
    x = make_lookup ('ensembl', 'hsapiens_gene_ensembl', 
                     attributes=['ensembl_gene_id', 'external_gene_id'])
