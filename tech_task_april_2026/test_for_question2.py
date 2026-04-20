#!/usr/bin/env python3
# Test a function with Pytest

import pandas as pd

def get_min_max_mutation_sampleids( df ) :

    df = df[[ "icgc_sample_id", "icgc_mutation_id" ]].drop_duplicates()
    
    counts = df.groupby( [ "icgc_sample_id" ] ).size()

    lowest_unique_mutations = counts.idxmin()
    highest_unique_mutations = counts.idxmax()

    return( ( lowest_unique_mutations, highest_unique_mutations ) )


# Testing with Pytest
def test_get_min_max_mutation_sampleids() :
     # Make up a dataframe with two samples, one with one mutation and one with three.
    df = pd.DataFrame({ "icgc_sample_id" : [ "sample1", "sample2", "sample2", "sample2" ],
                        "icgc_mutation_id" : [ "mut1", "mut2", "mut3", "mut4" ]
                     })
    
    # Check we get the expected results from the above function.
    assert get_min_max_mutation_sampleids( df ) == ( "sample1", "sample2" )

