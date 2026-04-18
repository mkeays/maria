#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

def count_unique_mutations( df ) :

    # Select the columns we need and remove duplicate rows caused by alternative transcripts.
    df = df[["icgc_mutation_id", "mutated_from_allele", "mutated_to_allele" ]].drop_duplicates()

    # Get the number of occurences of each combination in a Pandas series object.
    counts = df.groupby( [ "mutated_from_allele", "mutated_to_allele" ] ).size()

    return( counts )

def get_min_max_mutation_sampleids( df ) :

    df = df[[ "icgc_sample_id", "icgc_mutation_id" ]].drop_duplicates()
    
    counts = df.groupby( [ "icgc_sample_id" ] ).size()

    lowest_unique_mutations = counts.idxmin()
    highest_unique_mutations = counts.idxmax()

    return( ( lowest_unique_mutations, highest_unique_mutations ) )

if __name__ == "__main__" :

    filename = sys.argv[1]

    df = pd.read_csv( filename, sep = "\t" )

    mutation_counts = count_unique_mutations( df )
    print( mutation_counts )
    
    min_max_mutation_sampleids = get_min_max_mutation_sampleids( df )
    
    print( f"The sample ID with the lowest unique mutations is {min_max_mutation_sampleids[0]}" )
    print( f"The sample ID with the highest unique mutations is {min_max_mutation_sampleids[1]}" )


