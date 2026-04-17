#!/usr/bin/env python3

import sys
import pandas as pd

filename = sys.argv[1]

df = pd.read_csv( filename, sep = "\t" )

# Select the columns we need and remove duplicate rows caused by alternative transcripts.
df = df[["icgc_mutation_id", "mutated_from_allele", "mutated_to_allele" ]].drop_duplicates()

# Get the number of occurences of each combination in a Pandas series object.
counts = df.groupby( [ "mutated_from_allele", "mutated_to_allele" ] ).size()

print( counts )
