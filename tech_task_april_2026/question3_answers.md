### Question 3: Database

## Question 3.1

```
select count(*) from gene where id_biotype = '23';
```
Answer: 174

## Question 3.2

```
select ENSEMBL_GENE_ID from gene where GENE_SYMBOL = 'TTTY2';
```
Answer: ENSG00000212855

## Question 3.3

```
select CHROMOSOME, COUNT(ID_GENE) from gene GROUP BY CHROMOSOME;
```
Answer:

| CHROMOSOME | count(id_gene) |
-------------------------------
| 1	| 51 |
| 2	| 25 |
| 3	| 31 |
| 4	| 20 |
| 5	| 25 |
| 6	| 16 |
| 7	| 19 |
| 8	| 25 |
| 9	| 18 |
| 10 | 21 |
| 11 | 25 |
| 12 | 21 |
| 13 | 16 |
| 14 | 18 |
| 15 | 17 |
| 16 | 28 |
| 17 | 27 |
| 18 | 9 |
| 19 | 27 |
| 20 | 16 |
| 21 | 9 |
| 22 | 16 |
| 23 | 16 |
| 24 | 4 |

## Question 3.4

```
select count(distinct transcript.ID_TRANSCRIPT)
from transcript
inner join gene on gene.ID_GENE = transcript.ID_GENE
where gene.GENE_SYMBOL = 'RAI14';
```
Answer: 29

## Question 3.5

```
select transcript.ACCESSION
from transcript
inner join gene on gene.ID_GENE = transcript.ID_GENE
where gene.ENSEMBL_GENE_ID = 'ENSG00000266960'
and transcript.IS_CANONICAL = 'y';
```
Answer: ENST00000586416

## Question 3.6

```
select distinct transcript.ACCESSION
from transcript
inner join gene on gene.ID_GENE = transcript.ID_GENE
where gene.GENE_SYMBOL = 'AK1'
and gene.ID_BIOTYPE = '23'
and transcript.FLAGS = 'gencode_basic';
```
Answer: ENST00000373156, ENST00000373176, ENST00000223836

## Question 3.7

Outer join.

## Question 3.8

I would try creating indexes for the tables using `CREATE INDEX` statements.

## Question 3.9

Use the `UNIQUE` keyword to ensure no values are duplicated.

## Question 3.10

If the `ID_GENE` column is the primary key in the `gene` table I would make
`ID_GENE` in the `transcript` table a foreign key to the `gene` table.
