# Gene annotation tables

`inst/extdata/annotations/GRCh37.csv` and
`inst/extdata/annotations/GRCh38.csv` were manually curated by Malte Nick using
genomic information collected from multiple [NCBI](https://www.ncbi.nlm.nih.gov/)
webpages. The tables support variant-to-gene mapping for the human GRCh37 and
GRCh38 genome assemblies, respectively.

Each semicolon-separated CSV contains the following columns:

| Column | Description |
| --- | --- |
| `Begin` | Genomic start coordinate of the gene |
| `End` | Genomic end coordinate of the gene |
| `Chromosome` | Chromosome identifier |
| `Symbol` | Gene symbol |

clevRact compares variant positions with the inclusive `Begin`--`End` interval
in the table matching the selected genome assembly.
