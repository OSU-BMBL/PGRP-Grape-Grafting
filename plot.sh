snakemake --dag $@ | dot -Tpng -Gdpi=300 > images/dag.png
snakemake --rulegraph $@ | dot -Tpng -Gdpi=300 > images/rule-graph.png
