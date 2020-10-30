# PGRP-Grape-Grafting

## Set up your environment
```bash
bash setup.sh
```

## Running (dry-run)
```bash
snakemake -n --configfile config.json
```

## Create Rules' Graph
```bash
snakemake --rulegraph --configfile config.json | dot -Tpng -Gdpi=300 > images/rule-graph.png
```

## Create Directed Acyclic Graph (DAG)
```bash
snakemake --dag --configfile config.json | dot -Tpng -Gdpi=300 > images/dag.png
```
