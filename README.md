# PGRP-Grape-Grafting

## Running (dry-run)
```bash
snakemake -n --configfile config.json
```

## Create Rules' Graph
```bash
snakemake --rulegraph --configfile config.json | dot -Tpng -Gdpi=300 > images/rule-graph.png
```
