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

<img src="images/rule-graph.png" height="500" />

## Create Directed Acyclic Graph (DAG)
```bash
snakemake --dag --configfile config.json | dot -Tpng -Gdpi=300 > images/dag.png
```


## HPC cluster execution
### South Dakota State University: Roaring Thunder
```bash
snakemake -j JOBS  \ # maximum number of simultaneous jobs to spawn
	  --configfile config.json # configuration file
          --latency-wait 1000 \ # files latency in seconds
          --cluster-config cluster.sdsu.json \ # cluster configuration file
          --cluster "sbatch --job-name={cluster.name} 
                            --nodes={cluster.nodes} 
                            --ntasks-per-node={cluster.ntasks} 
                            --output={cluster.log} 
                            --partition={cluster.partition} 
                            --time={cluster.time}"
``` 
