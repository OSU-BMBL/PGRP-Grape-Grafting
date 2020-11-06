#!/usr/bin/env python

import json
import multiprocessing
from pathlib import Path
import subprocess
import sys

def main(argv):
    CPUS = str(min([16, multiprocessing.cpu_count() - 1]))
    GENOME_DIR = "GENOME"
    if(len(argv) != 2):
        print("error: 1 parameter expected, " + str(len(argv)-1) +" received")
        return -1
    if(not argv[1].endswith(".json")):
        print("error: expecting file with .json format")
        return -1
    with open(argv[1]) as genome_data:
        data = json.load(genome_data)
        for i in data:
            print("Downloading %s" % (i))
            if(not Path(GENOME_DIR + "/" + Path(i).stem).exists()):
                process = subprocess.Popen("curl " + data[i] + " -O", shell=True, stdout=subprocess.PIPE)
                output, error = process.communicate()
            else:
                print("File already downloaded")
            if(not Path(GENOME_DIR + "/" + Path(i).stem).exists()):
                process = subprocess.Popen("mkdir -p " + GENOME_DIR + " && mv " + i + " " + GENOME_DIR + " && yes n | gunzip " + GENOME_DIR + "/" + i + " && rm -rf " + GENOME_DIR + "/" + i, shell=True, stdout=subprocess.PIPE)
                output, error = process.communicate()
        print("Done")
    return 1

main(sys.argv)