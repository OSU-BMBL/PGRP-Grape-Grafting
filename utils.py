####### Libraries #######
import glob
import os

####### Util functions #######
def extractFilenames(fullnames,suffix):
        names = []
        for file in fullnames:
            names.append(os.path.basename(file).split(suffix)[0])
        return sorted(names)

def findLibraries(path,prefix,suffix):
	filenames_path = glob.glob(os.path.join(path,prefix) + "*" + suffix)
	names = []
	for file in filenames_path:
	    library = os.path.basename(file).split(suffix)[0]
	    if(library not in names):
		    names.append(library)
	return sorted(names)

def loadGenome(ref):
    if(not ref.endswith(".json")):
        raise ValueError("expecting file with .json format for the reference genome")
    FA = None
    GTF = None
    with open(ref) as genome_data:
        data = json.load(genome_data)
        for i in data.keys():
            if(i.endswith(".fa.gz") and (FA is None)):
                FA = Path(i).stem
            elif((i.endswith(".gtf.gz") or i.endswith(".gff3.gz")) and (GTF is None)):
                GTF = Path(i).stem

    if((FA is None) and (GTF is None)):
        raise ValueError("reference genome file wrongly formatted")
    elif(FA is None):
        raise ValueError("reference genome NOT found")
    elif(GTF is None):
        raise ValueError("gene annotation NOT found")
    return FA, GTF

def verifyGenome(ref,fa,gtf):
    existsFA = os.path.exists(fa)
    existsGTF = os.path.exists(gtf)

    if(not existsFA and not existsGTF):
        raise ValueError("the reference genome and gene annotation file don't exist\n"+
        "please run the following command: \n\tpython download.genome.py " + ref)
    elif(not existsFA):
        raise ValueError("the reference genome file doesn't exist\n"+
        "please run the following command: \n\tpython download.genome.py " + ref)
    elif(not existsGTF):
        raise ValueError("the gene annotation file doesn't exist\n"+
        "please run the following command: \n\tpython download.genome.py " + ref)

def which(file):
        for path in os.environ["PATH"].split(os.pathsep):
                if os.path.exists(os.path.join(path, file)):
                        return path
        return None
