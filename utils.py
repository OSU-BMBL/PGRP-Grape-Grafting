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

def which(file):
        for path in os.environ["PATH"].split(os.pathsep):
                if os.path.exists(os.path.join(path, file)):
                        return path
        return None
