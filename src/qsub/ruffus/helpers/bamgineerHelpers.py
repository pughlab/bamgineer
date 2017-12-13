import os
from helpers import parameters as params

configReader = params.GetConfigReader()
name = 'bamgineer'


def GetBamgineerMem(mem_type):
    """returns specified value from the configuration file"""
    mem_string = 'bamgineer_mem_' + mem_type
    return configReader.get('CLUSTER', mem_string)

def GetExons():
    exons_path  = configReader.get('REFERENCE', 'exons_path')
    return exons_path

def GetRef():
    reference_path = configReader.get('REFERENCE', 'reference_path')
    return reference_path

def GetVCF():
    vcf_path = configReader.get('REFERENCE', 'vcf_path')
    return vcf_path
    
def GetResultsPath():
    results_path = configReader.get('RESULTS', 'results_path')
    return results_path


def GetJavaPath():
    #if config file specfies path use it , otherwise search the system for the tool, if still not found
    try:    
        java_path = configReader.get('SOFTWARE', 'java_path')
        return java_path
    except:
        inpath=tool_loaded('java')
        if(inpath):
            print('User has not defined Java but the tool is in system path: '+ inpath)
            java_path = inpath
            return java_path
            pass
        else:
            print("Couldn't find tool in the path. Exiting the program!")
            return
   
def GetBeaglePath():
    beagle_path = configReader.get('SOFTWARE', 'beagle_path')
    return beagle_path
  
def GetSamtoolsPath():
    try:    
        samtools_path = configReader.get('SOFTWARE', 'samtools_path')
        return samtools_path
    except:
        inpath=tool_loaded('samtools')
        if(inpath):
            print('User has not defined Samtools but the tool is in system path: '+ inpath)
            samtools_path = inpath
            return samtools_path
            pass
        else:
            print("Couldn't find Samtools in the path. Exiting the program!")
            return   
   
def GetBedtoolsPath():
    try:    
        bedtools_path = configReader.get('SOFTWARE', 'bedtools_path')
        return bedtools_path
    except:
        inpath=tool_loaded('bedtools')
        if(inpath):
            print('User has not defined Bedtools but the tool is in system path: '+ inpath)
            bedtools_path = inpath
            return bedtools_path
            pass
        else:
            print("Couldn't find Bedtools in the path. Exiting the program!")
            return     
 
def GetVCFtoolsPath():
    try:    
        vcftools_path = configReader.get('SOFTWARE', 'vcftools_path')
        return vcftools_path
    except:
        inpath=tool_loaded('vcftools')
        if(inpath):
            print('User has not defined VCFtools but the tool is in system path: '+ inpath)
            vcftools_path = inpath
            return vcftools_path
            pass
        else: 
            print("Couldn't find VCFtools in the path. Exiting the program!")
            return     
     
def GetSambambaPath():
    try:    
        sambamba_path = configReader.get('SOFTWARE', 'sambamba_path')
        return sambamba_path
    except:
        inpath=tool_loaded('sambamba')
        if(inpath):
            print('User has not defined Sambamba but the tool is in system path: '+ inpath)
            sambamba_path = inpath
            return sambamba_path
            pass
        else:
            print("Couldn't find Sambamba in the path. Exiting the program!")
            return      
   
def tool_loaded(program):
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    def ext_candidates(fpath):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            for candidate in ext_candidates(exe_file):
                if is_exe(candidate):
                    return candidate

    return False