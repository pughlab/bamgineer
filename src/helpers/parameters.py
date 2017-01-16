
from ConfigParser import SafeConfigParser

### global variables ###
configReader = SafeConfigParser()
project_name = ''
infile = ''
gaincnv_path = ''
losscnv_path = ''
cancer_type  = ''
spltbams_path  = ''
het_path = ''
nonhet_path = ''
outbamfn = ''
results_path=''
#exons_path = ''
#current_path = ''
#program_path = ''
#num_samples = 1
#reference_path = ''


def InitConfigReader(configFile):
    """init the config file"""
    configReader.readfp(open(configFile))

def GetConfigReader():
    """return the configreader"""
    return configReader

#def GetVCF():
#    return configReader.get('REFERENCE','vcf_path')
#
#def GetRef():
#    return configReader.get('REFERENCE','reference_path')
#
#def GetExons():
#    return configReader.get('REFERENCE','exons_path')
#
def GetResultsPath():
    return results_path

def SetResultsPath(path):
    global results_path
    results_path= path 

def GetSplitBamsPath():
    return spltbams_path
def SetSplitBamsPath(spltbams):
    global spltbams_path
    spltbams_path= spltbams

def SetCancerType(can_type):
    global cancer_type
    cancer_type = can_type
    
def GetCancerType():
    return cancer_type

def SetGainCNV(cnv_gain):
    global gaincnv_path
    gaincnv_path = cnv_gain
    
def GetGainCNV():
    return gaincnv_path

def SetLossCNV(cnv_loss):
    global losscnv_path
    losscnv_path= cnv_loss
    
def GetLossCNV():
    return losscnv_path

def SetOutputFileName(out_bam_file):
    global outbamfn
    outbamfn = out_bam_file
    
def GetOutputFileName():
    return outbamfn

def SetLogPath(path):
    global log_path
    log_path = path

def GetLogPath():
    return log_path

def SetHetPath(path):
    global het_path
    het_path = path

def GetHetPath():
    return het_path

def SetNonHetPath(path):
    global nonhet_path
    nonhet_path = path
def GetNonHetPath():
    return nonhet_path





