
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
inbamfn = ''
results_path=''
java_path =''
beagle_path=''
samtools_path=''
bedtools_path='' 
vcftools_path=''
phase=True
ctDNA=False

def InitConfigReader(configFile):
    """init the config file"""
    configReader.readfp(open(configFile))

def GetConfigReader():
    """return the configreader"""
    return configReader

def GetctDNA():
    return ctDNA
def SetctDNA(ctdna):
    global ctDNA
    ctDNA= ctdna
def GetPhase():
    return phase
def SetPhase(ph):
    global phase
    phase= ph
    
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

def SetJavaPath(path):
    global java_path
    java_path = path
def GetJavaPath():
    return java_path
def SetBeaglePath(path):
    global java_path
    java_path = path
def GetBeaglePath():
    return java_path

def SetInputBam(in_bam_fn):
    global inbamfn
    inbamfn = in_bam_fn
def GetInputBam():
    return inbamfn



def SetSoftwarePath(j_path, b_path, s_path, bd_path, v_path,sb_path):
    
    configReader.set('SOFTWARE','java_path',str(j_path))
    configReader.set('SOFTWARE','beagle_path',str(b_path))
    configReader.set('SOFTWARE','samtools_path',str(s_path))
    configReader.set('SOFTWARE','bedtools_path',str(bd_path))
    configReader.set('SOFTWARE','vcftools_path',str(v_path))
    configReader.set('SOFTWARE','sambamba_path',str(sb_path))

def GetSoftwarePath():
    java_path = configReader.get('SOFTWARE', 'java_path')
    beagle_path = configReader.get('SOFTWARE', 'beagle_path')
    samtools_path = configReader.get('SOFTWARE', 'samtools_path')
    bedtools_path = configReader.get('SOFTWARE', 'bedtools_path')
    vcftools_path =  configReader.get('SOFTWARE', 'vcftools_path')
    sambamba_path =  configReader.get('SOFTWARE', 'sambamba_path')
    
    return java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path


