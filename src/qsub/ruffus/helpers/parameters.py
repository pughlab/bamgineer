
from ConfigParser import SafeConfigParser

### global variables ###
configReader = SafeConfigParser()
project_name = ''
project_path = ''
reference = ''
reference_version = ''
infile = ''
num_samples = 1
platform = ''
program_path = ''
#SOROUSH
exonbed = ''
current_path = '' 
pat_gain_event = ''
pat_loss_event = ''
mat_gain_event = ''
mat_loss_event = ''
cancer_type  = ''
vcf_path  = ''



def InitConfigReader(configFile):
    """init the config file"""
    configReader.readfp(open(configFile))


def GetConfigReader():
    """return the configreader"""
    return configReader


def SetProgramPath(prog_path):
    global program_path
    program_path = prog_path


def GetProgramPath():
    return program_path


def SetPlatform(pltfrm):
    global platform
    platform = pltfrm


def GetPlatform():
    return platform


def GetQsubStatement():
    return configReader.get('CLUSTER', 'qsub_statement')


def SetProjectName(pname):
    global project_name
    project_name = pname


def GetProjectName():
    return project_name


def SetProjectPath(ppath):
    global project_path
    project_path = ppath


def GetProjectPath():
    return project_path


def SetInfile(arg_file):
    global infile
    infile = arg_file


def GetInfile():
    return infile


def SetReference(build):
    global reference
    reference = configReader.get('REFERENCE', 'reference_path')


def GetReference():
    return reference


def GetBedFile():
    return configReader.get('REFERENCE','exonbed_path')
    

def SetReferenceVersion(ref_ver):
    global reference_version
    reference_version = ref_ver


def GetReferenceVersion():
    return reference_version


def GetJava():
    return configReader.get('SOFTWARE', 'java')


def GetGATK():
    return configReader.get('SOFTWARE', 'gatk')


def GetProjectName():
    return project_name

###########################
def SetVCF(vcf):
    global vcf_path
    vcf_path = vcf
def GetVCF():
    return vcf_path

def SetExonPath(exon):
    global exon_path
    exon_path = exon
def GetExonPath():
    return exon_path

def SetCancerType(can_type):
    global cancer_type
    cancer_type = can_type
    
def GetCancerType():
    return cancer_type

def SetPatGainCNV(pat_cnv_gain):
    global pat_gain_event
    pat_gain_event = pat_cnv_gain
    
def GetPatGainCNV():
    return pat_gain_event

def SetPatLossCNV(pat_cnv_loss):
    global pat_loss_event
    pat_loss_event = pat_cnv_loss
    
def GetPatLossCNV():
    return pat_loss_event

def SetMatGainCNV(mat_cnv_gain):
    global mat_gain_event
    mat_gain_event = mat_cnv_gain
    
def GetMatGainCNV():
    return mat_gain_event

def SetMatLossCNV(mat_cnv_loss):
    global mat_loss_event
    mat_loss_event = mat_cnv_loss
    
def GetMatLossCNV():
    return mat_loss_event

def SetOutputFileName(out_bam_file):
    global outbamfn
    outbamfn = out_bam_file
    
def GetOutputFileName():
    return outbamfn