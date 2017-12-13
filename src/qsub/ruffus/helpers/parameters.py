
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
java_path =''
beagle_path=''
samtools_path=''
bedtools_path='' 
vcftools_path=''
project_path=''
exons_path=''
current_path = ''
sambamba_path = ''
bamutil_path=''

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

def SetReferenceVersion(ref_ver):
    global reference_version
    reference_version = ref_ver

def GetReferenceVersion():
    return reference_version

def GetGATK():
    return configReader.get('SOFTWARE', 'gatk')

def GetProjectName():
    return project_name

def SetVCF(vcf):
    global vcf_path
    vcf_path = vcf
def GetVCF():
    return vcf_path

def SetExonPath(exon):
    global exons_path
    exons_path = exon
def GetExonPath():
    return exons_path

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

def GetSplitBamsPath():
    return spltbams_path
def SetSplitBamsPath(spltbams):
    global spltbams_path
    spltbams_path= spltbams



def GetResultPath():
    return configReader.get('RESULTS', 'results_path')
def SetResultPath(res_path):
    global results_path
    results_path = res_path
    
def GetSoftwarePath():

    java_path = configReader.get('SOFTWARE', 'java_path')
    beagle_path = configReader.get('SOFTWARE', 'beagle_path')
    samtools_path = configReader.get('SOFTWARE', 'samtools_path')
    bedtools_path = configReader.get('SOFTWARE', 'bedtools_path')
    vcftools_path =  configReader.get('SOFTWARE', 'vcftools_path')
    sambamba_path =  configReader.get('SOFTWARE', 'sambamba_path')
    bamutil_path= configReader.get('SOFTWARE', 'bamutil_path')
    return (java_path, beagle_path,samtools_path, bedtools_path, vcftools_path,sambamba_path,bamutil_path)
 
#def GetProjectPaths():
#    results_path=configReader.get('RESULTS', 'results_path')
#    cancer_type = GetCancerType() 
#    cancer_dir_path = "/".join([results_path, cancer_type])
#    haplotype_path = "/".join([cancer_dir_path, "haplotypedir"])
#    tmpbams_path = "/".join([cancer_dir_path, "tmpbams"])
#    finalbams_path = "/".join([cancer_dir_path, "finalbams"])
#     
#    return (results_path,haplotype_path,cancer_dir_path,tmpbams_path,finalbams_path)



