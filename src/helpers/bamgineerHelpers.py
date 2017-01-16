from helpers import parameters as params
configReader = params.GetConfigReader()

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
    #results_path = configReader.get('RESULTS', 'results_path')
    #return results_path
    return configReader.get('RESULTS', 'results_path')