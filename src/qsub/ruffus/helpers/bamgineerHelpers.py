from helpers import parameters as params

configReader = params.GetConfigReader()
name = 'bamgineer'


def GetBamgineerMem(mem_type):
    """returns specified value from the configuration file"""
    mem_string = 'bamgineer_mem_' + mem_type
    return configReader.get('CLUSTER', mem_string)


#def GetNumClusters():
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'num_clusters')
#
#
#def GetPseudoCounts():
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'pseudo_counts')
#
#
#def GetTxnExpLen():
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'txn_exp_len')
#
#
#def GetTxnZStrength():
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'txn_z_strength')
#
#
#def GetTitanPloidy():
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'titan_ploidy')
#
#
#def GetEstimatePloidy():
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'estimate_ploidy')
#
#
#def GetMaxIters():          
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'max_iters')
#
#
#def GetNormParamNZero():
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'normal_params_n0')
#
#
#def GetEstimateMethod():
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'normal_estimate_method')
#
#
#def GetGCandMapFiles():
#    """returns specified value from the configuration file"""
#    mapper = configReader.get('REFERENCE', 'map_path')
#    gc_val = configReader.get('REFERENCE', 'gc_path')
#    return mapper, gc_val
#
#
#def GetAlphaHigh():
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'alpha_high')
#
#
#def GetAlphaK():
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'alpha_k')
#
#
#def GetMaxCN():
#    """returns specified value from the configuration file"""
#    return configReader.get('TITAN', 'maxCN')
#
#
#def GetDbSnpRef():
#    """returns specified value from the configuration file"""
#    return configReader.get('REFERENCE', 'dbsnp_path')
