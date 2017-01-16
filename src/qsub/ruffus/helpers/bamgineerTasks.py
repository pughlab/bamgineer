
import bamgineerHelpers
import taskHelpers
import runIDHelpers as rid

(results_path, intermediate_path,
 sentinel_path,cancerDir,tmpdir,tmpbams,splittmpbams, finalbams,patfinalbams,matfinalbams) = taskHelpers.GetProjectPaths(bamgineerHelpers.name)

#SOROUSH
from helpers import parameters as params
#cancer_type = params.GetCancerType()
#cancerDir = "/".join([intermediate_path, cancer_type.upper()])
#tmpdir = "/".join([cancerDir, "tmpdir/"])
#tmpbams = "/".join([tmpdir, "tmpbams/"])
#splittmpbams = "/".join([tmpbams,"splittmpbams/"])
#finalbams = "/".join([cancerDir, "finalbams/"])
#patfinalbams = "/".join([finalbams, "pat/"])
#matfinalbams = "/".join([finalbams, "mat/"])
#tmpdirpat = "/".join([tmpdir, "pat/"])
#tmpdirmat = "/".join([tmpdir, "mat/"])
#haplotype_list=['pat','mat']
#event_list=['gain','loss']

pat_gain_event = params.GetPatGainCNV()
pat_loss_event = params.GetPatLossCNV()
mat_gain_event = params.GetMatGainCNV()
mat_loss_event = params.GetMatLossCNV()

def SplitBamsTaskList():
    """populates task inputs and outputs"""

    inputs = []
    outputs = []
    prev_sentinels = []

    prev_sentinels.append(
        taskHelpers.CreateFileList('None', -1, sentinel_path))
    
    sentinels = taskHelpers.CreateFileList(
        '{0}_split.sentinel', 1, sentinel_path)

    inputs.append(taskHelpers.CreateFileList(
        'bam', 1, intermediate_path))
    outputs.append(taskHelpers.CreateFileList(
        'chr{1}.bam', 2, intermediate_path))
      
    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')

    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                                sample_ids, prev_sentinels)
    for job in job_parameters:
        yield job

def SortByNameTaskList():
    """populates task inputs and outputs"""

    inputs = []
    outputs = []
    prev_sentinels = []
    
    prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_split.sentinel', 1, sentinel_path))
    
    sentinels = taskHelpers.CreateFileList(
        '{0}_sortn.sentinel', 1, sentinel_path)
    
    inputs.append(taskHelpers.CreateFileList(
        'chr{1}.bam', 3, intermediate_path))
    
    outputs.append(taskHelpers.CreateFileList(
        'chr{1}.byname.bam', 3, intermediate_path))

    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')

    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                                sample_ids, prev_sentinels)
    for job in job_parameters:
        yield job
        
def FindRoiBamTaskList():
    """populates task inputs and outputs"""

    inputs = []
    outputs = []
    prev_sentinels = []
    
    prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_sortn.sentinel', 1, sentinel_path))
    
    sentinels = taskHelpers.CreateFileList(
        '{0}_findroi.sentinel', 1, sentinel_path)
    
    inputs.append(taskHelpers.CreateFileList(
        'chr{1}.byname.bam', 3, intermediate_path))
    
    outputs.append(taskHelpers.CreateFileList(
            '{0}_{1}_het_roi.sorted.bam', 12, splittmpbams,"extractROI")) # max number of outputs chr*events*haplotypes (8 for 2 chromosomes)
    outputs.append(taskHelpers.CreateFileList(
        '{0}_{1}_non_het_roi.sorted.bam', 12, splittmpbams,"extractROI"))
    
    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')

  
    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                                sample_ids, prev_sentinels)
    for job in job_parameters:
        yield job
    
    
def MutateReadsTaskList():
    """populates task inputs and outputs"""

    inputs = []
    outputs = []
    prev_sentinels = []
    
    prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_findroi.sentinel', 1, sentinel_path))
    
    sentinels = taskHelpers.CreateFileList(
        '{0}_mutatereads.sentinel', 1, sentinel_path)
    
    inputs.append(taskHelpers.CreateFileList(
        '{0}_{1}_het_roi.sorted.bam', 12, splittmpbams,"extractROI"))
    inputs.append(taskHelpers.CreateFileList(
        '{0}_{1}_non_het_roi.sorted.bam', 12, splittmpbams,"extractROI"))
    
    outputs.append(taskHelpers.CreateFileList(
         '{0}_{1}_het_alt_roi.bam', 12, splittmpbams, "extractROI")) # max number of outputs chr*events*haplotypes (8 for 2 chromosomes)
        
    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')

    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                                sample_ids, prev_sentinels)
    for job in job_parameters:
        yield job
        

def RemoveOverlappingReadsTaskList():
     """removes overlappig reads"""

     inputs = []
     outputs = []
     prev_sentinels = []
     
     prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_mutatereads.sentinel', 1, sentinel_path))
    
     sentinels = taskHelpers.CreateFileList(
        '{0}_removeoverlap.sentinel', 1, sentinel_path)
    
     inputs.append(taskHelpers.CreateFileList(
        '{0}_{1}_het_roi.sorted.bam', 12, splittmpbams,"extractROI"))
     inputs.append(taskHelpers.CreateFileList(
        '{0}_{1}_non_het_roi.sorted.bam', 12, splittmpbams,"extractROI"))
    
     outputs.append(taskHelpers.CreateFileList(
         '{0}_{1}_only_het_alt_roi.sorted.bam', 12, splittmpbams, "extractROI")) # max number of outputs chr*events*haplotypes (8 for 2 chromosomes)
     outputs.append(taskHelpers.CreateFileList(
         '{0}_{1}_only_non_het_roi.sorted.bam', 12, splittmpbams, "extractROI")) # max number of outputs chr*events*haplotypes (8 for 2 chromosomes)
        
     sample_ids = taskHelpers.CreateFileList('{0}', 1, '')

     job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                                 sample_ids, prev_sentinels)
     for job in job_parameters:
        yield job
    
    
    
def ImplementGainLossTaskList():    
     """main gain and loss algorithm"""

     inputs = []
     outputs = []
     prev_sentinels = []
     
     prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_removeoverlap.sentinel', 1, sentinel_path))
    
     sentinels = taskHelpers.CreateFileList(
        '{0}_gainloss.sentinel', 1, sentinel_path)
    
     inputs.append(taskHelpers.CreateFileList(
         '{0}_{1}_HET_ALT.bam', 12, splittmpbams, "GAINLOSS"))
     
     inputs.append(taskHelpers.CreateFileList(
         '{0}_{1}_NONHET.bam', 12, splittmpbams, "GAINLOSS"))
    
     outputs.append(taskHelpers.CreateFileList(
         '{0}_{1}_FINAL.bam', 12, finalbams, "GAINLOSS"))
       
     sample_ids = taskHelpers.CreateFileList('{0}', 1, '')

     job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                                 sample_ids, prev_sentinels)
     for job in job_parameters:
        yield job
        
        
def SortMergeTaskList(): 
     inputs = []
     outputs = []
     prev_sentinels = []
     
     
     prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_gainloss.sentinel', 1, sentinel_path))
     
     sentinels = taskHelpers.CreateFileList(
        '{0}_sortmerge.sentinel', 1, sentinel_path)
    
     
     inputs.append(taskHelpers.CreateFileList(
         '{0}_{1}_{2}.bam', 12, finalbams, "FINAL"))
     #inputs.append(taskHelpers.CreateFileList(
     #    '{0}_{1}_{2.bam', 4, finalbams, "FINAL"))
     #inputs.append(taskHelpers.CreateFileList(
     #    '{0}_{1}_FINAL.bam', 4, finalbams, "FINAL"))
     #
     
     outputs.append(taskHelpers.CreateFileList(
         'FINALX.bam', 1, finalbams, "FINAL"))
     
     sample_ids = taskHelpers.CreateFileList('{0}', 1, '')

     job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                                 sample_ids, prev_sentinels)
     for job in job_parameters:
        yield job
     
     