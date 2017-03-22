
import bamgineerHelpers
import taskHelpers
import runIDHelpers as rid
from helpers import parameters as params
import utils 

from helpers import parameters as params

def split_bam_task_list():
    """populates task inputs and outputs"""
    
    (sentinel_path,results_path,haplotype_path,cancer_dir_path,tmpbams_path, finalbams_path) = taskHelpers.GetProjectNamePathRunID()
    inputs = []
    outputs = []
    prev_sentinels = []
    prev_sentinels.append(
        taskHelpers.CreateFileList('None', -1, sentinel_path))
    split_path = "/".join([results_path, "splitbams"])
    
    utils.createDirectory(split_path)
    sentinels = taskHelpers.CreateFileList(
        '{0}_split.sentinel', 1, sentinel_path)
    inputs.append(taskHelpers.CreateFileList(
        'bam', 1, split_path+"/"))
    outputs.append(taskHelpers.CreateFileList(
        'chr{1}.bam', 2, split_path+"/")) 
    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')
    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                                sample_ids, prev_sentinels)
    for job in job_parameters:
        yield job

def sort_by_name_task_list():
    """populates task inputs and outputs"""
    (sentinel_path,results_path,haplotype_path,cancer_dir_path,tmpbams_path, finalbams_path) = taskHelpers.GetProjectNamePathRunID()
    inputs = []
    outputs = []
    prev_sentinels = []
    split_path = "/".join([results_path, "splitbams"])
    
    prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_split.sentinel', 1, sentinel_path)) 
    sentinels = taskHelpers.CreateFileList(
        '{0}_sortn.sentinel', 1, sentinel_path)
    inputs.append(taskHelpers.CreateFileList(
        'chr{1}.bam', 3, split_path+"/"))
    outputs.append(taskHelpers.CreateFileList(
        'chr{1}.byname.bam', 3 ,split_path+"/"))
    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')
    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                                sample_ids, prev_sentinels)
    for job in job_parameters:
        yield job
        
def find_roi_bam_task_list():
    """populates task inputs and outputs"""
    (sentinel_path,results_path,haplotype_path,cancer_dir_path,tmpbams_path, finalbams_path) = taskHelpers.GetProjectNamePathRunID()
    inputs = []
    outputs = []
    prev_sentinels = []
    split_path = "/".join([results_path, "splitbams"])
    
    prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_sortn.sentinel', 1, sentinel_path))
    sentinels = taskHelpers.CreateFileList(
        '{0}_findroi.sentinel', 1, sentinel_path)
    inputs.append(taskHelpers.CreateFileList(
        '{0}.byname.bam', 3, split_path, "extractROI"))
    outputs.append(taskHelpers.CreateFileList(
            '{0}.{1}.roi.bam', 12, tmpbams_path,"extractROI")) # max number of outputs chr*events*haplotypes (8 for 2 chromosomes)
    
    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')
    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                                sample_ids, prev_sentinels)
    for job in job_parameters:
        yield job

def repair_task_list():    
    """main gain and loss algorithm"""
    (sentinel_path,results_path,haplotype_path,cancer_dir_path,tmpbams_path, finalbams_path) = taskHelpers.GetProjectNamePathRunID()
    inputs = []
    outputs = []
    prev_sentinels = []
    
    prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_findroi.sentinel', 1, sentinel_path))
    sentinels = taskHelpers.CreateFileList(
        '{0}_repair.sentinel', 1, sentinel_path)
     
    inputs.append(taskHelpers.CreateFileList(
        '{0}.{1}.roi.sorted.bam', 12, tmpbams_path,"gain")) 
    outputs.append(taskHelpers.CreateFileList(
        '{0}.{1}.re.paired.bam', 12, tmpbams_path, "gain"))

    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')
    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                               sample_ids, prev_sentinels)
    for job in job_parameters:
       yield job

def mutate_gain_task_list():
    
    (sentinel_path,results_path,haplotype_path,cancer_dir_path,tmpbams_path, finalbams_path) = taskHelpers.GetProjectNamePathRunID()
    inputs = []
    outputs = []
    prev_sentinels = []
    prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_repair.sentinel', 1, sentinel_path))
    sentinels = taskHelpers.CreateFileList(
        '{0}_mutate_gain.sentinel', 1, sentinel_path)
     
    inputs.append(taskHelpers.CreateFileList('{0}.{1}.roi.re.paired.sorted.bam', 12, tmpbams_path, "gain"))
    outputs.append(taskHelpers.CreateFileList(
        '{0}.{1}.mutated.merged.renamed.bam', 12, tmpbams_path, "gain"))

    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')
    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                               sample_ids, prev_sentinels)
    for job in job_parameters:
       yield job
 
def subsample_gain_task_list():
    (sentinel_path,results_path,haplotype_path,cancer_dir_path,tmpbams_path, finalbams_path) = taskHelpers.GetProjectNamePathRunID()
    inputs = []
    outputs = []
    prev_sentinels = []
    prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_mutate_gain.sentinel', 1, sentinel_path))
    sentinels = taskHelpers.CreateFileList(
        '{0}_subsample_gain.sentinel', 1, sentinel_path)
     
    inputs.append(taskHelpers.CreateFileList('{0}.{1}.roi.re.paired.mutated.merged.renamed.sorted.bam', 12, tmpbams_path, "gain"))
    outputs.append(taskHelpers.CreateFileList(
        '{0}{1}_GAIN.bam', 12, tmpbams_path, "gain"))

    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')
    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                               sample_ids, prev_sentinels)
    for job in job_parameters:
       yield job
  