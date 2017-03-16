
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
        'chr{1}.byname.bam', 3, split_path+"/"))
    outputs.append(taskHelpers.CreateFileList(
            '{0}_{1}_roi.sorted.bam', 12, tmpbams_path,"extractROI")) # max number of outputs chr*events*haplotypes (8 for 2 chromosomes)
    
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
        '{0}{1}_roi.sorted.bam', 12, tmpbams_path,"gain")) 
    outputs.append(taskHelpers.CreateFileList(
        '{0}{1}.re_paired.sorted.bam', 12, tmpbams_path, "gain"))

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
        '{0}_gain.sentinel', 1, sentinel_path))
    sentinels = taskHelpers.CreateFileList(
        '{0}_mutate_gain.sentinel', 1, sentinel_path)
     
    inputs.append(taskHelpers.CreateFileList('{0}{1}_roi.re_paired.sorted.bam', 12, tmpbams_path, "gain"))
    outputs.append(taskHelpers.CreateFileList(
        '{0}{1}.mutated_merged_renamed.sorted.bam', 12, tmpbams_path, "gain"))


    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')
    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                               sample_ids, prev_sentinels)
    for job in job_parameters:
       yield job
 
def subsample_gain_task_list():
    inputs = []
    outputs = []
    prev_sentinels = []
    prev_sentinels.append(taskHelpers.CreateFileList(
        '{0}_mutate_gain.sentinel', 1, sentinel_path))
    sentinels = taskHelpers.CreateFileList(
        '{0}_subsample_gain.sentinel', 1, sentinel_path)
     
    inputs.append(taskHelpers.CreateFileList('{0}{1}_roi.re_paired.mutated_merged_renamed.sorted.bam', 12, tmpbams_path, "gain"))
    outputs.append(taskHelpers.CreateFileList(
        '{0}{1}_GAIN.bam', 12, tmpbams_path, "gain"))

    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')
    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
                                               sample_ids, prev_sentinels)
    for job in job_parameters:
       yield job
  
   
#def mutate_reads_task_list():
#    """populates task inputs and outputs"""
#    inputs = []
#    outputs = []
#    prev_sentinels = []
#    prev_sentinels.append(taskHelpers.CreateFileList(
#        '{0}_findroi.sentinel', 1, sentinel_path))
#    sentinels = taskHelpers.CreateFileList(
#        '{0}_mutatereads.sentinel', 1, sentinel_path) 
#    inputs.append(taskHelpers.CreateFileList(
#        '{0}_{1}_het_roi.sorted.bam', 12, splittmpbams,"extractROI"))
#    inputs.append(taskHelpers.CreateFileList(
#        '{0}_{1}_non_het_roi.sorted.bam', 12, splittmpbams,"extractROI"))
#    outputs.append(taskHelpers.CreateFileList(
#         '{0}_{1}_het_alt_roi.bam', 12, splittmpbams, "extractROI")) # max number of outputs chr*events*haplotypes (8 for 2 chromosomes)   
#    sample_ids = taskHelpers.CreateFileList('{0}', 1, '')
#
#    job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
#                                                sample_ids, prev_sentinels)
#    for job in job_parameters:
#        yield job
      
#def remove_overlap_task_list():
#     """removes overlappig reads"""
#
#     inputs = []
#     outputs = []
#     prev_sentinels = []
#     prev_sentinels.append(taskHelpers.CreateFileList(
#        '{0}_mutatereads.sentinel', 1, sentinel_path))
#     sentinels = taskHelpers.CreateFileList(
#        '{0}_removeoverlap.sentinel', 1, sentinel_path)  
#     inputs.append(taskHelpers.CreateFileList(
#        '{0}_{1}_het_roi.sorted.bam', 12, splittmpbams,"extractROI"))
#     inputs.append(taskHelpers.CreateFileList(
#        '{0}_{1}_non_het_roi.sorted.bam', 12, splittmpbams,"extractROI"))   
#     outputs.append(taskHelpers.CreateFileList(
#         '{0}_{1}_only_het_alt_roi.sorted.bam', 12, splittmpbams, "extractROI")) # max number of outputs chr*events*haplotypes (8 for 2 chromosomes)
#     outputs.append(taskHelpers.CreateFileList(
#         '{0}_{1}_only_non_het_roi.sorted.bam', 12, splittmpbams, "extractROI")) # max number of outputs chr*events*haplotypes (8 for 2 chromosomes)   
#     sample_ids = taskHelpers.CreateFileList('{0}', 1, '')
#
#     job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
#                                                 sample_ids, prev_sentinels)
#     for job in job_parameters:
#        yield job
#    
#
#               
#def sort_merge_task_list(): 
#     inputs = []
#     outputs = []
#     prev_sentinels = []
#     prev_sentinels.append(taskHelpers.CreateFileList(
#        '{0}_gainloss.sentinel', 1, sentinel_path))
#     sentinels = taskHelpers.CreateFileList(
#        '{0}_sortmerge.sentinel', 1, sentinel_path)
#     inputs.append(taskHelpers.CreateFileList(
#         '{0}_{1}_{2}.bam', 12, finalbams, "FINAL"))
#     outputs.append(taskHelpers.CreateFileList(
#         'FINALX.bam', 1, finalbams, "FINAL"))
#     sample_ids = taskHelpers.CreateFileList('{0}', 1, '')
#
#     job_parameters = taskHelpers.CreateTaskList(inputs, sentinels, outputs,
#                                                 sample_ids, prev_sentinels)
#     for job in job_parameters:
#        yield job
# 
#    
