#####################################################################
# 
#   Helper functions used by all pipelines to genereate task list
#
#####################################################################

import os
import sys
import traceback
from helpers import runIDHelpers as rid
from helpers import parameters as params
import pipelineHelpers
import bamgineerHelpers
import itertools
import subprocess

chr_list = [12, 13,14]
event_list = ['gain','loss']
haplotype_list=['pat','mat']



##########################################
##########################################
#                                  TASK CREATION                                       #
##########################################
##########################################

#####################################################
#
#  GetProjectNamePath(): Returns <project path>/<name of project>_<runID>
#
#####################################################
def GetProjectNamePathRunID():
   return params.GetProjectPath() +'/'+ params.GetProjectName() +'_'+ rid.GetRunID() + '/'

def GetCancerType():
   return params.GetCancerType() 
   
#####################################################
#
#  ParseInfile(): return a list of strings, one for each line of the infile, with empty lines removed
#
#####################################################
def ParseInfile(infile):
   
    infile_list=[]
    with open(infile) as f:
        for line in f:
            if line.strip():
                line_list = line.rstrip().split()

                if len(line_list) != 2:
                    log = pipelineHelpers.GetLogFile('MAIN')
                    msg = 'Infile does not have 2 columns in line ' + str(line_list)
                    pipelineHelpers.Logging('ERROR', log, msg)
                    raise Exception(msg)
                
                infile_list.append(line_list)
    
    return infile_list
    

#####################################################
#
#  GetSampleIDs(): Returns tumour ids
#
#####################################################
def GetSampleIDs(job_record):
    return job_record[0]

#####################################################
#
#  CheckPath(): Check if path exists and create if it does not
#
#####################################################
def CheckPath(path):
   
   if not os.path.isdir(path):
      os.makedirs(path)
        
   return path
    
#####################################################
#
#  GetProjectPaths(): get all project paths for a given software dir
#
#####################################################
def GetProjectPaths(software):
    project_name_path = GetProjectNamePathRunID()
    #subprocess.call(['chmod', '-R', '+w', project_name_path + software])

    results_path = CheckPath( project_name_path + software + '/')
    intermediate_path = CheckPath( results_path + 'intermediate_' + rid.GetRunID() +'/')
    sentinel_path = CheckPath( project_name_path + 'sentinel_' + rid.GetRunID() +'/')
    
    cancerDir = GetCancerType()
    
    tmpdir = "/".join([intermediate_path,cancerDir, "tmpdir"])
    tmpbams = "/".join([tmpdir, "tmpbams"])
    splittmpbams = "/".join([tmpbams,"splittmpbams"])
    finalbams = "/".join([intermediate_path,cancerDir, "finalbams"])
    patfinalbams = "/".join([finalbams, "pat"])
    matfinalbams = "/".join([finalbams, "mat"])
    
    return (results_path, intermediate_path,sentinel_path,cancerDir,tmpdir,tmpbams,splittmpbams, finalbams,patfinalbams,matfinalbams)
#####################################################
#
#  CheckSentinel(): Checks if sentinel file exists
#
#####################################################
def CheckSentinel(sentinel_list):
    if len(sentinel_list) != 0:
        for sentinel_file in sentinel_list:
            for sent in sentinel_file:
                if sent[0] != '':
                    if os.path.isfile(sent[0]) != True:
                        return False
    return True


#####################################################
#
#  CreateFileList(): Returns filenames and paths in a list
#
#####################################################
def CreateFileList(file_type, num_files, path, flag= None):

    job_list = []
    infile_list = ParseInfile(params.GetInfile())

    for line in infile_list:
      file_list = []
      job_info = line
        
      ## Find the tumour_id for the sample
      sample_id = GetSampleIDs(line)
      
      if (flag is None):     
         if file_type == 'bam':
            file_list.append(line[1])
                
         elif num_files == 1:
               
            file_list.append(path + file_type.format(sample_id))
         elif num_files >= 2: 
                  
               for chr  in chr_list:
                  file_list.append(path + file_type.format('chr', str(chr)))
         else:
            file_list.append('')
         job_list.append(file_list)
         
      elif(flag=="extractROI"):
         pat_gain_event = params.GetPatGainCNV()
         pat_loss_event = params.GetPatLossCNV()
         mat_gain_event = params.GetMatGainCNV()
         mat_loss_event = params.GetMatLossCNV()
         for chr, event, hap   in itertools.product(chr_list, event_list, haplotype_list):
             
            hapev = eval(hap +'_' + event +'_event')
           
            if(not hapev is None):
              
               splittmpbams_hap = "/".join([path, hap])
               file_list.append(splittmpbams_hap + '/'+ 'chr'+file_type.format(chr,event))
               job_list.append(file_list)
               
      elif(flag=="GAINLOSS"):
         pat_gain_event = params.GetPatGainCNV()
         pat_loss_event = params.GetPatLossCNV()
         mat_gain_event = params.GetMatGainCNV()
         mat_loss_event = params.GetMatLossCNV()
         for chr, event, hap   in itertools.product(chr_list, event_list, haplotype_list):
             
            hapev = eval(hap +'_' + event +'_event')
           
            if(not hapev is None):
              
               splittmpbams_hap = "/".join([path, hap])
               file_list.append(splittmpbams_hap + '/'+ 'CHR'+ file_type.format(str(chr).upper(),event.upper()))
               job_list.append(file_list)
      elif(flag=="FINAL"):
         pat_gain_event = params.GetPatGainCNV()
         pat_loss_event = params.GetPatLossCNV()
         mat_gain_event = params.GetMatGainCNV()
         mat_loss_event = params.GetMatLossCNV()
         for chr, event, hap   in itertools.product(chr_list, event_list, haplotype_list):
            
            hapev = eval(hap +'_' + event +'_event')
            
            if(not hapev is None):
               (results_path, intermediate_path,sentinel_path,cancerDir,tmpdir,tmpbams,splittmpbams, finalbams,patfinalbams,matfinalbams) = GetProjectPaths(bamgineerHelpers.name)

               if(event == "gain"):
                  
                  finalbams_hap = "/".join([finalbams, hap])
                  nh= finalbams_hap + '/'+ 'CHR'+ file_type.format(str(chr).upper(),event.upper(),"NH")
                  h= finalbams_hap + '/'+ 'CHR'+ file_type.format(str(chr).upper(),event.upper(),"H")
                  if(os.path.isfile(nh)):
                     file_list.append(nh)
                  if(os.path.isfile(h)):
                     file_list.append(h)
                  job_list.append(file_list)
               elif(event == "loss"):
                  finalbams_hap = "/".join([finalbams, hap])
                  lf = finalbams_hap + '/'+ 'CHR'+ file_type.format(str(chr).upper(),event.upper(),"FINAL")
                  if(os.path.isfile(lf)):
                     file_list.append(lf)
                  job_list.append(file_list)                
    return job_list

########################################################
#
#  CreateTaskList(): Creates a list of tasks by sample

########################################################
def CreateTaskList(inputs, sentinels, outputs, sample_ids, prev_sentinels):
    job_list = []
    count = 0
    
    for sentinel, sample_id in zip(sentinels, sample_ids):
       input_list = []
       output_list = []
       prev_sentinel_list = []
       for inpt in inputs:
           input_list.append(inpt[count])
       for output in outputs:
           output_list.append(output[count])
       for prev_sentinel in prev_sentinels:
           prev_sentinel_list.append(prev_sentinel[count])
      
       mylist = [input_list, sentinel[0], output_list, sample_id[0], prev_sentinel_list]
       job_list.append(mylist)
       #print('appending ' + str(mylist) + '  to job list')
       count = count + 1
    
    return job_list

