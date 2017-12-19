import os
import sys
import traceback
from helpers import runIDHelpers as rid
from helpers import parameters as params
import itertools
import subprocess

chr_list = range(1,23)
event_list = ['gain','loss']


def GetProjectNamePathRunID():
   results_path = params.GetConfigReader().get('RESULTS', 'results_path')
   cancer_type = params.GetCancerType() 
   
   cancer_dir_path = "/".join([results_path, cancer_type])
   haplotype_path = "/".join([cancer_dir_path, "haplotypedir"])
   tmpbams_path = "/".join([cancer_dir_path, "tmpbams"])
   finalbams_path = "/".join([cancer_dir_path, "finalbams"])
   sentinel_path = CheckPath( cancer_dir_path+'/'+ params.GetProjectName() + "_" + rid.GetRunID() +'/sentinels/')


   return sentinel_path, results_path,haplotype_path,cancer_dir_path,tmpbams_path,finalbams_path


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

def GetSampleIDs(job_record):
    return job_record[0]

def CheckPath(path):
   if not os.path.isdir(path):
      os.makedirs(path)
   return path


def CheckSentinel(sentinel_list):
    if len(sentinel_list) != 0:
        for sentinel_file in sentinel_list:
            for sent in sentinel_file:
                if sent[0] != '':
                    if os.path.isfile(sent[0]) != True:
                        return False
    return True


def CreateFileList(file_type, num_files, path, flag= None):

    sentinel_path, results_path,haplotype_path,cancer_dir_path,tmpbams_path,finalbams_path = GetProjectNamePathRunID()
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
         for chr, event   in itertools.product(chr_list, event_list):
            exonsinroibed = "/".join([haplotype_path,   event + "_exons_in_roi_"+ 'chr'+str(chr) +'.bed'])
            if(os.path.isfile(exonsinroibed)):
               splittmpbams = "/".join([path])
               file_list.append(splittmpbams + '/'+ 'chr'+file_type.format(chr,event))
               job_list.append(file_list)
               
      elif(flag=="gain"):

         for chr  in  chr_list:
            splittmpbams = "/".join([path])
            if(os.path.isfile(splittmpbams+ '/'+ 'chr'+ str(chr)+ '.gain.roi.sorted.bam')):
               file_list.append(splittmpbams + '/'+ 'chr'+ file_type.format(chr,"gain"))
               job_list.append(file_list)
      
      elif(flag=="loss"):
        
         for chr  in  chr_list:
            splittmpbams = "/".join([path])
            if(os.path.isfile(splittmpbams+ '/'+ 'chr'+ str(chr)+ '.loss.roi.sorted.bam')):
               file_list.append(splittmpbams + '/'+ 'chr'+ file_type.format(chr,"loss"))
               job_list.append(file_list)
               
               
      elif(flag=="FINAL"):

            for chr,event in itertools.product(chr_list, event_list):
                chrbam="/".join([finalbams_path,   'CHR'+ str(chr) + '_' +event.upper()+'.bam'])
                sortbyCoord =  "/".join( [params.GetSplitBamsPath(),'chr'+str(chr)+'.bam' ])
                if(os.path.isfile(chrbam)):
                  file_list.append(chrbam)
                  job_list.append(file_list)
                elif(event == 'loss' and sortbyCoord):
                  
                   os.symlink(sortbyCoord, chrbam)
                  
    return job_list


def CreateTaskList(inputs, sentinels, outputs, sample_ids, prev_sentinels):
    job_list = []
    count = 0
    
    for sentinel, sample_id in zip(sentinels, sample_ids):
       input_list = []
       output_list = []
       prev_sentinel_list = []


       for inpt in inputs:
           print(sentinels+ '  ==== ' + str(inpt) + count)
           input_list.append(inpt[count])
           for output in outputs:
               output_list.append(output[count])
           for prev_sentinel in prev_sentinels:
               prev_sentinel_list.append(prev_sentinel[count])

           mylist = [input_list, sentinel[0], output_list, sample_id[0], prev_sentinel_list]
           job_list.append(mylist)
           count = count + 1
    
    return job_list

