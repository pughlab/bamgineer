import os

def create_directory(path):
   if not os.path.isdir(path):
        os.makedirs(path)

def init(cancerType, logger_local, terminating_local, splitbamdir=None):
    global tmpbams, splitbams, finalbams, logdir,logfile,splittmpbams ,RESULTS,cancerDir, haplotypedir,terminating,logger
    RESULTS = "/".join([os.path.abspath(os.path.dirname(__file__)), 'RESULTS'])
    cancerDir = "/".join([RESULTS, cancerType.upper()])
    haplotypedir = "/".join([cancerDir, "haplotypedir"])
    logdir = "/".join([cancerDir, "logs"])
    logfile = "/".join([logdir,  "DEBUG_LOG.log"])
    tmpbams = "/".join([cancerDir, "tmpbams"])
    splittmpbams = "/".join([tmpbams,"splittmpbams"])
    finalbams = "/".join([cancerDir, "finalbams"])
    
    if(splitbamdir == None):
        splitbams = "/".join([os.path.abspath(os.path.dirname(__file__)), "splitbams"])
    else:
        splitbams = splitbamdir
        
    create_directory(cancerDir)
    create_directory(haplotypedir)
    create_directory(logdir)
    create_directory(tmpbams)
    create_directory(splitbams)
    create_directory(splittmpbams)
    create_directory(finalbams)
    create_directory(RESULTS)

    logger = logger_local
    terminating = terminating_local        

    
    