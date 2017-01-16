"""
inputChecker.py
  check input file for:
        1) tab-delim / right num of columns;
        2) check for sample duplication;
        3) check presence of bam and bai files
  check config file:  files exist for all paths
"""

from ConfigParser import SafeConfigParser
import os


def CheckFileExists(path):
    """returns path of non-existent file or None when file exists"""

    if os.path.exists(path):
        return None
    else:
        return path


def CheckInfile(infile):
    """verify existence of bam and bai files within the input file"""

    result = True
    six_msg = 'PASS'
    ids_msg = 'PASS'
    bam_msg = 'PASS'
    missing_bams = []
    ids = {}

    with open(infile) as f:
        for line in f:
            if line.strip():
                line_list = line.rstrip().split()

                # check for correct number of columns
                if len(line_list) != 2:
                    six_msg = 'FAIL --- Infile does not have 2 columns in line ' + str(line_list)
                    result = False
                    continue

                # check for sample id duplication
                #id_quad = '_'.join([line_list[0], line_list[1], line_list[3], line_list[4]])
                #if id_quad in ids:
                #    ids_msg = 'FAIL --- found duplicate instance of ' + id_quad
                #    result = False
                #else:
                #    ids[id_quad] = True

                # check for existence of bam and bai files
                input_bam = line_list[1]
                #normal_bam = line_list[5]

                missing_bams.append(CheckFileExists(input_bam))
                #missing_bams.append(CheckFileExists(normal_bam))

                tumour_bai = input_bam + '.bai'
                #normal_bai = normal_bam + '.bai'

                missing_bams.append(CheckFileExists(tumour_bai))
                #missing_bams.append(CheckFileExists(normal_bai))


    # filter 'none' instances from the list
    missing_bams = [x for x in missing_bams if x is not None]
    if missing_bams:
        bam_msg = 'FAIL\n\tThe following bam/bai files are missing ' + str(missing_bams)
        result = False

    result_msg = 'Checking input file ---> ' + six_msg
    if six_msg == 'PASS':
        result_msg += '\nChecking for ID uniqueness ---> ' + ids_msg
        result_msg += '\nChecking for existence of bam and bai files ---> ' + bam_msg

    return (result, result_msg)


def CheckConfig(configfile):
    """check existence of all paths in config file"""

    # read in the config file
    config_rdr = SafeConfigParser()
    config_rdr.readfp(open(configfile))

    missing_items = []
    cfg_msg = 'PASS'
    result = True

    # get all params from config file and check their paths
    #  (if they are found in the config_params list
    for s in config_rdr.sections():
        for p in config_rdr.items(s):
            if '_path' in p[0]:
                missing_items.append(CheckFileExists(p[1]))

    missing_items = [x for x in missing_items if x is not None]
    if missing_items:
        cfg_msg = 'FAIL\n\tThe following config file items are missing ' + str(missing_items)
        result = False

    result_msg = 'Checking Config file paths ---> ' + cfg_msg
    return (result, result_msg)


