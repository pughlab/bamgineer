
import lockfile as lf

# global variable
run_id = 0


def GetRunID():
    """Returns the RunID"""
    return "R%04d" % (run_id)


def IncRunID(project_name, db_dir):
    """Increment the RunID and append new value with project name to the file"""
    database_file = db_dir + '/runID_database.txt'

    # lock the file
    lock = lf.FileLock(database_file)
    while not lock.i_am_locking():
        try:
            # wait up to 20 seconds
            lock.acquire(timeout=20)
        except lf.LockTimeout:
            raise Exception(
                'ERROR: Timed out waiting for file lock at ' + lock.path)

    # get the last run_id from the db file
    rundb = open(database_file, 'r')
    for line in rundb:
        (old_id, old_project) = line.split()

    rundb.close()
    global run_id
    run_id = int(old_id) + 1

    # write the incremented run_id with project name to the db file
    with open(database_file, 'a') as rundb:
        rundb.write(str(run_id) + '\t' + project_name + '\n')

    rundb.close()
    lock.release()

    return


def SetRunID(custom_id, db_dir):
    """set a custom run_id specified in the command line"""
    database_file = db_dir + '/runID_database.txt'

    # lookup the custom_id in the run db file
    rundb = open(database_file, 'r')
    found = False
    for line in rundb:
        (old_id, old_project_name) = line.split()
        if (custom_id == int(old_id)):
            found = True
            break

    rundb.close()

    if not (found):
        raise Exception('run ' + str(custom_id) +
                        ' does not exist -- check runID_database.txt')

    global run_id
    run_id = custom_id

    return
