import os
import subprocess

def read_from_subprocess( arglist ):
    ''' Read line by line from subprocess
    '''

    proc = subprocess.Popen( arglist, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    res = []
    while True:
        l = proc.stdout.readline()
        if l !=  '':
            res.append( l.rstrip() )
        else:
            break
    return res


def read_info_from_batchLine( line ):
    entries = line.split()
    return { "jobID":entries[0], "partition":entries[1], "title":entries[2], "user":entries[3], "status":entries[4], "time":entries[5], "nNodes":entries[6], "worker":entries[7] }

def format_batchInfo( batchOutput ):
    out = []
    for line in batchOutput:
        out.append( read_info_from_batchLine(line) )
    return out

def filter_with_wildcards( str, comp ):
    if comp.startswith("*") and comp.count("*") == 1:                          return str.endswith(comp[1:])
    elif comp.endswith("*") and comp.count("*") == 1:                          return str.startswith(comp[:-1])
    elif comp.endswith("*") and comp.startswith("*") and comp.count("*") == 2: return comp[1:-1] in str
    else:                                                                      return str == comp

def get_batchInfo( jobID=None, partition=None, title=None, user=None, status=None ):

    jobs = format_batchInfo( read_from_subprocess( ["squeue", "-u", os.getenv("USER")] )[1:] )

    if jobID and isinstance( jobID, str ):
        jobs = [job for job in jobs if filter_with_wildcards(job["jobID"], jobID)]
    if jobID and isinstance( jobID, list ):
        jobs = [job for job in jobs if any( [filter_with_wildcards(job["jobID"], j) for j in jobID] )]
    if partition and isinstance( partition, str ) and partition in ["c","m","g"]:
        jobs = [job for job in jobs if job["partition"] == partition]
    if partition and isinstance( partition, list ) and all( [p in ["c","m","g"] for p in partition] ):
        jobs = [job for job in jobs if job["partition"] in partition]
    if user and isinstance( user, str ):
        jobs = [job for job in jobs if user.startswith(job["user"])]
    if status and isinstance( status, str ) and status in ["R","PD","CG"]:
        jobs = [job for job in jobs if job["status"] == status]
    if status and isinstance( status, list ) and all( [s in ["R","PD","CG"] for s in status] ):
        jobs = [job for job in jobs if job["status"] in status]
    if title and isinstance( title, str ):
        jobs = [job for job in jobs if filter_with_wildcards(job["title"], title)]
    if title and isinstance( title, list ):
        jobs = [job for job in jobs if any( [filter_with_wildcards(job["title"], j) for j in title] )]
    if title and isinstance( title, list ):
        jobs = [job for job in jobs if any( [filter_with_wildcards(job["title"], j) for j in title] )]

    return jobs




