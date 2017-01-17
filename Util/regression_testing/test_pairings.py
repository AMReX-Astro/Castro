from __future__ import print_function

import os


class Run(object):
    def __init__(self, date=None, status=None, hashes=None):
        self.date = date
        self.status = status
        self.hashes = hashes

    def __lt__(self, other):
        return self.date < other.date
    

def get_run_history(top_dir, changelogs=None):
    """return the list of output directories run over the history of the
       suite and a separate list of the tests run (unique names)

    """

    valid_runs = []

    cwd = os.getcwd()
    os.chdir(top_dir)
    
    # start by finding the list of valid test directories
    for f in os.listdir(top_dir):

        # look for a directory of the form 20* (this will work up until 2099
        if f.startswith("20") and os.path.isdir(f):

            passed = False
            hashes = {}
            
            # look for the status file
            status_file = os.path.normpath("{}/{}.status".format(f, f))

            try:
                sf = open(status_file)
            except:
                #print("unable to open status file: {}".format(status_file))
                continue
            else:
                lines = sf.readlines()
                for line in lines:
                    if "PASSED" in line:
                        passed = True
                        break
                sf.close()

            # now get the hashes
            for cl in changelogs:
                _, code = cl.split(".")
                
                cfile = os.path.normpath("{}/{}".format(f, cl))
                try:
                    cf = open(cfile)
                except:
                    #print("unable to open changelog: {}".format(cfile))
                    continue
                else:
                    lines = cf.readlines()
                    for line in lines:
                        if line.startswith("commit"):
                            _, chash = line.split()
                            hashes[code] = chash
                            break
                    cf.close()

            if passed and len(hashes) == len(changelogs):
                valid_runs.append(Run(date=f, status="passed", hashes=hashes))


    valid_runs.sort()
    valid_runs.reverse()

    os.chdir(cwd)
    
    return valid_runs




top_dir = "/home/www/Castro/test-suite/test-suite-gfortran/"

changelogs_needed = ["ChangeLog.Castro", "ChangeLog.BoxLib", "ChangeLog.Microphysics"]

k = [q.split(".")[1] for q in changelogs_needed]

runs = get_run_history(top_dir, changelogs=changelogs_needed)

print("date, status, ", " ".join(k))

for r in runs:
    hash_out = ""
    for key in k:
        hash_out += "{} ".format(r.hashes[key])
    print(r.date, r.status, hash_out)

