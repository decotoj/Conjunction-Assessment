# Conjunction Assessment Demo
# Jake Decoto (decotoj@gmail.com)
# Last Update: October 4, 2022
# Original Version: November 10, 2021

from datetime import datetime
import pandas
from collections import defaultdict
import time

# try:
#     import catapult.util_conjunction_assessment as ca # private party library
# except:
import util_conjunction_assessment as ca # standalone version for demo (not in catapult library)

# All vs All Satellite Conjunction Assessment Demo

# NOTE: 10/4/2022, All vs All 24 hour, 1 second step size conjunction screening with ~20,000 satellites takes
# approximately 45 minutes on reference desktop class machine

if __name__ == "__main__":

    # config
    cat = '20220918_TLE.txt' # filepath to TLE catalog
    dt0 = datetime.strptime('2022-09-18 00:00:00.000', '%Y-%m-%d %H:%M:%S.%f').replace(tzinfo=None) # start epoch
    steps = 86400 # number of steps
    stepSize = 1 # step size (seconds)
    out_file = 'conjunction_report.csv'
    maxEphemBlockSize = 3600 # maximum number of steps processed at once, default = 3600, may run into memory issues if set too high
    rt = 10 # range threshold (km)

    # Perform Conjunction Run
    t0 = time.time()
    conj = ca.conjunction_run(cat, rt, dt0, steps, stepSize, maxEphemBlockSize)
    print(f'Total Conjunction Run Time {time.time()-t0}')

    # Reformat Output
    conjunctions = defaultdict(list)
    conjunctions['Epoch'] = [conj[k]['Epoch'] for k in conj.keys()]
    conjunctions['Object A'] = [k[0] for k in conj.keys()]
    conjunctions['Object B'] = [k[1] for k in conj.keys()]
    conjunctions['Range (km)'] = [conj[k]['Range (km)'] for k in conj.keys()]

    # Output CSV File for Humans
    out = pandas.DataFrame.from_dict(conjunctions)
    out.to_csv(out_file)
    print(f'Report Produced: {out_file}')
