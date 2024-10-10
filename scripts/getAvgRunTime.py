#!/usr/bin/python3
import htcondor2 as htcondor
import numpy as np
import sys
import matplotlib.pyplot as plt

schedd = htcondor.Schedd()

def seconds_to_hours(seconds):
  return seconds / 3600

def main():
  savedir = "/data/alice/gweelden/Fragmentation/pythia/SubmitScripts/plots/"
  savename = "gweelden"
  myConstraint = 'Owner == "gweelden"'
  if (len(sys.argv) > 1):
    myConstraint = 'ClusterId == ' + sys.argv[1]
    savename = sys.argv[1]

  print("Getting history with constraint:", myConstraint)
  history = schedd.history(
    constraint=myConstraint,
    projection=['ClusterId', 'ProcId', 'JobStatus', 'RemoteWallClockTime']
    )

  # pairs = [] # ['ProcId', 'RemoteWallClockTime']
  # proc_ids = []
  walltimes_s = []
  for job in history:
    if (job['JobStatus'] != 4): # Skip if job is not done
      continue
    walltimes_s.append(job['RemoteWallClockTime'])
    # proc_ids.append(job['ProcId'])
    # pair = [job['ProcId'], job['RemoteWallClockTime']]
    # pairs.append(pair)

  walltimes_h = [seconds_to_hours(w) for w in walltimes_s]

  min_walltime = round(min(walltimes_h), 1)
  avg_walltime = round(np.average(walltimes_h), 1)
  max_walltime = round(max(walltimes_h), 1)
  print(min_walltime, "(min)", avg_walltime, "(avg)", max_walltime, "(max)")

  var_walltime = round(np.var(walltimes_h), 2)
  stddev_walltime = round(np.std(walltimes_h), 2)
  print(var_walltime, "(var)", stddev_walltime, "(stddev)")

  five_walltime = round(np.quantile(walltimes_h, 0.05), 1)
  ten_walltime = round(np.quantile(walltimes_h, 0.1), 1)
  ninety_walltime = round(np.quantile(walltimes_h, 0.9), 1)
  ninetyfive_walltime = round(np.quantile(walltimes_h, 0.95), 1)
  print(five_walltime, "(5%)", ten_walltime, "(10%)", ninety_walltime, "(90%)", ninetyfive_walltime, "(95%)")

  savename = savedir + "walltimes_" + savename + ".png"
  print("Saving plot to", savename)
  plt.hist(walltimes_h, density=False, bins=2*int(np.ceil(max_walltime)))
  plt.title("Walltimes (" + myConstraint + "); " + str(avg_walltime) + u"\u00B1" + str(var_walltime) + "h")
  plt.xlabel("Walltime (h)")
  plt.ylabel("N jobs")
  plt.savefig(savename)


if __name__ == "__main__":
    main()