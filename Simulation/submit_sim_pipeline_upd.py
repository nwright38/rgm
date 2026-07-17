#!/usr/bin/env python3

from subprocess import Popen, PIPE
import sys

if len(sys.argv) != 2:
	print('Usage: ./submit_full_sim_pipeline.py <number_of_jobs>')
	sys.exit(1)

NJOBS = int(sys.argv[1])
NEVENTS_PER_JOB = 100000 # go for about 10% of the total events

SCRIPT = '/work/clas12/users/nwright/rgm_andrew/Simulation/submit_sim_pipeline_upd.sh'
START_IDX = 4001
for idx in range(START_IDX, NJOBS + START_IDX):

	command = f"""#!/bin/bash
#SBATCH --job-name=rgm_{idx}
#SBATCH --account=clas12
#SBATCH -p production
#SBATCH --constraint=farm19
#SBATCH --mem-per-cpu=10000
#SBATCH --time=13:00:00
#SBATCH --error=/farm_out/nwright/err_rgm_{idx}.txt
#SBATCH --output=/farm_out/nwright/out_rgm_{idx}.txt

time {SCRIPT} {idx} {NEVENTS_PER_JOB}
"""

	print(command)

	p = Popen(args=["sbatch"], stdin=PIPE)
	p.communicate(command.encode())