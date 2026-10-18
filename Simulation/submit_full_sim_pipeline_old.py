#!/usr/bin/env python3

from subprocess import Popen, PIPE
import sys

if len(sys.argv) != 2:
	print('Usage: ./submit_full_sim_pipeline.py <number_of_jobs>')
	sys.exit(1)

NJOBS = int(sys.argv[1])

SCRIPT = '/work/clas12/users/nwright/rgm_andrew/Simulation/submit_GCF_GEMC_RECON.sh'  # path to your bash script

for idx in range(1, NJOBS + 1):

	command="""#!/bin/sh
	#SBATCH --job-name=rgm_{1}
	#SBATCH --account=clas12
	#SBATCH -p production
	#SBATCH --constraint=farm19
	#SBATCH --mem-per-cpu=5000
	#SBATCH -t240
	#SBATCH --error=/farm_out/nwright/err_rgm_{1}.txt
	#SBATCH --output=/farm_out/nwright/out_rgm_{1}.txt

	time {0} {1}
	""".format(SCRIPT, idx)
	command = command.replace('\t', '')
	print(command)

	p = Popen(args=["sbatch"], stdin=PIPE)
	p.communicate(command.encode())