#!/bin/sh
#
#

# request bash shell as shell for job
#$ -S /bin/bash
#
# set working directory
cd ~/Thomas/hyperplan_server

# set local variables
outfile='hyperplane.out'.$$

# execute command(s)
R CMD BATCH run_server.R $outfile

# notify when script has finished
echo "Outfile is `pwd`/$outfile" | mail -v -s "hyperplane.sh finished" jiaqiyin@uw.edu
