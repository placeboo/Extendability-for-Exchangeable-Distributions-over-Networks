#!/bin/sh
#
#

# request bash shell as shell for job
#$ -S /bin/bash
#
# set working directory
cd ~/Thomas/monks_server

# set local variables
outfile='monk.out'.$$

# execute command(s)
R CMD BATCH monk.R $outfile

# notify when script has finished
echo "Outfile is `pwd`/$outfile" | mail -v -s "monk.sh finished" jiaqiyin@uw.edu
