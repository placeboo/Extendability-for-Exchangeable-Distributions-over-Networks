#!/bin/sh
#
#

# request bash shell as shell for job
#$ -S /bin/bash
#
# set working directory
cd ~/Thomas/monks_server

# set local variables
outfile='monk1_80.out'.$$

# execute command(s)
R CMD BATCH monk80.R $outfile

# notify when script has finished
echo "Outfile is `pwd`/$outfile" | mail -v -s "monk80.sh finished" jiaqiyin@uw.edu
