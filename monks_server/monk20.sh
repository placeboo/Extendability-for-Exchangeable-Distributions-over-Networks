#!/bin/sh
#
#

# request bash shell as shell for job
#$ -S /bin/bash
#
# set working directory
cd ~/Thomas/monks_server

# set local variables
outfile='monk1_20.out'.$$

# execute command(s)
R CMD BATCH monk20.R $outfile

# notify when script has finished
echo "Outfile is `pwd`/$outfile" | mail -v -s "monk1_20.sh finished" jiaqiyin@uw.edu
