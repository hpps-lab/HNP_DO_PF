#!/bin/sh
#$ -cwd

#./run > output.txt 2> errlog.txt
nohup ./run > output.txt 2> errlog.txt &