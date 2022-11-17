#!/bin/bash 

echo -e "\n#manifest file for basic mapping assembly with illumina data using MIRA 4\n\nproject = ${1}\n\njob=genome,mapping,accurate\n\nparameters = -NW:mrnl=0 -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no\n\nreadgroup\nis_reference\ndata = ${2}\nstrain = ${3}\n\nreadgroup = reads\ndata = ${4}\ntechnology = solexa\nstrain = ${1}\n" > $5
