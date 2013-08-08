#!/bin/bash

# get dust information for each kind of fungi
echo "1---> begin to dust step"
bash makedustfile_AUTO.bash
echo "finish dust step"

sleep 2

echo "2--->begin to makedb step"
# make database for each fungi species for backword blast preparation
bash makeblastdb_AUTO.bash
echo "finish makedb step"

sleep 2

echo "3--->begin to forward blast step"
# forward blast: species2uniprot
bash blastp_fwdblast.bash
echo "finish forward blast step"

sleep 2

echo "4--->begin the backward blast step"
# backward blast: uniprot2speices
bash blastp_bwdblast.bash
echo "FINNISH ALL"



