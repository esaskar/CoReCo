#!/bin/bash

set -e
wget -O gtg.zip https://www.dropbox.com/sh/sqlycf53kef1q8e/AACDSe3FihK79Gx775J1qjzEa/GTG?dl=1
unzip gtg.zip
rm gtg.zip
mv GTG/* .
rmdir GTG nrdb40 
