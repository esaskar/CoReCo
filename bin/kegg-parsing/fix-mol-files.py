#!/usr/bin/python  

##
## This script was previously used like this:
##
##    echo "  Download latest version 'mol files 28102013.zip' from "
##    echo "  http://myworks.vtt.fi/projektit/tk40/eu-bioledge/Documents/lists of reactions from Peter VERSION 4/"
##    echo "Fixing file names to all have capital letter in file name"
##    rename 's/cluster/Cluster/' $KDIR/mol/cluster*.mol
##
##    echo "Fixing mol files with R# or similar atoms"
##    mkdir -p $KDIR/mol/original_mol_before_fixing/
##    
##    FILES=$KDIR/mol/*.mol
##
##    for f in $FILES
##    do
##	echo "Processing $f file..."
##	cp $f $f.original
##	./fix-mol-files.py $f.original $f > $DIR/molfixing.log 
##	FILESIZE=$(stat -c%s "$DIR/molfixing.log")
##	if [ $FILESIZE == "0" ]; then
##	    echo "No change in the file: $f"
##	    diff $f $f.original
##	    rm $f.original
##	else
##	    echo "file changed: $f"
##	    diff $f $f.original	    
##	fi
##    done
##    mv $KDIR/mol/*.original $KDIR/mol/original_mol_before_fixing/


import sys, re          

f = open(sys.argv[1])      
outstring=""

changed = False

# turhat header-rivit
section = "HEADER"
outstring += f.readline()
outstring += f.readline()
outstring += f.readline()
outstring += f.readline()


reR = re.compile("R[1-9]+") # regulation expression patterns   					     for reEC; \s means whitespace   					     character;+

for line in f:
    newline = ""
    if line.startswith("M") or line.startswith("A") or line.startswith("R") :
            ## added "A" and "R" to handle a mol file with R#'s (biocyc 0301121217)
            ##A   65
            ##R2
            ##A   64
            ##R1
            section = "OTHER"
    elif len(line.strip().split()) in [6,7]:
            section = "BONDS"
    elif len(line.strip().split()) in [11,13,16]: ## Added 11 as additional option to allow parsing of biocyc .mol-files (Merja Oja, 29.10.2013)
            section = "ATOMS"

    if section == "ATOMS":
            words = line.strip().split()

            symbol = words[3]
            x = float(words[0])
            y = float(words[1])
            
            if symbol == ".":
                newline=re.sub(" \. "," * ",line)
                changed = True
                print ". symbol encountered ; replacing with *"
            elif symbol == "D":
                newline=re.sub(symbol,"H",line)
                changed = True
                print "D symbol encountered ; replacing with H"
            elif symbol.upper() == "X" or symbol == "A" or symbol == "0":
                newline=re.sub(symbol,"R",line)
                changed = "True"
                print ("%s symbol encountered ; replacing with R" % symbol)
            elif symbol == "R#":
                newline=re.sub(symbol,"R ",line)
                changed = True
                print "R# symbol encountered ; replacing with R"
            elif reR.match(symbol):
                newline=re.sub("R[1-9]+","R ",line)
                changed = True
                print ("%s symbol encountered ; replacing with R" % symbol)

            if newline != "":
                print "oldline: " + line
                print "newline: " + newline
                outstring += newline
            else:
                outstring += line


    else:
        outstring += line
        

f.close()

if changed:
    print "Printing new .mol file"
    o = open(sys.argv[2],"w")  
    o.write(outstring)
    o.close()
    

