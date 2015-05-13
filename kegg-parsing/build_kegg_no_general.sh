#!/bin/bash

echo "Building a version of KEGG with no general reactions"

RUNDIR=`pwd`
echo "Running the script from directory $RUNDIR"

SCRIPT=$(readlink -f $0)
SDIR=$(dirname $SCRIPT)

DIR=$RUNDIR/data/kegg/kegg-no-general 
KDIR=$RUNDIR/data/kegg #original files from KEGG
ADIR=$RUNDIR/data/kegg/aux 
MOLDIR=$KDIR/mol/ 

echo "Source: $KDIR"
echo "Target: $DIR"

mkdir -p $DIR
mkdir -p $ADIR

echo ""
if [ -e $DIR/reaction ]; then
    echo "Already created files enzyme, reaction, compound, and general-reactions to"
    echo "folder $DIR"
else
    echo "Filtering KEGG enzymes and reactions..."
    python $SDIR/filter_general_kegg_reactions.py $KDIR $DIR
    if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi
fi
    
echo ""
if [ -e $DIR/ec-list.txt ]; then
    echo "Already created ec-list.txt to"
    echo "folder $DIR"
else
    echo "Building EC to reaction map..."
    python $SDIR/extract-kegg-ecs-orthology.py $DIR/reaction $DIR/ec-list.txt
    if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi
fi    

echo ""
if [ -e $ADIR/reaction-pathways ]; then
    echo "Already created files ec-to-pathways pathway-names and reaction-pathways to"
    echo "folder $ADIR"
else
    echo "Extracting tmp1 pathway-names and tmp2 from KEGG"
    python $SDIR/extract-kegg-pathways.py $KDIR $ADIR/ec-to-pathways $ADIR/pathway-names $ADIR/reaction-pathways
    if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi
fi    

echo ""
if [ -e $ADIR/kegg-compounds ]; then
    echo "Already created kegg-compounds to"
    echo "folder $ADIR"
else
    echo "Extracting KEGG compound names..."
    python $SDIR/parse_compound_names.py $DIR/compound $ADIR/kegg-compounds
    if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi
fi    


echo ""
if [ -e $ADIR/kegg-metabolite-formulae ]; then
    echo "Already created kegg-metabolite-formulae to"
    echo "folder $ADIR"
else
    echo "Filtering kegg-metabolite-formulae..."
    python $SDIR/parse_kegg_compound_formulae.py $MOLDIR $ADIR/kegg-metabolite-formulae
    if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi
fi


echo ""
if [ -e $ADIR/cofactors ]; then
    echo "Already created cofactors to"
    echo "folder $ADIR"
else
    echo "creating cofactors empty file ..."
    echo "" > $ADIR/cofactors
    if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi
fi



echo ""
if [ -e $DIR/nooxygen.eqn ]; then
    echo "Already created files nooxygen.egn and nooxygen.status to"
    echo "folder $DIR"
else
    echo "Balancing KEGG reactions..."
    echo python $SDIR/check_element_balances.py -k $DIR -o $DIR/nooxygen
    python $SDIR/check_element_balances.py -k $DIR -o $DIR/nooxygen 
    if [ $? -ne 0 ]; then echo "Failed."; exit 1; fi
fi    

echo ""
if [ -e $DIR/nooxygen.mapperinput ]; then
    echo "Already created nooxygen.mapperinput to"
    echo "folder $DIR"
else
    echo "Converting stoichiometry to atom mapper input format..."
    python $SDIR/convert_balances_to_mapper.py $DIR/nooxygen.eqn $DIR/nooxygen.mapperinput
    if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi
fi

echo ""
if [ -e $DIR/nooxygen.mapperinput.metabolites ]; then
    echo "Already created nooxygen.mapperinput.metabolites to"
    echo "folder $DIR"
else
    echo "Extracting metabolites seen in the reactions"
    #
    # Extract the reaction part (columns 3 to end)
    # replace " p " with " "
    # replace " " with line break
    # keep each metabolite only once
    cut -f3- -d" " $DIR/nooxygen.mapperinput | perl -pe "s/ p / /g;" | perl -pe "s/ /\n/g;" | sort -u > $DIR/nooxygen.mapperinput.metabolites
fi

echo ""
if [ -e $KDIR/atommaps ]; then
    echo "Atom maps already computed copying"

    echo "Deleting atom maps of general reactions..."
    rm -f `cat $DIR/general-reactions | sed s/$/.txt/ | sed s/^/$DIR'\/atommaps\//'`

else

    ATOMFEATURES=$DIR/atomfeatures/
    mkdir -p $ATOMFEATURES
    
    ATOMMAPS=$DIR/atommaps/
    mkdir -p $ATOMMAPS
    
    KEGG_REACTIONS=$DIR/nooxygen.mapperinput
    KEGG_METABOLITES=$DIR/nooxygen.mapperinput.metabolites

    echo ""
    echo "Generating Atom Features.... This will take a while"
    echo "To track the progress check atomfeatures.log"
    echo "in $DIR".

    python atomfeatures/FeatureGenerator.py --inputfile $KEGG_METABOLITES --moldir $MOLDIR --output-dir $ATOMFEATURES -A -W -M -R --context 0-10 > $DIR/atomfeatures.log
    if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi
    
    echo ""
    echo "Generating Atom Maps.... This will take a long time (unless parallelized)"
    echo "To track the progress check atommaps.log"
    echo "in $DIR"
    cd atommapper2011/bin/
    java Mapper2000 -reacfile $KEGG_REACTIONS --moldir $MOLDIR -featdir $ATOMFEATURES -output $ATOMMAPS > $DIR/atommaps.log
    if [ $? -ne 0 ]; then 
	cd $RUNDIR
	echo "Failed"; 
	exit 1; 
    else 
	cd $RUNDIR
    fi
fi



echo "A version of KEGG without general reactions" >$DIR/README
echo "Built from \"$KDIR\" by $0 on `date`" >>$DIR/README
