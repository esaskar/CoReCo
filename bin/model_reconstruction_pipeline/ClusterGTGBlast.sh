#$ -S /bin/bash
#$ -N CLusterGTGBlast
#$ -j y
#$ -cwd
#$ -pe smp 4
#$ -R y

######################################################################  blastp  ##############################################################################

if [ -z "$1" ]; then
   echo "Specify blastp path";
   exit;

else
   blastp=$1
fi

###################################################################### blastDB ##############################################################################


if [ -z "$2" ]; then
   echo "Specify blast database";
   exit;

else
   blastDB=$2
fi

###################################################################### queryFastaFile ##############################################################################

if [ -z "$3" ]; then
   echo "Specify fasta file";
   exit;

else
   queryFastaFile=$3 
fi

##################################################################### outfmt ###############################################################################

if [ -z "$4" ]; then
   echo "Specify output format";
   exit;
else
   outfmt=$4
fi

###################################################################### blastResultFile ##############################################################################


if [ -z "$5" ]; then
   echo "Specify output file";
   exit;
else
   blastResultFile=$5
fi

####################################################################### eValue #############################################################################

if [ -z "$6" ]; then
   echo "Specify output file";
   exit;
else
   eValue=$6
fi

####################################################################### numThreads #############################################################################

#if [ -z "$7" ]; then
#   echo "Specify number of threads";
#   exit;
#else
#   numThreads=$7
#fi

###################################################################### databaseSize #############################################################################



if [ -z "$7" ]; then
	$blastp -db $blastDB -query $queryFastaFile -outfmt $outfmt -out $blastResultFile -evalue  $eValue -num_threads 4 
else
   uniprotDatabaseSize=$7
   $blastp -db $blastDB -query $queryFastaFile -outfmt $outfmt -out $blastResultFile -evalue $eValue -dbsize $uniprotDatabaseSize -num_threads 4
   
fi

####################################################################################################################################################




