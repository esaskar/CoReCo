#!/bin/bash
#
# Required:
# $1  target dir :result dir
# 
# Optional
# $2  EC score dir : .ecscore
# $3  CPD file :all.blastpos/all.blastneg/all.gtgpos/all.gtgneg

if [ -z "$1" ]; then
   echo "Specify data directory";
   echo "usage: $0 data-dir [ec-score-dir] [model-file]"
   exit;
fi

SOURCE_DIR=/mnt/nfs/SynBioModelling/ModelTraining

TARGET_DIR=$1

if [ $2 ]; then
    EC_SCORE_DIR=$SOURCE_DIR/$2
else
    EC_SCORE_DIR=$SOURCE_DIR/rawscore
fi

if [ $3 ]; then
    CPD_FILE=$SOURCE_DIR/$3
else
    CPD_FILE=$SOURCE_DIR/cpds_1/all
fi

ORG_LIST=$SOURCE_DIR/org_list.backup
TREE_FILE=$SOURCE_DIR/tree.txt

IPRSCAN_MODEL_DIR=$SOURCE_DIR/50fungi_modelecs

#YEASTNET_ECS=/home/esa/projects/50fungi/reconstruction/public-models/yeastnet/eclist
#ANIGER_ECS=/home/esa/projects/50fungi/reconstruction/public-models/aniger_nielsen/eclist

mkdir -p $TARGET_DIR

echo Importing data to $TARGET_DIR...
python ../model_training_scripts/import_data.py $ORG_LIST $EC_SCORE_DIR $CPD_FILE $TREE_FILE $TARGET_DIR
if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi

echo Copying default parameters...
cp data/parameters $TARGET_DIR/
if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi

# Write reaction-scores
echo Computing reaction scores...
python compute_reaction_scores.py $TARGET_DIR
if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi

# Write reaction-scores.full
echo Merging posterior, BLAST and GTG scores...
python merge_scores.py $TARGET_DIR $EC_SCORE_DIR
if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi

echo Plotting reaction scores...
python plot_reaction_scores.py $TARGET_DIR
if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi

echo Plotting ROC curves/IPR models...
python roc.py $TARGET_DIR $IPRSCAN_MODEL_DIR
if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi

echo Plotting ROC curves/Yeast consensus model...
python roc.py $TARGET_DIR $YEASTNET_ECS Scer rocs-yeastnet
if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi

echo Plotting ROC curves/A. niger model...
python roc.py $TARGET_DIR $ANIGER_ECS Anig rocs-aniger
if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi

echo Performing yeast model comparison...
python compare_yeast_prediction.py $TARGET_DIR $YEASTNET_ECS $IPRSCAN_MODEL_DIR/Scer $TARGET_DIR/yeast-comparison-scores.txt
if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi

echo Drawing EC score trees...
python visualize_ec_scores.py $TARGET_DIR
if [ $? -ne 0 ]; then echo "Failed"; exit 1; fi
