#!/bin/bash
for filename in /opt/static/data/*; 
do 
path="/opt/static/data/" 
path2="/opt/static/results/" 
xbase=${filename##*/}
shname=$path2${xbase%.*}
gapseq find -p all -b 100 -m Bacteria $filename
gapseq find-transport -b 100 $filename
gapseq draft -r $shname-all-Reactions.tbl -t $shname-Transporter.tbl -p $shname-all-Pathways.tbl -c $filename -u 100 -l 50
gapseq fill -m $shname-draft.RDS -n /opt/gapseq/dat/media/meerwasser.csv -c $shname-rxnWeights.RDS -b 100 -g $shname-rxnXgenes.RDS
done

