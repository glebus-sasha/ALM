#!/bin/bash
gapseq find -p all -b 200 -m Bacteria /opt/static/metaSPAdes_S3_Contigs.fasta
filename='/opt/static/tren.fasta'
gapseq find -p all -b 100 -m Bacteria $filename
gapseq find-transport -b 100 $filename
gapseq draft -r $shname-all-Reactions.tbl -t $shname-Transporter.tbl -p $shname-all-Pathways.tbl -c $filename -u 100 -l 50
gapseq fill -m $shname-draft.RDS -n /opt/gapseq/dat/media/meerwasser.csv -c $shname-rxnWeights.RDS -b 100 -g $shname-rxnXgenes.RDS

