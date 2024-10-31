#!/bin/bash
# Author: Charlotte Vavourakis

# stool
tr -d '\r' < archaeal_ASV.txt > archaeal_ASV_unix.txt
out1=$(paste -s -d ' ' archaeal_ASV_unix.txt)
singularity exec ./mafft.sif mafft --auto stool_S_ctrl_ASV.fasta > stool_S_ctrl_ASV.aligned.fasta
singularity exec ./fasttree.sif FastTree -nt -gtr -gamma -spr 4 -boot 100 stool_S_ctrl_ASV.aligned.fasta > stool_S_ctrl_ASV_treeML.nwk
#echo "singularity exec ./nw-utils.sif nw_reroot stool_S_ctrl_ASV_treeML.nwk $out1 > stool_S_ctrl_ASV_rootedarchaea_treeML.nwk"
singularity exec ./nw-utils.sif nw_reroot stool_S_ctrl_ASV_treeML.nwk $out1 > stool_S_ctrl_ASV_rootedarchaea_treeML.nwk

# re-root, alternative:
#command="singularity exec ./nw-utils.sif nw_reroot stool_S_ctrl_ASV_treeML.nwk"
#
#while IFS= read -r line; do
#    command="$command $line"
#done < archaeal_ASV_unix.txt
#
#command="$command > stool_S_ctrl_ASV_rootedarchaea_treeML.nwk"
#
#echo "$command"
#eval "$command"

# saliva
tr -d '\r' < cpr_ASV.txt > cpr_ASV_unix.txt
out2=$(paste -s -d ' ' cpr_ASV_unix.txt)
singularity exec ./mafft.sif mafft --auto saliva_S_ctrl_ASV.fasta > saliva_S_ctrl_ASV.aligned.fasta
singularity exec ./fasttree.sif FastTree -nt -gtr -gamma -spr 4 -boot 100 saliva_S_ctrl_ASV.aligned.fasta > saliva_S_ctrl_ASV_treeML.nwk
#echo "singularity exec ./nw-utils.sif nw_reroot saliva_S_ctrl_ASV_treeML.nwk $out2 > saliva_S_ctrl_ASV_rootedarchaea_treeML.nwk"
singularity exec ./nw-utils.sif nw_reroot saliva_S_ctrl_ASV_treeML.nwk $out2 > saliva_S_ctrl_ASV_rootedcpr_treeML.nwk
