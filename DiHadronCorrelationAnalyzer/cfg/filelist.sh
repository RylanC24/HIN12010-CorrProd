#!/bin/sh
PATH_TO_FILES=$1

rm filelist_noGplus_tmp.txt
rm filelist_noGplus.txt
rm filelist_Gplus_tmp.txt
rm filelist_Gplus.txt

find /mnt/hadoop/cms/store/user/davidlw/PPJet/PP2013_FlowCorr_PromptReco_HighPtTrk_NoGplus_v2/2a01db98f9c5d38ba9c3802e841ac721/*.root >> filelist_noGplus_tmp.txt
find /mnt/hadoop/cms/store/user/davidlw/PPJet/PP2013_FlowCorr_PromptReco_HighPtTrk_Gplus_v2/73d1f7d214283b4725d2fcb6e434a68d/*.root >> filelist_Gplus_tmp.txt

sed -e 's/\/mnt\/hadoop\/cms\/store/\/store/g' filelist_noGplus_tmp.txt >> filelist_noGplus.txt 
sed -e 's/\/mnt\/hadoop\/cms\/store/\/store/g' filelist_Gplus_tmp.txt >> filelist_Gplus.txt 
