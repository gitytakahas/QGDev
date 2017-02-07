#!/bin/bash

##python sendOnBatch.py \
##	--query \
##	-e /QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISummer16MiniAODv2-PUFlat0to70_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v4-v1/MINIAODSIM \
##	--put-in=/store/user/amarini/qg/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8 \
##	--follow \
##	-i qgMiniTupleForMiniAOD_cfg.py \
##	-d mysub -q 8nh -n 100
##


dirs=""
##giorgia
if [ "$USER" == "grauco" ];
then
	dirs+=" QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8"
fi

if [ "$USER" == "amarini" ];
then
#amarini
	dirs+=" QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8"
	dirs+=" QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8"
fi

DESTUSER=grauco
dest="/eos/user/${DESTUSER:0:1}/$DESTUSER/qg"
#eos root://eosuser mkdir /eos/cms/store/user/$DESTUSER/qg
	#--put-in=/store/user/$DESTUSER/qg/$dir \

for dir in $dirs; 
do
   #eos mkdir /eos/cms/store/user/$DESTUSER/qg/$dir
   mkdir /eos/user/${DESTUSER:0:1}/$DESTUSER/qg/$dir

   tags=$(eos ls /store/mc/RunIISummer16MiniAODv2/$dir/MINIAODSIM | grep PUMoriond17)

   for tag in $tags;
   do
   echo "-> donig $dirs/$tag: eos= /store/mc/RunIISummer16MiniAODv2/$dir/MINIAODSIM/$tag"
   python sendOnBatch.py \
	-e /store/mc/RunIISummer16MiniAODv2/$dir/MINIAODSIM/$tag \
	--instance root://eosuser \
	--put-in=/eos/user/${DESTUSER:0:1}/$DESTUSER/qg/$dir/$tag \
	--follow \
	-i ./qgMiniTupleForMiniAOD_cfg.py \
	-d mysub/$dir/$tag -q 8nh -n 25
   done
done
