#!/bin/sh

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
# https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#Pile_up 
PILEUP_LATEST=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt
JSON=$CMSSW_BASE/src/Samples/Tools/data/json/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
LUMI=41530

if [ ! -f "$PILEUP_LATEST" ]; then
    echo "File $PILEUP_LATEST does not exist on this site, copying from lxplus"
    if [[ $HOSTNAME == *"clip"* ]]; then
        echo "Running from CLIP, using CERN_USER variable $CERN_USER"
        export USR=$CERN_USER
    else
        export USR=$USER
    fi
    scp $USR@lxplus.cern.ch:$PILEUP_LATEST pileup_latest_2017.txt
    PILEUP_LATEST=pileup_latest_2017.txt
fi


echo "Calculating PU 2017 XSecVDown"
pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 62834 --maxPileupBin 200 --numPileupBins 200 PU_2017_${LUMI}_XSecVDown.root
echo "Calculating PU 2017 XSecDown"
pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 66017 --maxPileupBin 200 --numPileupBins 200 PU_2017_${LUMI}_XSecDown.root
echo "Calculating PU 2017 XSecCentral"
pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 69200 --maxPileupBin 200 --numPileupBins 200 PU_2017_${LUMI}_XSecCentral.root
echo "Calculating PU 2017 XSecUp"
pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 72383 --maxPileupBin 200 --numPileupBins 200 PU_2017_${LUMI}_XSecUp.root
echo "Calculating PU 2017 XSecVUp"
pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 75566 --maxPileupBin 200 --numPileupBins 200 PU_2017_${LUMI}_XSecVUp.root
echo "Calculating PU 2017 XSecVVUp"
pileupCalc.py -i $JSON --inputLumiJSON $PILEUP_LATEST --calcMode true --minBiasXsec 78750 --maxPileupBin 200 --numPileupBins 200 PU_2017_${LUMI}_XSecVVUp.root
