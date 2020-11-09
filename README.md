# IFCALongLivedAnalysis

This framework is devoted to get the relevant collections from MINIAOD and put them into plain trees in the context of the Long Lived analysis developed at IFCA.

Part of the code works with muon collections (reco::Track objects) that are only available in AOD format. This part of the analyzer must be run in modified MINIAOD formats that have these collections added.

## Installation

This analyzer works with 94X:

```
cmsrel CMSSW_9_4_4

cd CMSSW_9_4_4/src

cmsenv

mkdir MyAnalysis

cd MyAnalysis

git clone https://github.com/longlivedpeople/IFCALongLivedAnalysis.git
```

In order to compute the lumiweights, the histograms with the PU distributions must be copied in the working directory:

```
python init.py
```

Please, make sure that you have reading permissions to /eos/user/f/fernance/.


## PU reweighting

The configurations for computing the PU weights stored in the MC and Data pileup distributions are

**2016**
MC PU scenario: SimGeneral.MixingModule.mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi
GOLDEN Json: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt
PU Json: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt

Command to get the DATA distribution (with ```CMSSW_9_4_4```)
```
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100  2016DataPileupHistogram.root
```


**2017**
MC PU scenario: SimGeneral.MixingModule.mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi


**2018**
MC PU scenario: SimGeneral.MixingModule.mix_2018_25ns_ProjectedPileup_PoissonOOTPU_cfi
GOLDEN Json: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
PU Json: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt

Command to get the DATA distribution (with ```CMSSW_10_2_5```)
```
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100  2018DataPileupHistogram.root
```




