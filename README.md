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
