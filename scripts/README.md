# Create PU reweighting histograms

## 2016 

Note: Both pre-VFP and post-VFP have the same premixing scenario. 

### Data distribution

Golden json: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
PileUp json: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt

```
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 2016DataPileupHistogram.root
```

### Monte Carlo

The premixing library is /Neutrino_E-10_gun/RunIISummer20ULPrePremix-UL16_106X_mcRun2_asymptotic_v13-v1/PREMIX and the PU profile is found in https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/SimGeneral/MixingModule/python/mix_2016_25ns_UltraLegacy_PoissonOOTPU_cfi.py

```
python makeMCPileupHist.py SimGeneral.MixingModule.mix_2016_25ns_UltraLegacy_PoissonOOTPU_cfi --outputFilename 2016MCPileupHistogram.root
```





