import shutil

"""
IMPORTANT: This script should be run before running the analyzer
"""

####################################
##  PUreweighting initialization  ##
####################################
print('>> Accessing /eos/user/f/fernance/LLP_Analysis/...')

shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/2016DataPileupHistogram.root', 'test/2016/PUreweighting')
shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/2016MCPileupHistogram.root', 'test/2016/PUreweighting')

print('> PUreweighting histograms have been copied')

