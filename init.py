import shutil
import os

"""
IMPORTANT: This script should be run before running the analyzer
"""

####################################
##  PUreweighting initialization  ##
####################################
print('>> Accessing /eos/user/f/fernance/LLP_Analysis/...')

shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/2016DataPileupHistogram.root', 'test/2016/PUreweighting')
shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/2016DataPileupHistogram.root', './')
shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/2016MCPileupHistogram.root', './')

print('> PUreweighting histograms have been copied')

print('>> Creating the logs folder in test...')

if not os.path.exists('test/logs/'): 
    os.makedirs('test/logs/')
    print('> test/logs/ folder created')
else:
    print('> The folder already exists: Passing...')
