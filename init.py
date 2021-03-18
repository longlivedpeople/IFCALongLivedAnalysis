import shutil
import os

"""
IMPORTANT: This script should be run before running the analyzer
"""

if not os.path.exists('test/2016/'): 
    os.makedirs('test/2016/')
    print('> test/2016/ folder created')
if not os.path.exists('test/2017/'): 
    os.makedirs('test/2017/')
    print('> test/2017/ folder created')
if not os.path.exists('test/2018/'): 
    os.makedirs('test/2018/')
    os.makedirs('test/2018/PUreweighting/')
    print('> test/2018/ folder created')

####################################
##  PUreweighting initialization  ##
####################################
print('>> Accessing /eos/user/f/fernance/LLP_Analysis/...')

print('>> Copying 2016 PU histograms...')
shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/Final/2016/2016DataPileupHistogram.root', 'test/2016/PUreweighting/')
shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/Final/2016/2016MCPileupHistogram.root', 'test/2016/PUreweighting/')
shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/Final/2016/2016DataPileupHistogram.root', './')
shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/Final/2016/2016MCPileupHistogram.root', './')

print('>> Copying 2018 PU histograms...')
shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/Final/2018/2018DataPileupHistogram.root', 'test/2018/PUreweighting/')
shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/Final/2018/2018MCPileupHistogram.root', 'test/2018/PUreweighting/')
shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/Final/2018/2018DataPileupHistogram.root', './')
shutil.copy('/eos/user/f/fernance/LLP_Analysis/PUreweighting/Final/2018/2018MCPileupHistogram.root', './')

print('> PUreweighting histograms have been copied')

print('>> Creating the logs folder in test...')

if not os.path.exists('test/logs/'): 
    os.makedirs('test/logs/')
    print('> test/logs/ folder created')
else:
    print('> The folder already exists: Passing...')
