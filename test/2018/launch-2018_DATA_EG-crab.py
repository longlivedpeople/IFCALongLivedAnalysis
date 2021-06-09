import os, sys, optparse
from subprocess import Popen


template = '''
import CRABClient
from CRABClient.UserUtilities import config

config = config()

# All output/log files go in directory workArea/requestName/
config.General.workArea = 'crab_projects_TASKTAG'

config.General.requestName = 'crab_JOBTAG-NTuples'

config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

# Set pluginName = Analysis if you are reading a dataset, or to PrivateMC if not (so you are generating events)
config.JobType.pluginName = 'Analysis'
# CMSSW cfg file you wish to run
config.JobType.psetName = 'runLongLived2018_DATA_EG_cfg.py'
# Increase virtual memory limit (sum needed by all threads) from default of 2000 MB.
config.JobType.maxMemoryMB = 2500
# Number of threads to use.
#config.JobType.numCores = 2
# To allow use of SL7 CMSSW versions that were SL6 at time of original DATA/MC production.
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['PUreweighting/2018DataPileupHistogram.root',
                             'PUreweighting/2018MCPileupHistogram.root']
config.JobType.outputFiles = ['output.root']

config.Data.inputDataset = 'INPUTDATASET'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
config.Data.inputDBS = 'global'


# Units of "totalUnits" and "unitsPerJob" (e.g. files, events, lumi sections)
config.Data.splitting = 'FileBased'
# Total number of these units to be processed.
#config.Data.totalUnits = NTOTAL
# Requested number in each subjob.
config.Data.unitsPerJob = NSPLITTING
config.Data.outLFNDirBase = '/store/user/fernance/'
config.Data.publication = False
config.Data.outputDatasetTag = 'JOBTAG_2018_NTuples-Galapago'


config.Site.storageSite = 'T2_ES_IFCA'
'''    

def makeFile(inputdataset, nsplitting, jobtag, tasktag):


    text = template
    text = text.replace('INPUTDATASET', inputdataset)
    text = text.replace('NSPLITTING', nsplitting)
    text = text.replace('JOBTAG', jobtag)
    text = text.replace('TASKTAG', tasktag)

    #dir_ = 'crab_dir_' + tasktag + '/'
    dir_ = ''
    if not os.path.exists(dir_) and dir_ != '': os.makedirs(dir_) 
    name_ = 'crab_{0}-NTuple.py'.format(jobtag)
    f = open(dir_ + name_, 'w')
    f.write(text + '\n')
    f.close()

    return dir_,name_
 

if __name__=='__main__':

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-t', '--tasktag', action='store', type=str, dest='tasktag', default='jobtag', help='tasktag')
    parser.add_option('-v', '--verbose', action='store_true', dest='verbose', help='verbose')
    parser.add_option('--submit', action='store_true', dest='submit', help='Submit the job if true')
    (opts, args) = parser.parse_args()

    # List of datasets 

    Datasets = {}
    Datasets['Run2018A']  = '/EGamma/Run2018A-12Nov2019_UL2018-v2/MINIAOD'
    Datasets['Run2018B']  = '/EGamma/Run2018B-12Nov2019_UL2018-v2/MINIAOD'
    Datasets['Run2018C']  = '/EGamma/Run2018C-12Nov2019_UL2018-v2/MINIAOD'
    Datasets['Run2018D']  = '/EGamma/Run2018D-12Nov2019_UL2018-v6/MINIAOD'

    # Get datasets
    for key in Datasets.keys():
        dataset = Datasets[key]
        print(">>> Dataset: " + dataset)

        command = """/cvmfs/cms.cern.ch/common/dasgoclient --query="dataset={0} instance=prod/global | grep dataset.nfiles" -format string""".format(dataset)
        nfiles = os.popen(command).read()
        nfiles = int(nfiles)

        if opts.verbose:
            print('   > Dataset name:')
            print('   >   ' + dataset)
            print('   > Number of files:')
            print('   >   ' + str(nfiles))

        split = int(nfiles / 500)
        if split == 0: split = 1

        crab_dir,crab_file = makeFile(dataset, str(split), key, opts.tasktag)

        print('   > crab file to launch:')
        print('   >   ' + crab_dir + crab_file)

        ### Launch job if submit mode activated:
        if opts.submit:
            command = 'crab submit ' + crab_file
            os.system(command)
            os.system('rm ' + crab_file)
            #Popen(command, shell = True, cwd = crab_dir)




