import os, sys, optparse
from das_client import get_data



template = '''
import CRABClient
from CRABClient.UserUtilities import config

config = config()

# All output/log files go in directory workArea/requestName/
config.General.workArea = 'crab_signal_projects'

config.General.requestName = 'crab_NLO_HToSSTo4l_MHMASSHIGGS_MSMASSX_ctauSCTAU_13TeV-NTuples'

config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

# Set pluginName = Analysis if you are reading a dataset, or to PrivateMC if not (so you are generating events)
config.JobType.pluginName = 'Analysis'
# CMSSW cfg file you wish to run
config.JobType.psetName = 'runLongLived2017_MCAnalysis_cfg.py'
# Increase virtual memory limit (sum needed by all threads) from default of 2000 MB.
config.JobType.maxMemoryMB = 2500
# Number of threads to use.
#config.JobType.numCores = 2
# To allow use of SL7 CMSSW versions that were SL6 at time of original DATA/MC production.
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['PUreweighting/2016DataPileupHistogram.root',
                             'PUreweighting/2016MCPileupHistogram.root']
config.JobType.outputFiles = ['output.root']

config.Data.inputDataset = 'INPUTDATASET'
config.Data.inputDBS = 'phys03'


# Units of "totalUnits" and "unitsPerJob" (e.g. files, events, lumi sections)
config.Data.splitting = 'FileBased'
# Total number of these units to be processed.
#config.Data.totalUnits = NTOTAL
# Requested number in each subjob.
config.Data.unitsPerJob = NSPLITTING
config.Data.outLFNDirBase = '/store/user/fernance/'
config.Data.publication = False
config.Data.outputDatasetTag = 'NTuples-Galapago'

config.Site.storageSite = 'T3_CH_CERNBOX'
'''    

def makeFile(higgsmass, xmass, ctau, inputdataset, nsplitting):


    text = template
    text = text.replace('MASSX', xmass)
    text = text.replace('CTAU', ctau)
    text = text.replace('MASSHIGGS', higgsmass)
    text = text.replace('INPUTDATASET', inputdataset)
    text = text.replace('NSPLITTING', nsplitting)


    name = 'crab_NLO_HToSSTo4l_MH{0}_MS{1}_ctauS{2}_13TeV-NTuple.py'.format(higgsmass, xmass, ctau)
    f = open(name, 'w')
    f.write(text + '\n')
    f.close()

    return name
 

if __name__=='__main__':

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    #parser.add_option('-i', '--input', action='store', type=str, dest='input', default='input', help='txt file with datasets')
    parser.add_option('-v', '--verbose', action='store_true', dest='verbose', help='verbose')
    parser.add_option('--submit', action='store_true', dest='submit', help='Submit the job if true')
    (opts, args) = parser.parse_args()

    mass_points = []
    #mass_points.append(['300', '150'])
    #mass_points.append(['300', '50'])
    #mass_points.append(['300', '20'])
    mass_points.append(['400', '150'])


    # Get datasets
    for mp in mass_points:
        print(">>> Mass point mH = {0} GeV and mS = {1} GeV".format(mp[0], mp[1]))
        dcommand = """/cvmfs/cms.cern.ch/common/dasgoclient --query="dataset=/ggH_HToSSTo4l_MH-{0}_MS-{1}_ctauS-*_Tune*_13TeV-powheg-pythia8*/*RunIIFall17MiniAODv2ext*/USER instance=prod/phys03 " """.format(mp[0], mp[1])
        dataset_outstring = os.popen(dcommand).read()
        datasets = dataset_outstring.split('\n')

        for n in range(0, len(datasets) - 1):
            dataset_name = datasets[n]
            command = """/cvmfs/cms.cern.ch/common/dasgoclient --query="dataset={0} instance=prod/phys03 | grep dataset.nfiles" -format string""".format(dataset_name)
            nfiles = os.popen(command).read()
            nfiles = int(nfiles)

            mH = (dataset_name.split('_MH')[1]).split('_')[0]
            mS = (dataset_name.split('_MS')[1]).split('_')[0]
            ctau = (dataset_name.split('_ctauS')[1]).split('_')[0]

            if opts.verbose:
                print('   > Dataset name:')
                print('   >   ' + dataset_name)
                print('   > Number of files:')
                print('   >   ' + str(nfiles))

            crab_file = makeFile(mH, mS, ctau, dataset_name, str(nfiles))

            print('   > crab file to launch:')
            print('   >   ' + crab_file)

            ### Launch job if submit mode activated:
            if opts.submit:
                os.system('crab submit ' + crab_file)




