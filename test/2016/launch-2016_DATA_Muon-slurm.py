import os, sys, optparse
from subprocess import Popen

runningfile = os.path.abspath(__file__)

WORKPATH = ''
for d in runningfile.split('/'):
    WORKPATH += d
    WORKPATH += '/'
    if d == 'IFCALongLivedAnalysis': break

CMSSWPATH = ''
for d in runningfile.split('/'):
    CMSSWPATH += d
    CMSSWPATH += '/'
    if 'CMSSW_' in d: break


print('Running analyzer from: ' + WORKPATH)

if __name__=='__main__':

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-t', '--tasktag', action='store', type=str, dest='tasktag', default='jobtag', help='tasktag')
    parser.add_option('-v', '--verbose', action='store_true', dest='verbose', help='verbose')
    parser.add_option('--submit', action='store_true', dest='submit', help='Submit the job if true')
    (opts, args) = parser.parse_args()

    # List of datasets 

    Datasets = {}
    Datasets['Run2016B_ver2_HIPM'] = '/gpfs/projects/cms/fernance/DoubleMuon/Run2016B_ver2-RunIISummer20UL16MiniAODAPVext/210216_134625'
    #Datasets['Run2016C_HIPM']      = '/gpfs/projects/cms/fernance/DoubleMuon/Run2016C-RunIISummer20UL16MiniAODAPVext/210224_085105'
    #Datasets['Run2016D_HIPM']      = '/gpfs/projects/cms/fernance/DoubleMuon/Run2016D-RunIISummer20UL16MiniAODAPVext/210217_220429'
    Datasets['Run2016E_HIPM']      = '/gpfs/projects/cms/fernance/DoubleMuon/Run2016E-RunIISummer20UL16MiniAODAPVext/210222_110506'
    #Datasets['Run2016F_HIPM']      = '/gpfs/projects/cms/fernance/DoubleMuon/Run2016F-RunIISummer20UL16MiniAODAPVext/210305_092358'
    Datasets['Run2016F_noHIPM']    = '/gpfs/projects/cms/fernance/DoubleMuon/Run2016F-RunIISummer20UL16MiniAODAPVext_noHIPM/210308_121101'
    Datasets['Run2016G_noHIPM']    = '/gpfs/projects/cms/fernance/DoubleMuon/Run2016G-RunIISummer20UL16MiniAODAPVext_noHIPM/210312_090055'
    #Datasets['Run2016H_noHIPM_1']  = '/gpfs/projects/cms/fernance/DoubleMuon/Run2016H-RunIISummer20UL16MiniAODAPVext_noHIPM/210318_145853'
    #Datasets['Run2016H_noHIPM_2']  = '/gpfs/projects/cms/fernance/DoubleMuon/Run2016H-RunIISummer20UL16MiniAODAPVext_noHIPM/210324_171842'


    for dataset in Datasets.keys():

        # Get path:
        inputpath = Datasets[dataset]

        # Create environment:
        if not os.path.exists(WORKPATH + 'test/2016/' + dataset): os.makedirs(WORKPATH + 'test/2016/' + dataset)
        logs_dir = WORKPATH + 'test/2016/' + dataset + '/logs/'
        confs_dir = WORKPATH + 'test/2016/' + dataset + '/confs/'
        scripts_dir = WORKPATH + 'test/2016/' + dataset + '/scripts/'
        if not os.path.exists(logs_dir): os.makedirs(logs_dir)
        if not os.path.exists(confs_dir): os.makedirs(confs_dir)
        if not os.path.exists(scripts_dir): os.makedirs(scripts_dir)

        # Define command:
        command = 'python ' + WORKPATH + 'scripts/scriptForSlurm.py'

        # cfg template
        if 'noHIPM' in dataset:
            command += ' -c ' + WORKPATH + 'test/2016/runLongLived2016postVFP_DATA_Muon_cfg.py'
        else:
            command += ' -c ' + WORKPATH + 'test/2016/runLongLived2016preVFP_DATA_Muon_cfg.py'

        # input
        command += ' -i ' + inputpath

        # release path
        command += ' -p ' + CMSSWPATH

        # work path
        command += ' -w ' + WORKPATH

        # splitting (to be calibrated)
        command += ' -n 30'

        # Output
        command += ' -o /gpfs/projects/cms/fernance/DoubleMuon/' + dataset + '_SLURM'

        # logs, cfs and scripts dir
        command += ' -l ' + logs_dir
        command += ' -f ' + confs_dir
        command += ' -s ' + scripts_dir


        #print(command)
        os.system(command)


