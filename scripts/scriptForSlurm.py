from optparse import OptionParser
import os
import stat
import glob
import sys
import os


###########################################################################################
###########################################################################################
### WARNING: For things to go smother and not overload the working dir this script should 
### be launched following the instructions in the README.md

###########################################################################################
theExecutable = """#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd CMSSWRELEASE
eval `scramv1 runtime -sh`
cd WORKAREA
cmsRun CONF
"""
###########################################################################################


if __name__ == "__main__":


    parser = OptionParser(usage="%prog --help")

    parser.add_option("-c", "--conf",           dest="conf",             type="string",      default='conf.py',               help="CMSSW configuration file.")
    parser.add_option("-i", "--input",          dest="inputDir",         type="string",      default='./',                    help="Name of the input directory.")
    parser.add_option("-p", "--path",           dest="path",             type="string",      default='./CMSSW_X_Y_Z/src/',    help="Path of the CMSSW release.")
    parser.add_option("-w", "--workingpath",    dest="workingpath",      type="string",      default='./',                    help="Path of the working area.")
    parser.add_option("-n", "--nfilejob",       dest="nfilejob",         type="int",         default=5,                       help="Number of files per job.")
    parser.add_option("-N", "--neventsjob",     dest="neventsjob",       type="int",         default=-1,                      help="Number of events per job.")
    parser.add_option("-o", "--output",         dest="output",           type="string",      default='.',                     help="Output directory.")
    parser.add_option("-l", "--logs",           dest="logs",             type="string",      default='.',                     help="Logs directory.")
    parser.add_option("-f", "--cfgs",           dest="cfgs",             type="string",      default='.',                     help="Cfgs directory.")
    parser.add_option("-s", "--scripts",        dest="scripts",          type="string",      default='.',                     help="Scripts directory.")

    (options, args) = parser.parse_args()

    #We take the list of files and put in a stream
    listOfFiles = glob.glob(options.inputDir + '/*/*.root')
    if len(listOfFiles) < 1:
        print 'No root files were found in the given input directory'
        sys.exit()
 
    #Here we do the arithmetics of the jobs
    nFiles = len(listOfFiles)
    nJobs = int(float(nFiles)/float(options.nfilejob))
    listOfFilesPerJob = []
    for i in range(0, nJobs+1):
        thefilesforthisjob = ''
        for j in range(0, options.nfilejob):
            index = i * options.nfilejob + j
            if index == nFiles:
                break
            thefilesforthisjob += "'file:" + (listOfFiles[index]) + "',"
        thefilesforthisjob = thefilesforthisjob[:-1]
        listOfFilesPerJob.append(thefilesforthisjob)
                          
    #Create the output dir if it does not exist yet
    if not os.path.exists(options.output):
        os.makedirs(options.output)


    #Read the configuration file: make sure that the configuration file has the labels INPUT, OUTPUT and MAXEVENTS to be replaced by the code 
    confText = open(options.conf).read()
    confText = confText.replace("MAXEVENTS", str(options.neventsjob))
    
    #Read the bash file to be executed by condor
    theExecutable = theExecutable.replace('CMSSWRELEASE', options.path)
    theExecutable = theExecutable.replace('WORKAREA', options.workingpath)
    
    bash_path = options.scripts if options.scripts[-1] == '/' else options.scripts + '/'
    thefile = open(bash_path + 'run.sh', 'w')
    #Running on jobs
    for i, inputfiles in enumerate(listOfFilesPerJob):

        conf_path = options.cfgs if options.cfgs[-1] == '/' else options.cfgs + '/'
        confFileName = conf_path + 'conf_' + str(i) + ".py"
        outputFileName = "'" + options.output + "/" + "output_" + str(i) + ".root'"
        bashFileName = bash_path + 'script_' + str(i) + ".sh"
        myConf = confText.replace('INPUT', inputfiles)
        myConf = myConf.replace('OUTPUT', outputFileName)
        myBash = theExecutable.replace('CONF', confFileName)
        fConf = open(confFileName, 'w')
        fConf.write(myConf)
        fConf.close()
        fBash = open(bashFileName, 'w')
        fBash.write(myBash)
        fBash.close()
        st = os.stat(bashFileName)
        os.chmod(bashFileName, st.st_mode | stat.S_IEXEC)
        logFileName = options.logs + '/log_' + str(i) + '.out'
        errFileName = options.logs + '/log_' + str(i) + '.err'
        addendum = 'sbatch ' + ' -o ' + logFileName + ' -e ' + errFileName + ' --qos=gridui_medium --partition=cloudcms ' + bashFileName + '\n'
        thefile.write(addendum)

    thefile.close()

  


