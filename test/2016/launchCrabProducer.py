import os

########################################
##   Declaration of global variables  ##
########################################

global dirs
global listOfBackgrs
global listOfSignals
global listOfData

##########################################
##  Initialization of global variables  ##
##########################################

dirs = os.listdir('./')


listOfBackgrounds = ['DYJetsToLL_M-50',
                     'DYJetsToLL_M-10to50',
                     'WJetsToLNu',
#                     'TT',
#                    'QCD_Pt-30to40',
#                    'QCD_Pt-40toInf',
                     'WW',
                     'WZ',
                     'ZZ']




listOfMuonData = ['DoubleMuon-Run2016B-07Aug17_ver2-v1',
                  'DoubleMuon-Run2016C-07Aug17-v1',
                  'DoubleMuon-Run2016D-07Aug17-v1',
                  'DoubleMuon-Run2016E-07Aug17-v1',
                  'DoubleMuon-Run2016F-07Aug17-v1',
                  'DoubleMuon-Run2016G-07Aug17-v1',
                  'DoubleMuon-Run2016H-07Aug17-v1'
                  ]

listOfEGData = ['DoubleEG-Run2016B-17Jul2018_ver2-v1',
                  'DoubleEG-Run2016C-17Jul2018-v1',
                  'DoubleEG-Run2016D-17Jul2018-v1',
                  'DoubleEG-Run2016E-17Jul2018-v1',
                  'DoubleEG-Run2016F-17Jul2018-v1',
                  'DoubleEG-Run2016G-17Jul2018-v1',
                  'DoubleEG-Run2016H-17Jul2018-v1'
                  ]

listOfSignals = [
#                 '1000-148-60_DisplacedSUSY',
#                 '120-48-165_DisplacedSUSY',
#                 '1500-494-160_DisplacedSUSY',
#                 '350-148-173_DisplacedSUSY',
                 '1000-150-100_ScalarBosons',
                 '1000-150-10_ScalarBosons',
                 '1000-350-350_ScalarBosons',
                 '1000-350-35_ScalarBosons',
                 '400-150-400_ScalarBosons',
                 '400-50-400_ScalarBosons',
                 '400-50-40_ScalarBosons',
                 '400-50-4_ScalarBosons'
                 ]


###########################
##  Function definition  ##
###########################

def hasBeenLaunched(sample):

    for d in dirs:
        if 'crab_'+sample in d: 
            return True

    return False



if __name__=="__main__":

    print(dirs)


    for sample in (listOfBackgrounds + listOfSignals + listOfEGData + listOfMuonData):

        if hasBeenLaunched(sample): continue
        os.system('crab submit ' + sample + '_crabSubmission.py')




