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

"""
listOfBackgrs = ['DYJetsToLL_M-50',
                 'DYJetsToLL_M-10to50',
                 'WJetsToLNu',
                 'TTJets_DiLept',
                 'QCD_Pt-30to40',
                 'QCD_Pt-40toInf',
                 'WW',
                 'WZ',
                 'ZZ']

"""


listOfBackgrounds = []

listOfSignals = [
#                 '1000-148-60_DisplacedSUSY',
#                 '120-48-165_DisplacedSUSY',
                 '1500-494-160_DisplacedSUSY',
#                 '350-148-173_DisplacedSUSY',
                 '1000-150-100_ScalarBosons',
#                 '1000-150-10_ScalarBosons',
                 '1000-350-350_ScalarBosons',
#                 '1000-350-35_ScalarBosons',
                 '400-150-400_ScalarBosons',
                 '400-150-40_ScalarBosons'
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


    for sample in (listOfBackgrounds + listOfSignals):

        if hasBeenLaunched(sample): continue
        os.system('crab submit ' + sample + '_crabSubmision.py')




