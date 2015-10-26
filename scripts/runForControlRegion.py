#Quick script to loop over all channels to make fitting plots

#from ROOT import *

import subprocess
import sys

channels = {"1":"eee","2":"eemu","4":"emumu","8":"mumumu","16":"eeeInv","32":"eemuInv","64":"emumuInv","128":"mumumuInv"}

channelList = {"eee":"17","eemu":"34","emumu":"68","mumu":"136"}

metCut = sys.argv[1]
metStr = metCut.split(".")[0]

mtwCut = sys.argv[2]
mtwStr = mtwCut.split(".")[0]

dirPostfix = "ctrl"
if len(sys.argv) > 3:
    dirPostfix = sys.argv[3]

#Make the skim directory
subprocess.call("mkdir mvaDirs/skims/met"+metStr+"mtw"+mtwStr+dirPostfix,shell=True)

for chanName in channelList.keys():
#    print "bin/analysisMain.exe -c configs/"+chanName+"Conf.txt -u -t -k "+str(channelList[chanName])+" --jetRegion 1,1,2,4 -v 4095 -z --mvaDir mvaDirs/skims/met"+metStr+"mtw"+mtwStr + "/ --metCut " + str(metCut) + " --mtwCut " + str(mtwCut)
    subprocess.call("bin/analysisMain.exe -c configs/"+chanName+"Conf.txt -u -t -k "+str(channelList[chanName])+" --jetRegion 1,1,3,4 -v 4095 -z --mvaDir mvaDirs/skims/met"+metStr+"mtw"+mtwStr+dirPostfix + "/ --metCut " + str(metCut) + " --mtwCut " + str(mtwCut),shell=True)

#subprocess.call("bin/analysisMain.exe -c configs/wzSystsConf.txt -u -t -k 15 --jetRegion 1,1,3,4 -z --mvaDir mvaDirs/skims/met"+metStr+"mtw"+mtwStr+dirPostfix + "/ --metCut " + str(metCut) + " --mtwCut " + str(mtwCut),shell=True)

#Make the mvaInput directory
subprocess.call("mkdir mvaDirs/inputs/met"+metStr+"mtw"+mtwStr+dirPostfix,shell=True)

print "python scripts/makeMVAInput.py [\\\"eee\\\",\\\"eemu\\\",\\\"emumu\\\",\\\"mumumu\\\"] mvaDirs/skims/met"+metStr+"mtw"+mtwStr+dirPostfix+"/ mvaDirs/inputs/met"+metStr+"mtw"+mtwStr+dirPostfix+"/ -s"
subprocess.call("python scripts/makeMVAInput.py [\\\"eee\\\",\\\"eemu\\\",\\\"emumu\\\",\\\"mumumu\\\"] mvaDirs/skims/met"+metStr+"mtw"+mtwStr+dirPostfix+"/ mvaDirs/inputs/met"+metStr+"mtw"+mtwStr+dirPostfix+"/ -s",shell=True)

