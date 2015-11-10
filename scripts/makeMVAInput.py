#A tool to pull a load of information from the mva trees. Should be all easy jazz...

from ROOT import *

import numpy as n
import sys
import os
from array import array
from jetCorrectionUncertainty import JetCorrectionUncertainty

def sortOutLeptons(tree,channel):
    ###Returns three TLorentzVectors containing the two z leptons and the w lepton. This will be VERY useful for making all of the plots.
    #Reads the position of the w and z leptons from variables stored at mvaTree making time, because I'm great and finally got around to doing it.
    zMass = 100
    zLep1,zLep2,wLep = 0,0,0
    #Let's try commenting this out and see if everything breaks? Hopefully it won't do...
    #if tree.numElePF2PAT < 3:
    #    return (0,0,0)
    if channel == "eee":
        wLep = TLorentzVector(tree.elePF2PATGsfPx[tree.wLepIndex],tree.elePF2PATGsfPy[tree.wLepIndex],tree.elePF2PATGsfPz[tree.wLepIndex],tree.elePF2PATGsfE[tree.wLepIndex])
        zLep1 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep1Index],tree.elePF2PATGsfPy[tree.zLep1Index],tree.elePF2PATGsfPz[tree.zLep1Index],tree.elePF2PATGsfE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep2Index],tree.elePF2PATGsfPy[tree.zLep2Index],tree.elePF2PATGsfPz[tree.zLep2Index],tree.elePF2PATGsfE[tree.zLep2Index])
    if channel == "eemu":
        zLep1 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep1Index],tree.elePF2PATGsfPy[tree.zLep1Index],tree.elePF2PATGsfPz[tree.zLep1Index],tree.elePF2PATGsfE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.elePF2PATGsfPx[tree.zLep2Index],tree.elePF2PATGsfPy[tree.zLep2Index],tree.elePF2PATGsfPz[tree.zLep2Index],tree.elePF2PATGsfE[tree.zLep2Index])
        wLep = TLorentzVector(tree.muonPF2PATPx[tree.wLepIndex],tree.muonPF2PATPy[tree.wLepIndex],tree.muonPF2PATPz[tree.wLepIndex],tree.muonPF2PATE[tree.wLepIndex])
    if channel == "emumu":
        zLep1 = TLorentzVector(tree.muonPF2PATPx[tree.zLep1Index],tree.muonPF2PATPy[tree.zLep1Index],tree.muonPF2PATPz[tree.zLep1Index],tree.muonPF2PATE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.muonPF2PATPx[tree.zLep2Index],tree.muonPF2PATPy[tree.zLep2Index],tree.muonPF2PATPz[tree.zLep2Index],tree.muonPF2PATE[tree.zLep2Index])
        wLep = TLorentzVector(tree.elePF2PATGsfPx[tree.wLepIndex],tree.elePF2PATGsfPy[tree.wLepIndex],tree.elePF2PATGsfPz[tree.wLepIndex],tree.elePF2PATGsfE[tree.wLepIndex])
    if channel == "mumumu":
        zLep1 = TLorentzVector(tree.muonPF2PATPx[tree.zLep1Index],tree.muonPF2PATPy[tree.zLep1Index],tree.muonPF2PATPz[tree.zLep1Index],tree.muonPF2PATE[tree.zLep1Index])
        zLep2 = TLorentzVector(tree.muonPF2PATPx[tree.zLep2Index],tree.muonPF2PATPy[tree.zLep2Index],tree.muonPF2PATPz[tree.zLep2Index],tree.muonPF2PATE[tree.zLep2Index])
        wLep = TLorentzVector(tree.muonPF2PATPx[tree.wLepIndex],tree.muonPF2PATPy[tree.wLepIndex],tree.muonPF2PATPz[tree.wLepIndex],tree.muonPF2PATE[tree.wLepIndex])
    return (zLep1,zLep2,wLep)

def getJets(tree,syst,jetUnc,met):
    #Makes a short list of indices of the jets in the event
    jetList = []
    jetVecList = []
    for i in range(15):
        if tree.jetInd[i] > len(tree.jetPF2PATEta): continue
        if tree.jetInd[i] > -.5:
            jetList.append(tree.jetInd[i])
            jetVecList.append(getJetVec(tree,tree.jetInd[i],syst,jetUnc,met))
        else: continue
    return (jetList,jetVecList)

def getBjets(tree,syst,jetUnc,met,jets):
    #Return a list of the indices of the b-jets in the event
    bJetList = []
    bJetVecList = []
    for i in range(10):
        if tree.bJetInd[i] > len(tree.jetPF2PATEta): continue
        if tree.bJetInd[i] > -0.5:
            bJetList.append(tree.bJetInd[i])
            bJetVecList.append(getJetVec(tree,jets[tree.bJetInd[i]],syst,jetUnc,met))
        else:continue
#    print len(bJetList)
    return (bJetList,bJetVecList)

def getJetVec(tree, index, syst, jetUnc, metVec):
    #Gets a vector for a jet and applies jet corrections.
    oldSmearCorr = 0.;
    if (fabs(tree.jetPF2PATEta[index]) <= 1.1):
        oldSmearCorr = 0.066
    elif (fabs(tree.jetPF2PATEta[index]) <= 1.7):
        oldSmearCorr = 0.191
    elif (fabs(tree.jetPF2PATEta[index]) <= 2.3):
        oldSmearCorr = 0.096;
    else:
        oldSmearCorr = 0.166;
    oldSmearValue = 1.0;
    if (jetUnc and tree.genJetPF2PATPT[index] > -990.) :
        oldSmearValue = max(0.0, tree.jetPF2PATPtRaw[index] + (tree.jetPF2PATPtRaw[index] - tree.genJetPF2PATPT[index]) * oldSmearCorr)/tree.jetPF2PATPtRaw[index];
    newJECCorr = 0.;
    if (fabs(tree.jetPF2PATEta[index]) <= 0.5) :
        newJECCorr = 1.079;
        if (syst == 16): newJECCorr = 1.105;
        elif (syst == 32): newJECCorr = 1.053;
    elif (fabs(tree.jetPF2PATEta[index]) <= 1.1) :
        newJECCorr = 1.099;
        if (syst == 16): newJECCorr = 1.127;
        elif (syst == 32): newJECCorr = 1.071;
    elif (fabs(tree.jetPF2PATEta[index]) <= 1.7) :
        newJECCorr = 1.121;
        if (syst == 16): newJECCorr = 1.150;
        elif (syst == 32): newJECCorr = 1.092
    elif (fabs(tree.jetPF2PATEta[index]) <= 2.3) :
        newJECCorr = 1.208;
        if (syst == 16): newJECCorr = 1.254;
        elif (syst == 32): newJECCorr = 1.162;
    elif (fabs(tree.jetPF2PATEta[index]) <= 2.8):
        newJECCorr = 1.254;
        if (syst == 16): newJECCorr = 1.316;
        elif (syst == 32): newJECCorr = 1.192;
    elif (fabs(tree.jetPF2PATEta[index]) <= 3.2):
        newJECCorr = 1.395;
        if (syst == 16): newJECCorr = 1.458;
        elif (syst == 32): newJECCorr = 1.332;
    else :
        newJECCorr = 1.056;
        if (syst == 16): newJECCorr = 1.247;
        elif (syst == 32): newJECCorr = 0.865;
    newSmearValue = 1.0;
    if (jetUnc and tree.genJetPF2PATPT[index] > -990.): newSmearValue = max(0.0, tree.jetPF2PATPtRaw[index] + (tree.jetPF2PATPtRaw[index] - tree.genJetPF2PATPT[index]) * newJECCorr)/tree.jetPF2PATPtRaw[index];
    returnJet = TLorentzVector(0);
    if (newSmearValue < 0.01): returnJet.SetPxPyPzE(0.0001,0.0001,0.0001,0.0001);
    else :
        #Propogate through the met. But only do it if the smear jet isn't 0.
        metVec.SetPx(metVec.Px()+tree.jetPF2PATPx[index])
        metVec.SetPy(metVec.Py()+tree.jetPF2PATPy[index])
        returnJet.SetPxPyPzE(newSmearValue*tree.jetPF2PATPx[index]/oldSmearValue,newSmearValue*tree.jetPF2PATPy[index]/oldSmearValue,newSmearValue*tree.jetPF2PATPz[index]/oldSmearValue,newSmearValue*tree.jetPF2PATE[index]/oldSmearValue);
    if syst == 16:
        returnJet *= 1+ jetUnc.getUncertainty(returnJet.Pt(), returnJet.Eta(),1)
    elif syst == 32:
        returnJet *= 1+ jetUnc.getUncertainty(returnJet.Pt(), returnJet.Eta(),2)
    metVec.SetPx(metVec.Px()-returnJet.Px())
    metVec.SetPy(metVec.Py()-returnJet.Py())
    return returnJet

def doUncMet(tree,met,zLep1,zLep2,wLep,jetVecs,syst):
    #Subtracts all items from met and then varies what's left by 10% for systematic purposes.
    uncMetX = met.Px() + zLep1.Px() + zLep2.Px() + wLep.Px()
    uncMetY = met.Py() + zLep1.Py() + zLep2.Py() + wLep.Py()
    for i in range(len(jetVecs)):
        uncMetX += jetVecs[i].Px()
        uncMetY += jetVecs[i].Py()
    if syst == 1028:
        met.SetPx(met.Px() + 0.1*uncMetX)
        met.SetPy(met.Py() + 0.1*uncMetY)
    else:
        met.SetPx(met.Px() - 0.1*uncMetX)
        met.SetPy(met.Py() - 0.1*uncMetY)
    return met

def setupInputVars():
    #Make the variables we want to save
    inputVars = {}
    inputVars["eventWeight"] = array('f',[0.])
    inputVars["mTW"] = array('f',[0.])
    inputVars["leptWPt"] = array('f',[0.])
    inputVars["leadJetPt"] = array('f',[0.])
    inputVars["totPt"] = array('f',[0.])
    inputVars["totEta"] = array('f',[0.])
    inputVars["totPtVec"] = array('f',[0.])
    inputVars["totVecM"] = array('f',[0.])
    inputVars["chan"] = array('i',[0])
    inputVars["nJets"] = array('f',[0.])
    inputVars["nBjets"] = array('f',[0.])
    inputVars["met"] = array('f',[0.])
    inputVars["lepPt"] = array('f',[0.])
    inputVars["lepMetPt"] = array('f',[0.])
    inputVars["totPt2Jet"] = array('f',[0.])
    inputVars["leadJetbTag"] = array('f',[0.])
    inputVars["leadJetEta"] = array('f',[0.])
    inputVars["secJetbTag"] = array('f',[0.])
    inputVars["secJetPt"] = array('f',[0.])
    inputVars["secJetEta"] = array('f',[0.])
    inputVars["bTagDisc"] = array('f',[0.])
    inputVars["topMass"] = array('f',[0.])
    inputVars["topPt"] = array('f',[0.])
    inputVars["topEta"] = array('f',[0.])
    inputVars["zPt"] = array('f',[0.])
    inputVars["zEta"] = array('f',[0.])
    inputVars["wLepEta"] = array('f',[0.])
    inputVars["wZdelR"] = array('f',[0.])
    inputVars["j1j2delR"] = array('f',[0.])
    inputVars["minZJetR"] = array('f',[0.])
    inputVars["wLepjminR"] = array('f',[0.])
    inputVars["zWLepdelR"] = array('f',[0.])
    inputVars["zmetdelPhi"] = array('f',[0.])
    inputVars["zWLepdelPhi"] = array('f',[0.])
    inputVars["lbDelR"] = array('f',[0.])
    inputVars["lbDelPhi"] = array('f',[0.])
    inputVars["zlb1DelR"] = array('f',[0.])
    inputVars["zlb1DelPhi"] = array('f',[0.])
    inputVars["zlb2DelR"] = array('f',[0.])
    inputVars["zlb2DelPhi"] = array('f',[0.])
    inputVars["totHt"] = array('f',[0.])
    inputVars["lepHt"] = array('f',[0.])
    inputVars["jetHt"] = array('f',[0.])
    inputVars["lepMetHt"] = array('f',[0.])
    inputVars["totHtOverPt"] = array('f',[0.])
    inputVars["zMass"] = array('f',[0.])
    return inputVars

def setupBranches(tree,varMap):
    tree.Branch("EvtWeight", varMap["eventWeight"], "EvtWeight/F")
    tree.Branch("mTW", varMap["mTW"], "mTW/F")
    tree.Branch("leptWPt", varMap["leptWPt"], "leptWPt/F")
    tree.Branch("leadJetPt",varMap["leadJetPt"],"leadJetPt/F")
    tree.Branch("totPt",varMap["totPt"],"totPt/F")
    tree.Branch("totEta",varMap["totEta"],"totEta/F")
    tree.Branch("totPtVec",varMap["totPtVec"],"totPtVec/F")
    tree.Branch("totVecM",varMap["totVecM"],"totVecM/F")
    tree.Branch("Channel",varMap["chan"],"Channel/I")
    tree.Branch("NJets",varMap["nJets"],"NJets/F")
    tree.Branch("NBJets",varMap["nBjets"],"NBJets/F")
    tree.Branch("met",varMap["met"],"met/F")
    tree.Branch("lepPt",varMap["lepPt"],"lepPt/F")
    tree.Branch("lepMetPt",varMap["lepMetPt"],"lepMetPt/F")
    tree.Branch("totPt2Jet",varMap["totPt2Jet"],"totPt2Jet/F")
    tree.Branch("btagDiscri",varMap["bTagDisc"],"btagDiscri/F")
    tree.Branch("leadJetbTag",varMap["leadJetbTag"],"leadJetbTag/F")
    tree.Branch("leadJetEta",varMap["leadJetEta"],"leadJetEta/F")
    tree.Branch("secJetbTag",varMap["secJetbTag"],"secJetbTag/F")
    tree.Branch("secJetPt",varMap["secJetPt"],"secJetPt/F")
    tree.Branch("secJetEta",varMap["secJetEta"],"secJetEta/F")
    tree.Branch("topMass",varMap["topMass"],"topMass/F")
    tree.Branch("topPt",varMap["topPt"],"topPt/F")
    tree.Branch("topEta",varMap["topEta"],"topEta/F")
    tree.Branch("Zpt",varMap["zPt"],"Zpt/F")
    tree.Branch("Zeta",varMap["zEta"],"Zeta/F")
    tree.Branch("leptWEta",varMap["wLepEta"],"leptWEta/F")
    tree.Branch("wzdelR",varMap["wZdelR"],"wzdelR/F")
    tree.Branch("jjdelR",varMap["j1j2delR"],"jjdelR/F")
    tree.Branch("zjminR",varMap["minZJetR"],"zjminR/F")
    tree.Branch("wLepjminR",varMap["wLepjminR"],"wLepjminR/F")
    tree.Branch("ZlepWdelPhi",varMap["zWLepdelPhi"],"ZlepWdelPhi/F")
    tree.Branch("ZmetdelPhi",varMap["zmetdelPhi"],"ZmetdelPhi/F")
    tree.Branch("ZlepWdelR",varMap["zWLepdelR"],"ZlepWdelR/F")
    tree.Branch("lbDelR",varMap["lbDelR"],"lbDelR/F")
    tree.Branch("lbDelPhi",varMap["lbDelPhi"],"lbDelPhi/F")
    tree.Branch("zlb1DelR",varMap["zlb1DelR"],"zlb1DelR/F")
    tree.Branch("zlb1DelPhi",varMap["zlb1DelPhi"],"zlb1DelPhi/F")
    tree.Branch("zlb2DelR",varMap["zlb2DelR"],"zlb2DelR/F")
    tree.Branch("zlb2DelPhi",varMap["zlb2DelPhi"],"zlb2DelPhi/F")
    tree.Branch("totHt",varMap["totHt"],"totHt/F")
    tree.Branch("lepHt",varMap["lepHt"],"lepHt/F")
    tree.Branch("jetHt",varMap["jetHt"],"jetHt/F")
    tree.Branch("lepMetHt",varMap["lepMetHt"],"lepMetHt/F")
    tree.Branch("totHtOverPt",varMap["totHtOverPt"],"totHtOverPt/F")
    tree.Branch("zMass",varMap["zMass"],"zMass/F")



def fillTree(outTree, outTreeSdBnd, varMap, tree, label, channel, jetUnc, overRideWeight = -1., zPtEventWeight = 0.):
    #Fills the output tree. This is a new function because I want to access data and MC in different ways but do the same thing to them in the end.

    syst = 0
    if "__jer__plus" in label:
        syst = 4
    if "__jer__minus" in label:
        syst = 8
    if "__jes__plus" in label:
        syst = 16
    if "__jes__minus" in label:
        syst = 32
    if "__met__plus" in label:
        syst = 1024
    if "__met__minus" in label:
        syst = 2048
    if channel == "eee":
        varMap["chan"][0] = 3
    if channel == "eemu":
        varMap["chan"][0] = 2
    if channel == "emumu":
        varMap["chan"][0] = 1
    if channel == "mumumu":
        varMap["chan"][0] = 0
    #        topMass
    #Set up those variables as branches
    for event in range(tree.GetEntries()):
            #Fill some plots here. Let's make an example mTW plot.
            #Make a config that'll do this for me? I've done these before so should be easy. Fill expressions could be a pain?
        tree.GetEntry(event)
#        if tree.eventWeight < 0.000000001: continue
        (zLep1,zLep2,wLep) = sortOutLeptons(tree,channel)
        metVec = TLorentzVector(tree.metPF2PATPx,tree.metPF2PATPy,0,tree.metPF2PATEt)
        (jets,jetVecs) = getJets(tree,syst,jetUnc,metVec)
        (bJets,bJetVecs) = getBjets(tree,syst,jetUnc,metVec,jets)
        #Do unclustered met stuff here now that we have all of the objects, all corrected for their various SFs etc.
        if syst == 1024 or syst == 2048:
            metVec = doUncMet(tree,metVec,zLep1,zLep2,wLep,jetVecs,syst)
        if wLep == 0:
            continue
        if tree.eventWeight == tree.eventWeight:
            varMap["eventWeight"][0] = tree.eventWeight
            if overRideWeight > 0.:
                varMap["eventWeight"][0] = tree.eventWeight * overRideWeight
            if zPtEventWeight > 0.1:
                varMap["eventWeight"][0] *= tree.eventWeight
            if zPtEventWeight < -0.1:
                varMap["eventWeight"][0] = 1.
        else:
            varMap["eventWeight"][0] = 0.

        if varMap["eventWeight"][0] < 0.:
            varMap["eventWeight"][0] = 0.
        if varMap["eventWeight"][0] > 100:
            varMap["eventWeight"][0] = 0.
        varMap["leptWPt"][0] = wLep.Pt()
        varMap["wLepEta"][0] = wLep.Eta()
        varMap["leadJetPt"][0] = jetVecs[0].Pt()
        varMap["leadJetEta"][0] = jetVecs[0].Eta()
        #Make all the random Pt variables I'm saving for some reason
        totPx,totPy = 0,0
        totPx += zLep1.Px() + zLep2.Px() + wLep.Px()
        totPy += zLep1.Py() + zLep2.Py() + wLep.Py()
        varMap["lepPt"][0] = sqrt(totPx * totPx + totPy * totPy)
        totPx += metVec.Px()
        totPy += metVec.Py()
        varMap["lepMetPt"][0] = sqrt(totPx * totPx + totPy * totPy)
        totPx += jetVecs[0].Px()
        totPy += jetVecs[0].Py()
        if len(jetVecs) > 1:
            totPx += jetVecs[1].Px()
            totPy += jetVecs[1].Py()
        varMap["totPt2Jet"][0] = sqrt(totPx * totPx + totPy * totPy)
        for i in range(2,len(jets)):
            totPx+=jetVecs[i].Px()
            totPy+=jetVecs[i].Py()
        varMap["totPt"][0] = sqrt(totPx * totPx + totPy * totPy)
        totVec = (wLep + zLep1+zLep2)
        for i in range(len(jetVecs)):
            totVec += jetVecs[i]
        varMap["totEta"][0] = totVec.Eta()
        varMap["totPtVec"][0] = totVec.Pt()
        varMap["totVecM"][0] = totVec.M()
        varMap["mTW"][0] = sqrt(2*metVec.Pt()*wLep.Pt() * (1-cos(metVec.Phi() - wLep.Phi())))
        varMap["nJets"][0] = float(len(jets))
        varMap["nBjets"][0] = float(len(bJets))
        varMap["met"][0] = metVec.Pt()
        varMap["bTagDisc"][0] = 0.
        if len(bJets) > 0: varMap["bTagDisc"][0] = tree.jetPF2PATBDiscriminator[jets[bJets[0]]]
        varMap["leadJetbTag"][0] = tree.jetPF2PATBDiscriminator[jets[0]]
        varMap["secJetbTag"][0] = -10.
        varMap["secJetPt"][0] = -1.
        varMap["secJetEta"][0] = -500.
        if len(jetVecs) > 1:
            varMap["secJetbTag"][0] = tree.jetPF2PATBDiscriminator[jets[1]]
            varMap["secJetPt"][0] = jetVecs[1].Pt()
            varMap["secJetEta"][0] = jetVecs[1].Eta()

#        print bTagDisc[0], bJets[0], tree.jetPF2PATBDiscriminator[jets[bJets[0]]], len(bJets), nBjets[0]
        varMap["topMass"][0] = 0.
        varMap["topPt"][0] = 0.
        varMap["topEta"][0] = 0.
        if len(bJetVecs) > 0:
            varMap["topMass"][0] = (bJetVecs[0] + metVec + wLep).M()
            varMap["topPt"][0] = (bJetVecs[0] + metVec + wLep).Pt()
            varMap["topEta"][0] = (bJetVecs[0] + metVec + wLep).Eta()
        varMap["zPt"][0] = (zLep2 + zLep1).Pt()
        varMap["zEta"][0] = (zLep2 + zLep1).Eta()
        varMap["wZdelR"][0] = (zLep2 + zLep1).DeltaR(metVec + wLep)
        varMap["j1j2delR"][0] = -1.
        if len(jetVecs) > 1:
            varMap["j1j2delR"][0] = jetVecs[0].DeltaR(jetVecs[1])
        varMap["minZJetR"][0] = 3.0
        varMap["wLepjminR"][0] = 3.0
        jetHt = 0.
        for i in range(len(jetVecs)):
            jetHt += jetVecs[i].Pt()
            if jetVecs[i].DeltaR(wLep) < varMap["wLepjminR"][0]:
                varMap["wLepjminR"][0] = jetVecs[i].DeltaR(wLep)
            if jetVecs[i].DeltaR(zLep2 + zLep1) < varMap["minZJetR"][0]:
                varMap["minZJetR"][0] = jetVecs[i].DeltaR(zLep2 + zLep1)
        varMap["zWLepdelR"][0] = (zLep2 + zLep1).DeltaR(wLep)
        varMap["zmetdelPhi"][0] = (zLep2+zLep1).DeltaPhi(metVec)
        varMap["zWLepdelPhi"][0] = (zLep2 + zLep1).DeltaPhi(wLep)        
        varMap["lbDelR"][0] = 0.
        if len(bJets) > 0: varMap["lbDelR"][0] = wLep.DeltaR(jetVecs[bJets[0]])
        varMap["lbDelPhi"][0] = 0.
        if len(bJets) > 0: varMap["lbDelPhi"][0] = wLep.DeltaPhi(jetVecs[bJets[0]])
        varMap["zlb1DelR"][0] = 0.
        if len(bJets) > 0: varMap["zlb1DelR"][0] = zLep1.DeltaR(jetVecs[bJets[0]])
        varMap["zlb1DelPhi"][0] = 0.
        if len(bJets) > 0: varMap["zlb1DelPhi"][0] = zLep1.DeltaPhi(jetVecs[bJets[0]])
        varMap["zlb2DelR"][0] = 0.
        if len(bJets) > 0: varMap["zlb2DelR"][0] = zLep2.DeltaR(jetVecs[bJets[0]])
        varMap["zlb2DelPhi"][0] = 0.
        if len(bJets) > 0: varMap["zlb2DelPhi"][0] = zLep2.DeltaPhi(jetVecs[bJets[0]])
        ht = 0.
        ht += zLep1.Pt() + zLep2.Pt() + wLep.Pt()
        varMap["lepHt"][0] = ht
        ht += metVec.Pt()
        varMap["jetHt"][0] = jetHt
        varMap["lepMetHt"][0] = ht
        ht += jetHt
        varMap["totHt"][0] = ht
        varMap["totHtOverPt"][0] = ht /  sqrt(totPx * totPx + totPy * totPy) 
        varMap["zMass"][0] = (zLep1+zLep2).M()

        if varMap["nJets"][0] == 1 and outTreeSdBnd:
            outTreeSdBnd.Fill()
        else:
            outTree.Fill()


def main():

    zEnrichWeights = {"mumumu":-1.,"emumu":-1.,"eemu":-1.,"eee":-1.}
#    zEnrichWeights = {"mumumu":0.0001,"emumu":0.496,"eemu":0.0283,"eee":0.211}
    #lepton selection stage weights
#    zEnrichWeights = {"mumumu":0.057,"emumu":0.612,"eemu":0.0637,"eee":0.216}

    #Mapping of our mc names to IPHC names
    listOfMCs = {"WW2l2nu":"WW","WZ3l1nu":"WZ","ZZ4l":"ZZ","sChannel":"TsChan","sbarChannel":"TbarsChan","tChannel":"TtChan","tbarChannel":"TbartChan","tWInclusive":"TtW","tbarWInclusive":"TbartW","tZq":"tZq","tZq4Flavour3Lepton":"tZq4f","ttW":"TTW","ttZ":"TTZ","ttbarDilepton":"TT","wPlusJets":"Wjets", "zPlusJets10To50Filter":"DYToLL_M10-50","zPlusJetsTuneZ2Star":"Zjets"}
#    listOfMCs = {"WW2l2nu":"WW","WZ3l1nu":"WZ","sChannel":"TsChan","sbarChannel":"TbarsChan","tChannel":"TtChan","tbarChannel":"TbartChan","tWInclusive":"TtW","tbarWInclusive":"TbartW","ttZ":"TTZ","wPlusJets":"Wjets", "zPlusJets10To50Filter":"DYToLL_M10-50","zPlusJetsTuneZ2Star":"Zjets"}
    #Reduced sets of data for testing and running over restricted samples
#    listOfMCs = {"WW2l2nu":"WW","WZ3l1nu":"WZ","WZ2l2nu":"WZ","ZZ2l2q":"ZZ","ZZ4l":"ZZ","sChannel":"TsChan","sbarChannel":"TbarsChan","tChannel":"TtChan","tbarChannel":"TbartChan","tWInclusive":"TtW","tbarWInclusive":"TbartW","ttW":"TTW","ttZ":"TTZ","ttbarDilepton":"TT","wPlusJets":"Wjets", "zPlusJets10To50Filter":"DYToLL_M10-50","zPlusJetsTuneZ2Star":"Zjets"}
    #listOfMCs = {"ZZ4l":"ZZ"}
#    listOfMCs = {"ttW":"TTW","ZZ4l":"ZZ"}
#    listOfMCs = {"WZ3l1nu":"WZ"}
#    listOfMCs = {"WZ3l1nu":"WZ","zPlusJetsTuneZ2Star":"Zjets"}
#    listOfMCs = {"zPlusJets10To50Filter":"DYToLL_M10-50"}
 #   listOfMCs = {}

    #Set-up JEC corrections
    jetUnc = JetCorrectionUncertainty("scripts/Summer13_V5_MC_Uncertainty_AK5PFchs.txt")

    #mapping of channels to dataTypes
    channelToDataset = {"eee":"DataEG","eemu":"DataMuEG","emumu":"DataMuEG","mumumu":"DataMu"}

    #systematics list
    systs = ["","__trig__plus","__trig__minus","__jer__plus","__jer__minus","__jes__plus","__jes__minus","__pileup__plus","__pileup__minus","__met__plus","__met__minus","__bTag__plus","__bTag__minus","__pdf__plus","__pdf__minus","__zPt__plus","__zPt__minus"]
    #Reduced syst list for testing
#    systs = ["" ,"__bTag__plus"]
    #read what channel we're using here - changing this so that multiple things can be stored in the same file. i.e. should now be a list of channels to run over
    channels = eval(sys.argv[1])

    #Might make this customisable later, but don't really feel like I need to yet
    inputDir = "mvaTest/"
    if len(sys.argv) > 2:
        inputDir = sys.argv[2]
    
    outputDir = "mvaInputs/"
    if len(sys.argv) > 3:
        outputDir = sys.argv[3]
    inputVars = setupInputVars()

    useSidebandRegion = False
    if len(sys.argv) > 4 and sys.argv[4] == "-s":
        useSidebandRegion = True
    treeNamePostfixSig = ""
    treeNamePostfixSB = ""
    if useSidebandRegion:
        print "Using control region"
        treeNamePostfixSig = "sig_"
        treeNamePostfixSB = "ctrl_"

    #Loop over samples
    for sample in listOfMCs.keys():
        overrideWeight = -1.
        if sample == "WZ2l2nu" or sample == "ZZ2l2p":
            continue
        print "Doing " + sample + ": ",
        sys.stdout.flush()
        if "WZ" in sample:
            overrideWeight = 1286770./2133868.
        
        outFile = 0
        #update the appropriate root file
        outFile = TFile(outputDir+"histofile_"+listOfMCs[sample] + ".root","RECREATE")

        for syst in systs:
            #We now define the outtree out here, coz it seems like a more sensible option.
            outTreeSig = TTree("Ttree_"+treeNamePostfixSig+listOfMCs[sample]+syst, "Ttree_"+treeNamePostfixSig+listOfMCs[sample]+syst)
            outTreeSdBnd = 0
            if useSidebandRegion:
                outTreeSdBnd = TTree("Ttree_"+treeNamePostfixSB+listOfMCs[sample]+syst, "Ttree_"+treeNamePostfixSB+listOfMCs[sample]+syst)
                setupBranches(outTreeSdBnd,inputVars)
            setupBranches(outTreeSig,inputVars)
            for channel in channels:
                inFile = TFile(inputDir+sample+channel+"mvaOut.root","READ")
                if "met" in syst or "zPt" in syst:
                    tree = inFile.Get("tree")
                else:
                    tree = inFile.Get("tree"+syst)
                print syst + ": " + str(tree.GetEntriesFast()),
                sys.stdout.flush()
                #Various stuff needs to be saved in the same trees. Create new one if it doesn't exist, open current one if it does            
                fillTree(outTreeSig, outTreeSdBnd, inputVars, tree, listOfMCs[sample]+syst, channel, jetUnc, overRideWeight = overrideWeight)
                inFile.Close()
            outFile.cd()
            outFile.Write()
            outTreeSig.Write()
            if useSidebandRegion:
                outTreeSdBnd.Write()
        #if tree exists just update that.
        #        if outFile.GetListOfKeys().Contains("Ttree_"+listOfMCs[sample]):
        #            outTree = outFile.Get("Ttree_"+listOfMCs[sample])
        #        else:
    #next do the data files
        #and now grab the WZ systematics, which have to be done separately (but are required for all
        theorySysts = {"matching":"Match","scale":"Scale"}
        plusMinus = {"plus":"Up","minus":"Down"}
        for pm in plusMinus.keys():
            for syst in theorySysts.keys():
                #if not WZ, just clone the tree we already made here.
                if not ("WZ" in sample or ("ttZ" in sample and "scale" in syst)):
                    outTree = outFile.Get("Ttree_"+treeNamePostfixSig+listOfMCs[sample]).Clone("Ttree_"+treeNamePostfixSig+listOfMCs[sample]+"__"+syst+"__"+pm)
                    outTreeSB = 0
                    if useSidebandRegion: outTreeSB = outFile.Get("Ttree_"+treeNamePostfixSB+listOfMCs[sample]).Clone("Ttree_"+treeNamePostfixSB+listOfMCs[sample]+"__"+syst+"__"+pm)
                    outFile.cd()
                    outTree.Write()
                    if useSidebandRegion: outTreeSB.Write()
                    outFile.Write()
                    continue
                else:
                    outTree = TTree("Ttree_"+treeNamePostfixSig+listOfMCs[sample]+"__"+syst+"__"+pm, "Ttree_"+treeNamePostfixSig+listOfMCs[sample]+"__"+syst+"__"+pm)
                    outTreeSB = 0
                    if useSidebandRegion:
                        outTreeSdBnd = TTree("Ttree_"+treeNamePostfixSB+listOfMCs[sample]+"__"+syst+"__"+pm, "Ttree_"+treeNamePostfixSB+listOfMCs[sample]+"__"+syst+"__"+pm)
                    setupBranches(outTree,inputVars)
                    if useSidebandRegion:
                        setupBranches(outTreeSdBnd,inputVars)
                    for channel in channels:
                        inFile = TFile(inputDir+sample+theorySysts[syst]+plusMinus[pm]+channel+"mvaOut.root","READ")
                        tree = inFile.Get("tree")
                        print syst + "__" + pm + ": " + str(tree.GetEntriesFast()),
                        sys.stdout.flush()
                        fillTree(outTree,outTreeSdBnd,inputVars,tree,listOfMCs[sample]+"__"+syst+"__"+pm,channel,jetUnc)
                        inFile.Close()
                    outFile.cd()
                    outFile.Write()
                    outTree.Write()
                    if useSidebandRegion:
                        outTreeSdBnd.Write()
            #Now do the tZq thing.
        #outTree = TTree("Ttree_"+listOfMCs[sample]+"__ttZModel__"+pm, "Ttree_"+listOfMCs[sample]+"__ttZModel__"+pm)
        #for channel in channels:
        #    if not "ttZ" in sample:
        #        #if not one of the samples that has the syst uncerts, 
        #        inFile = TFile(inputDir+sample+channel+"mvaOut.root","READ")
        #        tree = inFile.Get("tree")
        #    else:
        #        inFile = TFile(inputDir+"ttZModel"+channel+"mvaOut.root","READ")
        #        tree = inFile.Get("tree")
        #    fillTree(outTree,inputVars,tree,listOfMCs[sample]+"__ttZModel__plus",channel,jetUnc)
        #    inFile.Close()
        #outFile.cd()
        #outFile.Write()
        #outTree.Write()
            
        outFile.Write()
        outFile.Close()
        print
    chanMap = {"eee":"eeRun2012","eemu":"emuRun2012","emumu":"emuRun2012","mumumu":"mumuRun2012"}

    outChannels = ["DataEG","DataMuEG","DataMu"]
    outChanToData = {}
    outChanToData["DataEG"] = ["eee"]
    outChanToData["DataMuEG"] = ["eemu","emumu"]
    outChanToData["DataMu"] = ["mumumu"]

    for outChan in outChannels:
        print "Data ",outChan
        outTree = TTree("Ttree_"+treeNamePostfixSig+outChan,"Ttree_"+treeNamePostfixSig+outChan)
        outTreeSB = 0
        if useSidebandRegion:
            outTreeSB = TTree("Ttree_"+treeNamePostfixSB+outChan,"Ttree_"+treeNamePostfixSB+outChan)
            setupBranches(outTreeSB,inputVars)
        setupBranches(outTree,inputVars)
        outFile = TFile(outputDir+"histofile_"+outChan+".root","RECREATE")
        for chan in outChanToData[outChan]:
            dataChain = TChain("tree")    
            for run in ["A","B","C","D"]:
                dataChain.Add(inputDir+chanMap[chan]+run+chan+"mvaOut.root")
            fillTree(outTree,outTreeSB, inputVars, dataChain, outChan, chan, 0)
        outFile.cd()
        outFile.Write()
        outTree.Write()
        if useSidebandRegion:
            outTreeSB.Write()
        outFile.Close()

    zEnrichSyst = ["","__zPt__plus","__zPt__minus"]
    for outChan in outChannels:
        print "And finally z-enriched data ",outChan
        outFileZ = TFile(outputDir+"histofile_"+outChan+"Zenriched.root","RECREATE")
        for systPost in zEnrichSyst:
            outTreeZ = TTree("Ttree_"+treeNamePostfixSig+outChan+"Zenriched"+systPost,"Ttree_"+treeNamePostfixSig+outChan+"Zenriched"+systPost)
            setupBranches(outTreeZ,inputVars)
            outTreeZSB = 0
            if useSidebandRegion:
                outTreeZSB = TTree("Ttree_"+treeNamePostfixSB+outChan+"Zenriched"+systPost,"Ttree_"+treeNamePostfixSB+outChan+"Zenriched"+systPost)
                setupBranches(outTreeZSB,inputVars)
            for chan in outChanToData[outChan]:
                dataChainZ = TChain("tree")
                for run in ["A","B","C","D"]:
                    dataChainZ.Add(inputDir+chanMap[chan]+run+chan+"invIsomvaOut.root")
                zPtSyst = 0.
                if "plus" in systPost:
                    zPtSyst = 1.
                if "minus" in systPost:
                    zPtSyst = -1.
                fillTree(outTreeZ,outTreeZSB, inputVars, dataChainZ, outChan, chan, 0, zEnrichWeights[chan],zPtEventWeight = zPtSyst)
            outFileZ.cd()
            outFileZ.Write()
            outTreeZ.Write()
            if useSidebandRegion:
                outTreeZSB.Write()
        outFileZ.Close()        


if __name__ == "__main__":
    main()
                               

