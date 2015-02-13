from CMGTools.TTHAnalysis.treeReAnalyzer import *
from ROOT import TLorentzVector, TVector2
import ROOT
import time
import itertools

def mt_2(p4one, p4two):
    return sqrt(2*p4one.Pt()*p4two.Pt()*(1-cos(p4one.Phi()-p4two.Phi())))


class EventVars1L:
    def __init__(self):
        self.branches = [ "METCopyPt", "DeltaPhiLepW", "minDPhiBMET", "idxMinDPhiBMET", "mTClBPlusMET", "mTBJetMET", "mTLepMET", "mLepBJet", "METtoTopProjection",
                         ("nTightLeps25","I"),("nTightMu25","I"),("nTightEl25","I"),("nVetoLeps10","I"),
                         ("tightLeps25idx","I",10,"nTightLeps25"),("tightMu25idx","I",10,"nTightMu25"),("tightEl25idx","I",10,"nTightEl25"),
                         ("vetoLeps10idx","I",10,"nVetoLeps10"),
                         ("nCentralJet30","I"),("centralJet30idx","I",100,"nCentralJet30"),
                         ("nBJetCMVAMedium30","I"),("BJetCMVAMedium30idx","I",50,"nBJetCMVAMedium30"),
                         "nGoodBJets_allJets", "nGoodBJets",
                         "LSLjetptGT80", "htJet30j", "htJet30ja"
                         ]

    
    def listBranches(self):
        return self.branches[:]
    

    def __call__(self,event):

        # make python lists as Collection does not support indexing in slices
        leps = [l for l in Collection(event,"LepGood","nLepGood")]
        jets = [j for j in Collection(event,"Jet","nJet")]

        njet = len(jets); nlep = len(leps)
        
        metp4 = ROOT.TLorentzVector(0,0,0,0)
        metp4.SetPtEtaPhiM(event.met_pt,event.met_eta,event.met_phi,event.met_mass)
        
        #isolation criteria as defined for PHYS14 1l synchronisation exercise
        ele_relisoCut = 0.14
        muo_relisoCut = 0.12
        #ele tight id --> >2
        #muo tight id ==1
        centralEta = 2.4
        
        tightLeps25 = []
        tightLeps25idx = []
        tightMu25 = []
        tightMu25idx = []
        tightEl25 = []
        tightEl25idx = []
        vetoLeps10 = []
        vetoLeps10idx = []
        for i,l in enumerate(leps):
            tightMu = l.pt>25 and l.relIso03<muo_relisoCut and abs(l.pdgId)==13 and l.tightId==1
            tightEl = l.pt>25 and l.relIso03<ele_relisoCut and abs(l.pdgId)==11 and l.tightId >2 
            if tightMu:
                tightMu25.append(l); tightMu25idx.append(i)
                tightLeps25.append(l); tightLeps25idx.append(i)
            elif tightEl:
                tightEl25.append(l); tightEl25idx.append(i)
                tightLeps25.append(l); tightLeps25idx.append(i)
            elif l.pt>10:
                vetoLeps10.append(l); vetoLeps10idx.append(i)

        ret = { 'nTightLeps25'   : len(tightLeps25) } #initialize the dictionary with a first entry
        ret['nTightMu25']  = len(tightMu25)
        ret['nTightEl25']  = len(tightEl25)
        ret['nVetoLeps10'] = len(vetoLeps10)
        
        ret['tightLeps25idx'] = tightLeps25idx
        ret['tightMu25idx']   = tightMu25idx
        ret['tightEl25idx']   = tightEl25idx
        ret['vetoLeps10idx']  = vetoLeps10idx
        

        centralJet30 = []
        centralJet30idx = []
        for i,j in enumerate(jets):
            if j.pt>30 and abs(j.eta)<centralEta:
                centralJet30.append(j)
                centralJet30idx.append(i)
        
        ret['nCentralJet30']   = len(centralJet30)
        ret['centralJet30idx'] = centralJet30idx
        
        ret['LSLjetptGT80'] = 1 if sum([j.pt>80 for j in centralJet30])>=2 else 0

        ret['htJet30j']  = sum([j.pt for j in centralJet30])
        ret['htJet30ja'] = sum([j.pt for j in jets if j.pt>30])

        BJetCMVAMedium30 = []
        BJetCMVAMedium30idx = []
        for i,j in enumerate(centralJet30):
            if j.btagCMVA>0.732:
                BJetCMVAMedium30.append(j)
                BJetCMVAMedium30idx.append(centralJet30idx[i])

        ret['nBJetCMVAMedium30']    = len(BJetCMVAMedium30)
        ret['BJetCMVAMedium30idx']  = BJetCMVAMedium30idx

        ret['nGoodBJets']    = sum([j.btagCMVA>0.732 for j in centralJet30])
        ret['nGoodBJets_allJets']    = sum([j.btagCMVA>0.732 and j.pt>30 and abs(j.eta)<centralEta for j in jets]) # where is the working point defined?


        #plain copy of MET pt (just as an example and cross-check for proper friend tree production)
        ret["METCopyPt"] = metp4.Pt()
        
        # deltaPhi between the (single) lepton and the reconstructed W (lep + MET)
        dPhiLepW = -999 # set default value to -999 to spot "empty" entries
        if ret['nTightLeps25']>=1:
            recoWp4 =  tightLeps25[0].p4() + metp4
            dPhiLepW = tightLeps25[0].p4().DeltaPhi(recoWp4)
        ret["DeltaPhiLepW"] = dPhiLepW

        ##################################################################
        # The following variables need to be double-checked for validity #
        ##################################################################

        # min deltaPhi between a (CMVA) b-jet and MET; needs to be double-checked
        minDPhiBMET    = 100
        idxMinDPhiBMET = -999
        for i, jet in enumerate(jets):
            if jet.btagCMVA>0.732:
                dPhiBMET = abs(jet.p4().DeltaPhi(metp4))
                if dPhiBMET<minDPhiBMET:
                    minDPhiBMET=dPhiBMET
                    idxMinDPhiBMET = i

        ret["idxMinDPhiBMET"] = idxMinDPhiBMET
        ret["minDPhiBMET"] = minDPhiBMET


        # transverse mass of (closest (to MET) BJet, MET), (closest (to MET) BJet, lepton), 
        # mass of (closest (to MET) BJet, lepton); need to be double-checked
        mTBJetMET      = -999
        mTLepMET       = -999
        mLepBJet       = -999
        if(idxMinDPhiBMET>=0):
            SumMetClosestBJet = jets[idxMinDPhiBMET].p4() + metp4
            ret["mTClBPlusMET"] = SumMetClosestBJet.Mt()
            mTBJetMET = mt_2(jets[idxMinDPhiBMET].p4(),metp4)
            if nlep>=1:
                mLepBJet = (jets[idxMinDPhiBMET].p4() + leps[0].p4()).M()
                mTLepMET = mt_2(leps[0].p4(),metp4)
        else:
            ret["mTClBPlusMET"] = -999

        ret["mTBJetMET"] = mTBJetMET
        ret["mTLepMET"]  = mTLepMET
        ret["mLepBJet"]  = mLepBJet
        

        # projection of MET along (MET + lepton + (closest (to MET) BJet)) sum; needs to be double-checked...
        METtoTopProjection = -999
        if(idxMinDPhiBMET>=0 and nlep >=1) :
            metV2  = ROOT.TVector2(0,0)
            lepV2  = ROOT.TVector2(0,0)
            bJetV2 = ROOT.TVector2(0,0)
            metV2   .SetMagPhi(event.met_pt, event.met_phi)
            lepV2   .SetMagPhi(leps[0].pt, leps[0].phi)
            bJetV2  .SetMagPhi(jets[idxMinDPhiBMET].pt, jets[idxMinDPhiBMET].phi)
            METtoTopProjection = (metV2*(metV2+lepV2+bJetV2))/(metV2*metV2)
        ret["METtoTopProjection"]  = METtoTopProjection

            
        return ret

if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree") 
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.sf = EventVars1L()
        def analyze(self,ev):
            print "\nrun %6d lumi %4d event %d: leps %d" % (ev.run, ev.lumi, ev.evt, ev.nLepGood)
            print self.sf(ev)
    el = EventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = 50)

        
