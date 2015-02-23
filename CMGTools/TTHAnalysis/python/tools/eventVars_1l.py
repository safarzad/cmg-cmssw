from CMGTools.TTHAnalysis.treeReAnalyzer import *
from ROOT import TLorentzVector, TVector2, std
#, std
import ROOT
import time
import itertools
import PhysicsTools.Heppy.loadlibs
import array

ROOT.gInterpreter.GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector") #need this to be able to use topness code

mt2wSNT = ROOT.heppy.mt2w_bisect.mt2w()
topness = ROOT.Topness.Topness()

def getPhysObjectArray(j): # https://github.com/HephySusySW/Workspace/blob/72X-master/RA4Analysis/python/mt2w.py
  px = j.pt*cos(j.phi )
  py = j.pt*sin(j.phi )
  pz = j.pt*sinh(j.eta )
  E = sqrt(px*px+py*py+pz*pz) #assuming massless particles...
  return array.array('d', [E, px, py,pz])


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
                         "LSLjetptGT80", "htJet30j", "htJet30ja",
                          "MT2W", "Topness"
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
        pmiss  =array.array('d',[event.met_pt * cos(event.met_phi), event.met_pt * sin(event.met_phi)] )
        
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
        NonBJetCMVAMedium30 = []
        for i,j in enumerate(centralJet30):
            if j.btagCMVA>0.732:
                BJetCMVAMedium30.append(j)
                BJetCMVAMedium30idx.append(centralJet30idx[i])
            else:
                NonBJetCMVAMedium30.append(j)
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
            if ret['nTightLeps25']>=1:
                mLepBJet = (jets[idxMinDPhiBMET].p4() + tightLeps25[0].p4()).M()
                mTLepMET = mt_2(tightLeps25[0].p4(),metp4)
        else:
            ret["mTClBPlusMET"] = -999

        ret["mTBJetMET"] = mTBJetMET
        ret["mTLepMET"]  = mTLepMET
        ret["mLepBJet"]  = mLepBJet
        

        # projection of MET along (MET + lepton + (closest (to MET) BJet)) sum; needs to be double-checked...
        METtoTopProjection = -999
        if(idxMinDPhiBMET>=0 and ret['nTightLeps25']>=1) :
            metV2  = ROOT.TVector2(0,0)
            lepV2  = ROOT.TVector2(0,0)
            bJetV2 = ROOT.TVector2(0,0)
            metV2   .SetMagPhi(event.met_pt, event.met_phi)
            lepV2   .SetMagPhi(tightLeps25[0].pt, tightLeps25[0].phi)
            bJetV2  .SetMagPhi(jets[idxMinDPhiBMET].pt, jets[idxMinDPhiBMET].phi)
            METtoTopProjection = (metV2*(metV2+lepV2+bJetV2))/(metV2*metV2)
        ret["METtoTopProjection"]  = METtoTopProjection


        #add topness and mt2W-variable (timing issue with topness: slows down the friend tree production by a factor of ~3)
        ret['Topness']=-999
        mt2w_values=[]
        if ret['nTightLeps25']>=1: #topness and mt2w only make sense for 
            lep = getPhysObjectArray(tightLeps25[0])
            if ret['nBJetCMVAMedium30']==0 and ret['nCentralJet30']>=3: #All combinations from the highest three light (or b-) jets
                consideredJets = [ getPhysObjectArray(jet) for jet in NonBJetCMVAMedium30[:3] ] # only throw arrays into the permutation business
                ftPerms = itertools.permutations(consideredJets, 2)
                for perm in ftPerms:
                    mt2wSNT.set_momenta(lep, perm[0], perm[1], pmiss)
                    mt2w_values.append(mt2wSNT.get_mt2w())
            elif ret['nBJetCMVAMedium30']==1 and ret['nCentralJet30']>=2: #All combinations from one b and the highest two light jets
                consideredJets = [ getPhysObjectArray(jet) for jet in NonBJetCMVAMedium30[:2] ] # only throw arrays into the permutation business
                consideredJets.append(getPhysObjectArray(BJetCMVAMedium30[0]))
                ftPerms = itertools.permutations(consideredJets, 2)
                for perm in ftPerms:
                    mt2wSNT.set_momenta(lep, perm[0], perm[1], pmiss)
                    mt2w_values.append(mt2wSNT.get_mt2w())
            elif ret['nBJetCMVAMedium30']==2: 
                consideredJets = [ getPhysObjectArray(jet) for jet in BJetCMVAMedium30[:2] ] # only throw arrays into the permutation business
                ftPerms = itertools.permutations(consideredJets, 2)
                for perm in ftPerms:
                    mt2wSNT.set_momenta(lep, perm[0], perm[1], pmiss)
                    mt2w_values.append(mt2wSNT.get_mt2w())
            elif ret['nBJetCMVAMedium30']>=3: #All combinations from the highest three b jets
                consideredJets = [ getPhysObjectArray(jet) for jet in BJetCMVAMedium30[:3] ] # only throw arrays into the permutation business
                ftPerms = itertools.permutations(consideredJets, 2)
                for perm in ftPerms:
                    mt2wSNT.set_momenta(lep, perm[0], perm[1], pmiss)
                    mt2w_values.append(mt2wSNT.get_mt2w())

            p4_jets = std.vector(TLorentzVector)();
            bdisc_jets = std.vector('float')();

            for jet in centralJet30:
                jetTLorentz = ROOT.TLorentzVector(0,0,0,0)
                jetTLorentz.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass)
                p4_jets.push_back(jetTLorentz)
                bdisc_jets.push_back(jet.btagCMVA)

            lepTLorentz = ROOT.TLorentzVector(0,0,0,0)
            lepTLorentz.SetPtEtaPhiM(tightLeps25[0].pt, tightLeps25[0].eta, tightLeps25[0].phi, tightLeps25[0].mass)

            if ret['nCentralJet30']>=3: # does not seem to work for njet =3 ??! # need to edit btag working point in the code...!! did not quickly find a twiki with official phys14 cmva working points
                tempTopness = topness.GetTopness(p4_jets,bdisc_jets,lepTLorentz,metp4) #this is really slow!
                if tempTopness <=0:
                    print tempTopness, "this will fail"
                ret['Topness'] = log(tempTopness) #this is really slow!
        if len(mt2w_values)>0:
            ret["MT2W"]=min(mt2w_values)
        else:
            ret["MT2W"]=-999
            
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

        
