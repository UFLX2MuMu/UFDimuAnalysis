import ROOT


def inv_mass(pt1,pt2,eta1,eta2,phi1,phi2,m1,m2):
  v1 = ROOT.TLorentzVector()
  v1.SetPtEtaPhiM(pt1,eta1,phi1,m1)
  v2 = ROOT.TLorentzVector()
  v2.SetPtEtaPhiM(pt2,eta2,phi2,m2)
  return (v1+v2).M()

def deltaR(ele,muons):
  deltaR=100
  el = ROOT.TLorentzVector(ele.pt,ele.eta,ele.phi,0.)
  for mu in muons:
    mu = ROOT.TLorentzVector(mu.pt,mu.eta,mu.phi,0.)
    deltaR = el.DeltaR(mu)
  return deltaR 

def jetPUid(j_pt, j_eta, puId,  wp):
   '''Function to calculate the pile jet ID as a function of eta and pt'''
   aeta= abs(j_eta);
   pt = j_pt; # no syst
   jet_id_pass = 1;
   if(wp == 100): #Loose 
       if ( aeta< 2.5):
            if (pt >=30 and pt< 50 and puId <-0.89):  return 0;
            if (pt >=10 and pt< 30 and puId <-0.97):  return 0;
       elif (aeta < 2.75):
            if (pt >=30 and pt< 50 and puId <-0.52):  return 0;
            if (pt >=10 and pt< 30 and puId <-0.68):  return 0;
       elif (aeta < 3.00):
            if (pt >=30 and pt< 50 and puId <-0.38):  return 0;
            if (pt >=10 and pt< 30 and puId <-0.53):  return 0;
       elif (aeta < 5.00):
            if (pt >=30 and pt< 50 and puId <-0.30):  return 0;
            if (pt >=10 and pt< 30 and puId <-0.47):  return 0;
   if(wp == 200 ): # Medium
       if ( aeta< 2.5):
            if (pt >=30 and pt< 50 and puId <+0.61):  return 0;
            if (pt >=10 and pt< 30 and puId <+0.18):  return 0;
       elif (aeta < 2.75):
            if (pt >=30 and pt< 50 and puId <-0.35):  return 0;
            if (pt >=10 and pt< 30 and puId <-0.55):  return 0;
       elif (aeta < 3.00):
            if (pt >=30 and pt< 50 and puId <-0.23):  return 0;
            if (pt >=10 and pt< 30 and puId <-0.42):  return 0;
       elif (aeta < 5.00):
            if (pt >=30 and pt< 50 and puId <-0.17):  return 0;
            if (pt >=10 and pt< 30 and puId <-0.36):  return 0;
   if(wp == 300 ): # Tight
       if ( aeta< 2.5):
            if (pt >=30 and pt< 50 and puId <+0.86):  return 0;
            if (pt >=10 and pt< 30 and puId <+0.69):  return 0;
       elif (aeta < 2.75):
            if (pt >=30 and pt< 50 and puId <-0.10):  return 0;
            if (pt >=10 and pt< 30 and puId <-0.35):  return 0;
       elif (aeta < 3.00):
            if (pt >=30 and pt< 50 and puId <-0.05):  return 0;
            if (pt >=10 and pt< 30 and puId <-0.26):  return 0;
       elif (aeta < 5.00):
            if (pt >=30 and pt< 50 and puId <-0.01):  return 0;
            if (pt >=10 and pt< 30 and puId <-0.21):  return 0;
   return jet_id_pass;




infile = ROOT.TFile("../Ntupliser/DiMuons/SingleMu2017B_output_test.root")
tree = infile.Get("dimuons/tree")

#infile = ROOT.TFile("/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/SingleMuon/SingleMu_2017B/180607_161047_tuple_0.root")
#tree = infile.Get("dimuons/tree")

print("Total entries")
print tree.GetEntries()

ele_veto = "Sum$(eles.isMediumID && eles.pt>10 && (abs(eles.eta)<1.4442 || (abs(eles.eta)>1.566 && abs(eles.eta)<2.5)))<1"
jet_selection = "jets.pt>20 && abs(jets.eta)<4.7"
preselection ="Sum$("+jet_selection+")>1 && jets[0].pt>40 && jets[1].pt>20" 
muon_selection = "muons.pt > 20 && abs(muons.eta)<2.4 && muons.isMediumID==1 && muons.relIso < 0.25"
event_selection = "( Sum$(jets.pt>20 && abs(jets.eta)<2.4 && jets.CSV>0.4941)>=0 )  && muons.pt>30 &&  muPairs.mass > 80 && muPairs.mass < 85 && ( muons[0].isHltMatched[2]==1 || muons[0].isHltMatched[3]==1 || muons[1].isHltMatched[2]==1 || muons[1].isHltMatched[3]==1 )"

nev_vbfTightSelection=0
nev_vbfTightSelectionFail=0
nev_ggFTightSelection=0
nev_ggFTightSelectionFail=0
nev_preselectionFail=0
nev_01jetTight=0
nev_01jetLoose=0
nev_preselection=0

nev_prepreselection=0
nev_nJetsAfterSelection=0
nev_bjets=0
nev_leadingJetPt40=0
nev_secondLeadingJetPt20=0 
nev_dijet_mass=0
 
file_vbf_tight= open("vbf_tight.txt","w") 

#gROOT.ProcessLine(
#  "Struct SlimJetInfo{ \
#     jets_	jets_
#     jets.pt	pt[jets_]
#     jets.eta	eta[jets_]
#     jets.phi	phi[jets_]
#     jets.mass	mass[jets_]
#     jets.partonID	partonID[jets_]
#     jets.jecFactor	jecFactor[jets_]
#     jets.jecUnc	jecUnc[jets_]
#     jets.CSV	CSV[jets_]
#     jets.puID	puID[jets_] 
#  };" 
#)

for ev in tree:
  nJetsAfterSelection=0
  dijet_mass=0
  dijet_deta=0
  nbjets=0
  muon_selection=0
  nmupairs=0
  leadingJetPt=0
  secondLeadingJetPt=0
  eleVeto=0
  #muon loop
  for mu in ev.muons:
    if(mu.isGlobal < 1 or mu.isTracker < 1 or mu.isPF < 1):
      continue
    muonPUIso = 1.
    muonPUIso = ( mu.sumChargedHadronPtR04 + max(0., mu.sumNeutralHadronEtR04 + mu.sumPhotonEtR04 - 0.5 * mu.sumPUPtR04 ) ) / mu.pt
    if(  mu.pt>20 and abs(mu.eta)<2.4 and mu.isMediumID==1 and muonPUIso < 0.25 ):
      if( mu.pt > 30 and (mu.isHltMatched[2]==1 or mu.isHltMatched[3]==1) ):
        muon_selection=1 
  #electon loop
  for ele in ev.eles:
    if (eleVeto==1):
      continue
    if ( ele.pt > 10 and (abs(ele.eta)<1.4442 or (abs(ele.eta)>1.566 and abs(ele.eta)<2.5))):
      if( deltaR(ele,ev.muons) < 0.4) :
         eleVeto=1  
  #dimuon loop
  for mupair in ev.muPairs:
    if ( abs(mupair.charge) < 1 ): # and  mupair.mass > 80 and mupair.mass < 85 ): 
      if( ev.muons[mupair.iMu1].isGlobal < 1 or ev.muons[mupair.iMu1].isTracker < 1 or ev.muons[mupair.iMu1].isPF < 1 or ev.muons[mupair.iMu1].pt < 30 ): continue
      if( ev.muons[mupair.iMu2].isGlobal < 1 or ev.muons[mupair.iMu2].isTracker < 1 or ev.muons[mupair.iMu2].isPF < 1 or ev.muons[mupair.iMu2].pt < 20 ): continue
      nmupairs+=1
  if( muon_selection != 1 or nmupairs!=1 ):
    continue
  #jet loop
  for j in range(0,ev.nJets):
   # print('jets {0} {1} {2}'.format(ev.jets[j].pt, ev.jets[j].eta , ev.jets[j].puID))
   # print(jetPUid(ev.jets[j].pt, ev.jets[j].eta , ev.jets[j].puID, 300))
   # jet selection # deltaR between muons and jet is already applied during ntupliser (see JetHelper.cc L179)
    if ( ev.jets[j].pt > 30 and abs(ev.jets[j].eta)<4.7 and jetPUid(ev.jets[j].pt, ev.jets[j].eta , ev.jets[j].puID, 100) > 0. ):
      nJetsAfterSelection+=1
      if( ev.jets[j].pt > leadingJetPt ):
        leadingJetPt = ev.jets[j].pt       
      elif ( ev.jets[j].pt > secondLeadingJetPt ):
          secondLeadingJetPt = ev.jets[j].pt
      # count bjets
      if( ev.jets[j].pt>30 and abs(ev.jets[j].eta)<2.4 and ev.jets[j].CSV>0.4941):
        nbjets+=1
      # build di-jet variables
      for k in range(j+1,ev.nJets):
        if ( ev.jets[k].pt > 30 and abs(ev.jets[k].eta)<4.7 and jetPUid(ev.jets[k].pt, ev.jets[k].eta , ev.jets[k].puID, 100) > 0. ):
          if(dijet_deta < abs(ev.jets[j].eta-ev.jets[k].eta)):
            dijet_deta = abs(ev.jets[j].eta-ev.jets[k].eta)
          if( dijet_mass < inv_mass(ev.jets[j].pt, ev.jets[k].pt, ev.jets[j].eta, ev.jets[k].eta, ev.jets[j].phi, ev.jets[k].phi, ev.jets[j].mass, ev.jets[k].mass)):
            dijet_mass = inv_mass(ev.jets[j].pt, ev.jets[k].pt, ev.jets[j].eta, ev.jets[k].eta, ev.jets[j].phi, ev.jets[k].phi, ev.jets[j].mass, ev.jets[k].mass)
  # preselection
#  if (ev.nJets>2 and ev.jets[0].pt>40 and ev.jets[1].pt>20 and nJetsAfterSelection > 1 and nbjets==0): 
  nev_prepreselection +=1
  if (nJetsAfterSelection > 0): nev_nJetsAfterSelection += 1
  if (nbjets > 0) : nev_bjets +=1
  if (leadingJetPt > 40): nev_leadingJetPt40 +=1
  if (secondLeadingJetPt > 20 ) : nev_secondLeadingJetPt20 +=1 
  if (dijet_mass > 100.): nev_dijet_mass +=1

  if( nJetsAfterSelection > 1 and nbjets==0 and leadingJetPt>40 and secondLeadingJetPt>20 ):
    nev_preselection+=1
    # VBF tight selection
    if(dijet_mass > 650 and dijet_deta > 3.5):  
      nev_vbfTightSelection+=1
      file_vbf_tight.write("{0}\n".format(ev.event.event))
    else:
      nev_vbfTightSelectionFail+=1
      mupt_pass=0
      for mup in ev.muPairs:
        if(mup.pt>50):
          mupt_pass=1 
      # ggF Tight selection
      if(dijet_mass > 250 and mupt_pass>0):
          nev_ggFTightSelection+=1
      # VBF Loose selection
      else:
          nev_ggFTightSelectionFail+=1
  else:
    nev_preselectionFail+=1
    # 0-1 jet Tight
    if(ev.muPairs>25):
      nev_01jetTight+=1
    # 0-1 jet Loose
    else:
      nev_01jetLoose+=1

print("nev_prepreselection = {0}".format(nev_prepreselection))
print("nev_nJetsAfterSelection = {0}".format(nev_nJetsAfterSelection))
print("nev_bjets = {0}".format(nev_bjets))
print("nev_dijet_mass = {0}".format(nev_dijet_mass))
print("nev_leadingJetPt40 = {0}".format(nev_leadingJetPt40))
print("nev_secondLeadingJetPt20 = {0}".format(nev_secondLeadingJetPt20 ))
 


print("nev_preselection = {0}".format(nev_preselection) )
print("nev_vbfTightSelection = {0}".format(nev_vbfTightSelection) )
print("nev_vbfTightSelectioniFail = {0}".format(nev_vbfTightSelectionFail) )
print("nev_ggFTightSelection = {0}".format(nev_ggFTightSelection) )
print("nev_ggFTightSelectionFail = {0}".format(nev_ggFTightSelectionFail))
print("nev_preselectionFail = {0}".format(nev_preselectionFail))
print("nev_01jetTight = {0}".format(nev_01jetTight))
print("nev_01jetLoose = {0}".format(nev_01jetLoose))

file_vbf_tight.close()
