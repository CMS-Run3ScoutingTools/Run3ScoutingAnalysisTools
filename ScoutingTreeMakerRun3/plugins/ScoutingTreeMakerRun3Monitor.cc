// -*- C++ -*-
//
// Package:    Run3ScoutingAnalysisTools/ScoutingTreeMakerRun3Monitor
// Class:      ScoutingTreeMakerRun3Monitor
//
/**\class ScoutingTreeMakerRun3Monitor ScoutingTreeMakerRun3Monitor.cc Run3ScoutingAnalysisTools/ScoutingTreeMakerRun3/plugins/ScoutingTreeMakerRun3Monitor.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  David Sperka
//         Created:  Sat, 11 Feb 2023 14:15:08 GMT
//
//

// system include files
#include <memory>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

#include "DataFormats/PatCandidates/interface/Muon.h"
//
// class declaration
//

class ScoutingTreeMakerRun3Monitor : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit ScoutingTreeMakerRun3Monitor(const edm::ParameterSet&);
  ~ScoutingTreeMakerRun3Monitor() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>                 triggerResultsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon> >      muonsToken;
  const edm::EDGetTokenT<std::vector<pat::Muon> >             offlineMuonsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> >  electronsToken;
  //const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >  verticesToken;
  const edm::EDGetTokenT<double>                              rhoToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> >    photonsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle> >  pfcandsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPFJet> >     pfjetsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingTrack> >     tracksToken;
  
  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;
  
  bool doL1;
  triggerExpression::Data triggerCache_;
  
  edm::InputTag                algInputTag_;
  edm::InputTag                extInputTag_;
  edm::EDGetToken              algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;
  std::vector<std::string>     l1MonitorSeeds_;
  std::vector<bool>            l1Result_;
  std::vector<bool>            l1Result_mon_;
  
  TTree* tree;
  float mass;
  float pt;
  float dr;
  float pt1, pt2;
  float eta1, eta2;
  float phi1, phi2;
  int   id1, id2;
  float rho;
  float pfIso;
  float dr_matching_1;
  float dr_matching_2;
  int   nScoutingMuons;
  int   nScoutingMuons_matched;
  float pt1_scout, pt2_scout;
  float mass_scout;
  float ptmm_scout;


  std::vector<bool> pfIsoWP;
  /* std::vector<std::string> pfIsoWP;
  pfIsoWP.push_back("PFIsoLoose");
  pfIsoWP.push_back("PFIsoMedium");
  pfIsoWP.push_back("PFIsoTight");
  pfIsoWP.push_back("PFIsoVeryTight"); 
  std::map<std::string, bool> pfIsoWPMap; */
  std::vector<bool> muID;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ScoutingTreeMakerRun3Monitor::ScoutingTreeMakerRun3Monitor(const edm::ParameterSet& iConfig):
    triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
    triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),
    muonsToken               (consumes<std::vector<Run3ScoutingMuon> >         (iConfig.getParameter<edm::InputTag>("muons"))),
    offlineMuonsToken        (consumes<std::vector<pat::Muon> >                (iConfig.getUntrackedParameter<edm::InputTag>("offlineMuons"))),
    electronsToken           (consumes<std::vector<Run3ScoutingElectron> >     (iConfig.getParameter<edm::InputTag>("electrons"))),
    //verticesToken          (consumes<std::vector<Run3ScoutingVertex> >       (iConfig.getParameter<edm::InputTag>("vertices"))),
    rhoToken                 (consumes<double>                                 (iConfig.getParameter<edm::InputTag>("rho"))), 
    photonsToken             (consumes<std::vector<Run3ScoutingPhoton> >       (iConfig.getParameter<edm::InputTag>("photons"))),
    pfcandsToken             (consumes<std::vector<Run3ScoutingParticle> >     (iConfig.getParameter<edm::InputTag>("pfcands"))),
    pfjetsToken              (consumes<std::vector<Run3ScoutingPFJet> >        (iConfig.getParameter<edm::InputTag>("pfjets"))),
    tracksToken              (consumes<std::vector<Run3ScoutingTrack> >        (iConfig.getParameter<edm::InputTag>("tracks"))),
    doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1") : false)
{
    usesResource("TFileService");
    if (doL1) {
        algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
        extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
        algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
        l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
        l1MonitorSeeds_ = iConfig.getParameter<std::vector<std::string> >("l1MonitorSeeds");
        l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
    }
    else {
        l1Seeds_ = std::vector<std::string>();
        l1MonitorSeeds_ = std::vector<std::string>();
        l1GtUtils_ = 0;
    }
}

ScoutingTreeMakerRun3Monitor::~ScoutingTreeMakerRun3Monitor() {
}

//
// member functions
//

// ------------ method called for each event  ------------
void ScoutingTreeMakerRun3Monitor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;

  /*
  Handle<vector<ScoutingVertex> > verticesH;
  iEvent.getByToken(verticesToken, verticesH);
  */

  Handle<vector<pat::Muon> > offlineMuonsH;
  iEvent.getByToken(offlineMuonsToken, offlineMuonsH);
  
  if (offlineMuonsH->size()<2) return;

  int nOfflineMuons=0;
  vector<int> idx;
  int j=0;

  for (auto muons_iter = offlineMuonsH->begin(); muons_iter != offlineMuonsH->end(); ++muons_iter) {
    if (muons_iter->pt()>3 &&  abs(muons_iter->eta())<1.9) {
    //if (muons_iter->pt()>3 && muons_iter->isMediumMuon() && abs(muons_iter->eta())<1.9) {

      // Compute the pfIso for the muon. Note: PUCorr = 0.5*muons_iter->puChargedHadronIso()
      // -----------------------------------------------------------------------------------
      pfIso = (muons_iter->chargedHadronIso() + std::max(muons_iter->photonIso() + muons_iter->neutralHadronIso() - 0.5*muons_iter->puChargedHadronIso(), 0.0))/muons_iter->pt();
      //std::cout << "----- pfIso: " << pfIso << endl;

      // Filling a vector of boolean to save the decision for various WP:
      // 0 PFIsoLoose, 1 PFIsoMedium, 2 PFIsoTight, 3 PFIsoVeryTight
      // ----------------------------------------------------------------
      //pfIsoWP.push_back(muons_iter->PFIsoLoose);
      //pfIsoWP.push_back(muons_iter->PFIsoMedium);
      //pfIsoWP.push_back(muons_iter->PFIsoTight);
      //pfIsoWP.push_back(muons_iter->PFIsoVeryTight);

      // Saving Muon ID info for various WP:
      // 0 isLooseMuon, 1 isMediumMuon, 2 isTightMuon
      // ----------------------------------------------------------------
      muID.push_back(muons_iter->isLooseMuon());
      muID.push_back(muons_iter->isMediumMuon());
      //muID.push_back(muons_iter->isTightMuon());

      //std::cout << "----- ID Loose: "  << muons_iter->isLooseMuon() << endl;
      //std::cout << "----- ID Medium: " << muons_iter->isMediumMuon() << endl;
      //std::cout << "----- ID Tight: "  << muons_iter->isTightMuon() << endl;

      //cout<<"offline muon pt: "<<muons_iter->pt()<<" id: "<<muons_iter->isMediumMuon()
      //                         <<" eta: "<<muons_iter->eta()<<" pfIso: " << pfIso <<endl;

      nOfflineMuons+=1;
      idx.push_back(j);
    }
    j+=1;
  }
  
  if (idx.size()<2) {/*cout<<"failed offline muons"<<endl;*/ return;}
  //cout << "idx = " << idx.size() << endl;
  edm::Handle<edm::TriggerResults> triggerResultsH;
  iEvent.getByToken(triggerResultsToken, triggerResultsH);

  bool passDST=false;
  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
      if (triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) passDST=true;
  }

  // Scouting path true => scouting reconstruction enabled
  if (!passDST) {/*cout<<"failed DST"<<endl;*/ return;}

  l1Result_.clear();

  // -------------------------------------------
  // NOTE: No selection on the Monitoring seeds
  // -------------------------------------------
  //bool passMonitor=false; 
  //bool passScouting=false; 
  if (doL1) {
      l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
      /*for(unsigned int r = 0; r<100; r++){
        string name ("empty");
                bool algoName_ = false;
                        algoName_ = l1GtUtils_->getAlgNameFromBit(i,name);
                        cout << "getAlgNameFromBit = " << algoName_  << endl;
                        cout << "L1 bit number = " << i << " ; L1 bit name = " << name << endl;
                        }*/
      for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
          bool l1htbit = 0;
          l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
          l1Result_.push_back( l1htbit );
          //if (l1htbit) passScouting=true; // No selection on the Scouting seeds
      }
      for( unsigned int iseed = 0; iseed < l1MonitorSeeds_.size(); iseed++ ) {
          bool l1htbit = 0;
          l1GtUtils_->getFinalDecisionByName(string(l1MonitorSeeds_[iseed]), l1htbit);
	  l1Result_mon_.push_back( l1htbit );
          //if (l1htbit) passMonitor=true; // No selection on the Monitoring seeds
      }
  }
  //if (!passMonitor) {/*cout<<"failed L1 seed"<<endl;*/ return;} // No selection on the Monitoring seeds
  //if (!passScouting) {/*cout<<"failed L1 seed"<<endl;*/ return;} // No selection on the Scouting seeds
    
  Handle<double> rhoH;
  iEvent.getByToken(rhoToken, rhoH);

  //cout << "------------------------------------------" << endl;
  //cout<<"...PASSED ALL SELECTIONS!"<<endl;
  //cout << "------------------------------------------" << endl;

  Handle<vector<Run3ScoutingMuon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);

  nScoutingMuons=0;
  nScoutingMuons_matched=0;
  vector<int> idx_scout;
  int s=0;
  for (auto muons_scout_iter = muonsH->begin(); muons_scout_iter != muonsH->end(); ++muons_scout_iter) {
    //cout<<"scouting muon pt: "<<muons_scout_iter->pt()<<" eta: "<<muons_scout_iter->eta()<<endl;                                 
    if (muons_scout_iter->pt()>3 and abs(muons_scout_iter->eta())<1.9) {
      nScoutingMuons+=1;
      idx_scout.push_back(s);
      
    }
    s+=1;
  }
  if (idx_scout.size()<2) {/*cout<<"failed scouting muons"<<endl;*/ return;}
  //cout << "idx_scout = " << idx_scout.size() << "  and nScoutingMuons = " << nScoutingMuons << endl;

  // -----------------------------------------
  // Checking the muon matching based on dR //                                                                                      
  // -----------------------------------------
  if (muonsH->at(idx_scout[0]).pt() > muonsH->at(idx_scout[1]).pt()) {
    dr_matching_1 = sqrt( (offlineMuonsH->at(idx[0]).eta() - muonsH->at(idx_scout[0]).eta())*(offlineMuonsH->at(idx[0]).eta() - muonsH->at(idx_scout[0]).eta()) + (offlineMuonsH->at(idx[0]).phi() - muonsH->at(idx_scout[0]).phi())*(offlineMuonsH->at(idx[0]).phi() - muonsH->at(idx_scout[0]).phi()) );
    dr_matching_2 = sqrt( (offlineMuonsH->at(idx[1]).eta() - muonsH->at(idx_scout[1]).eta())*(offlineMuonsH->at(idx[1]).eta() - muonsH->at(idx_scout[1]).eta()) + (offlineMuonsH->at(idx[1]).phi() - muonsH->at(idx_scout[1]).phi())*(offlineMuonsH->at(idx[1]).phi() - muonsH->at(idx_scout[1]).phi()) );
    //cout << "dr_matching_1 = " << dr_matching_1 << " with offline mu eta: " << offlineMuonsH->at(idx[0]).eta() << " and scout mu eta: " << muonsH->at(idx_scout[0]).eta() << endl;
    //cout << "dr_matching_2 = " << dr_matching_2 << " with offline mu eta: " << offlineMuonsH->at(idx[1]).eta() << " and scout mu eta: " << muonsH->at(idx_scout[1]).eta() << endl;

  }
  else if (muonsH->at(idx_scout[0]).pt() < muonsH->at(idx_scout[1]).pt()) {
    dr_matching_1 = sqrt( (offlineMuonsH->at(idx[0]).eta() - muonsH->at(idx_scout[1]).eta())*(offlineMuonsH->at(idx[0]).eta() - muonsH->at(idx_scout[1]).eta()) + (offlineMuonsH->at(idx[0]).phi() - muonsH->at(idx_scout[1]).phi())*(offlineMuonsH->at(idx[0]).phi() - muonsH->at(idx_scout[1]).phi()) );
    dr_matching_2 = sqrt( (offlineMuonsH->at(idx[1]).eta() - muonsH->at(idx_scout[0]).eta())*(offlineMuonsH->at(idx[1]).eta() - muonsH->at(idx_scout[0]).eta()) + (offlineMuonsH->at(idx[1]).phi() - muonsH->at(idx_scout[0]).phi())*(offlineMuonsH->at(idx[1]).phi() - muonsH->at(idx_scout[0]).phi()) );
    //cout << "dr_matching_1 = " << dr_matching_1 << " with offline mu eta: " << offlineMuonsH->at(idx[0]).eta() << " and scout mu eta: " << muonsH->at(idx_scout[1]).eta() << endl;
    //cout << "dr_matching_2 = " << dr_matching_2 << " with offline mu eta: " << offlineMuonsH->at(idx[1]).eta() << " and scout mu eta: " << muonsH->at(idx_scout[0]).eta() << endl;
  }
  else {
    cout << "NOTE: dR for the matching has not been computed!" << endl;
  }

  //cout << "Filling offline muon quantities!" << endl;
    pt1=offlineMuonsH->at(idx[0]).pt();
    pt2=offlineMuonsH->at(idx[1]).pt();
    eta1=offlineMuonsH->at(idx[0]).eta();
    eta2=offlineMuonsH->at(idx[1]).eta();
    phi1=offlineMuonsH->at(idx[0]).phi();
    phi2=offlineMuonsH->at(idx[1]).phi();
    id1=offlineMuonsH->at(idx[0]).pdgId();
    id2=offlineMuonsH->at(idx[1]).pdgId();
    
    TLorentzVector mu1;
    mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.105658);
    TLorentzVector mu2;
    mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.105658);
    TLorentzVector dimu = mu1+mu2;
    mass=dimu.M();
    pt=dimu.Pt();
    dr=mu1.DeltaR(mu2);
    
    rho=*rhoH;
    
    //cout << "Filling scouting muon quantities!" << endl;
    TLorentzVector mu1_scout;
    TLorentzVector mu2_scout;
    
    //cout << "%% muonsH->at(idx_scout[0]).pt() = " << muonsH->at(idx_scout[0]).pt() << endl;
    //cout << "%% muonsH->at(idx_scout[1]).pt() = " << muonsH->at(idx_scout[1]).pt() << endl;
    
    if (muonsH->at(idx_scout[0]).pt() > muonsH->at(idx_scout[1]).pt()) {
      //cout << "%%% Entering the first if"  << endl;
      pt1_scout=muonsH->at(idx_scout[0]).pt();
      pt2_scout=muonsH->at(idx_scout[1]).pt();
      mu1_scout.SetPtEtaPhiM(pt1_scout,muonsH->at(idx_scout[0]).eta(),muonsH->at(idx_scout[0]).phi(),0.105658);
      mu2_scout.SetPtEtaPhiM(pt2_scout,muonsH->at(idx_scout[1]).eta(),muonsH->at(idx_scout[1]).phi(),0.105658);
      //cout << "%%% pt1_scout = " << muonsH->at(idx_scout[0]).pt() << endl;
      //cout << "%%% pt2_scout = " << muonsH->at(idx_scout[1]).pt() << endl;
    }
    else if (muonsH->at(idx_scout[0]).pt() < muonsH->at(idx_scout[1]).pt()) {
      //cout << "%%% Entering the second if"  << endl;
      pt1_scout=muonsH->at(idx_scout[1]).pt();
      pt2_scout=muonsH->at(idx_scout[0]).pt();
      mu1_scout.SetPtEtaPhiM(pt1_scout,muonsH->at(idx_scout[1]).eta(),muonsH->at(idx_scout[1]).phi(),0.105658);
      mu2_scout.SetPtEtaPhiM(pt2_scout,muonsH->at(idx_scout[0]).eta(),muonsH->at(idx_scout[0]).phi(),0.105658);
      //cout << "%%% pt1_scout = " << muonsH->at(idx_scout[1]).pt() << endl;
      //cout << "%%% pt2_scout = " << muonsH->at(idx_scout[0]).pt() << endl;
    }
    
    TLorentzVector dimu_scout = mu1_scout+mu2_scout;
    mass_scout=dimu_scout.M();
    ptmm_scout=dimu_scout.Pt();
    
    if (dr_matching_1 < 0.2 and dr_matching_2 < 0.2){
      nScoutingMuons_matched+=2;
      tree->Fill();
      return;
    }

    /*if (dr_matching_1 < 0.2 and dr_matching_2 < 0.2){
    nScoutingMuons_matched+=2;
    cout << "------------------------------------------" << endl;
    cout<<"...MATCHING!"<<endl;
    cout << "------------------------------------------" << endl;
    tree->Fill();
    //idx.clear();
    //idx_scout.clear();
    }*/
    //###### tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void ScoutingTreeMakerRun3Monitor::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"      , "tree");
    tree->Branch("mass"                , &mass                         , "mass/F");
    tree->Branch("pt"                  , &pt                           , "pt/F");
    tree->Branch("dr"                  , &dr                           , "dr/F");
    tree->Branch("pt1"                 , &pt1                          , "pt1/F");
    tree->Branch("pt2"                 , &pt2                          , "pt2/F");
    tree->Branch("eta1"                , &eta1                         , "eta1/F");
    tree->Branch("eta2"                , &eta2                         , "eta2/F");
    tree->Branch("id1"                 , &id1                          , "id1/I");
    tree->Branch("id2"                 , &id2                          , "id2/I");
    tree->Branch("rho"                 , &rho                          , "rho/F");
    tree->Branch("pfIso"               , &pfIso                        , "pfIso/F");
    tree->Branch("muID", "std::vector<bool>"                           ,&muID, 32000, 0);
    tree->Branch("l1Result", "std::vector<bool>"                       ,&l1Result_, 32000, 0);
    tree->Branch("l1Result_mon", "std::vector<bool>"                   ,&l1Result_mon_, 32000, 0);
    tree->Branch("nScoutingMuons", &nScoutingMuons                     , "nScoutingMuons/I");
    //tree->Branch("nScoutingMuons_matched", &nScoutingMuons_matched     , "nScoutingMuons_matched/I");                                 
    tree->Branch("pt1_scout"           , &pt1_scout                    , "pt1_scout/F");                                    
    tree->Branch("pt2_scout"           , &pt2_scout                    , "pt2_scout/F");                                      
    tree->Branch("mass_scout"          , &mass_scout                   , "mass_scout/F");                                          
    tree->Branch("ptmm_scout"          , &ptmm_scout                   , "ptmm_scout/F");                           

}

// ------------ method called once each job just after ending the event loop  ------------
void ScoutingTreeMakerRun3Monitor::endJob() {
  // please remove this method if not needed
}

void ScoutingTreeMakerRun3Monitor::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    // HLT paths
    triggerPathsVector.push_back("DST_Run3_PFScoutingPixelTracking_v*");

    HLTConfigProvider hltConfig;
    bool changedConfig = false;
    hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);

    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        triggerPathsMap[triggerPathsVector[i]] = -1;
    }

    for(size_t i = 0; i < triggerPathsVector.size(); i++){
        TPRegexp pattern(triggerPathsVector[i]);
        for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
            std::string pathName = hltConfig.triggerNames()[j];
            if(TString(pathName).Contains(pattern)){
                triggerPathsMap[triggerPathsVector[i]] = j;
            }
        }
    }
}


void ScoutingTreeMakerRun3Monitor::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ScoutingTreeMakerRun3Monitor::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void ScoutingTreeMakerRun3Monitor::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingTreeMakerRun3Monitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingTreeMakerRun3Monitor);
