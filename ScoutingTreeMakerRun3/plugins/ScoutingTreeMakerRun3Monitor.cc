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
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >    primaryVerticesToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >    verticesToken;
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
  //Defining offline variables
  int   nScoutingMuons;
  float pt1, pt2;
  float eta1, eta2;
  float phi1, phi2;
  int   id1, id2;
  float rho;
  float pfIso1;
  float pfIso2;
  std::vector<bool> mu1_ID;
  std::vector<bool> mu2_ID;
  float mass;
  float ptmm;
  float drmm;
  //Defining scouting variables: muons
  int   nOfflineMuons;
  float pt1_scout, pt2_scout;
  float mass_scout;
  float ptmm_scout;
  float drmm_scout;
  //Defining scouting variables: vertices
  int nvtx;
  int ndvtx;
  bool vtxMatch;
  float avgPrimaryX;
  float avgPrimaryY;
  float lxy;
  float vtxErr;
  float IP;
  std::vector<float> vtxX;
  std::vector<float> vtxY;
  std::vector<float> vtxXError;
  std::vector<float> vtxYError;

  //Defining matching variables
  float dr_matching_1;
  float dr_matching_2;
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
  primaryVerticesToken     (consumes<std::vector<Run3ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  verticesToken            (consumes<std::vector<Run3ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("displacedVertices"))),
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


  // ------------- //
  // OFFLINE MUONS //
  // ------------- //
  Handle<vector<pat::Muon> > offlineMuonsH;
  iEvent.getByToken(offlineMuonsToken, offlineMuonsH);
  
  if (offlineMuonsH->size()<2) return;

  nOfflineMuons=0;
  vector<int> idx;
  int j=0;

  // Minimal selection on the offline muons
  // --------------------------------------
  for (auto muons_iter = offlineMuonsH->begin(); muons_iter != offlineMuonsH->end(); ++muons_iter) {
    if (muons_iter->pt()>3 &&  abs(muons_iter->eta())<1.9) { //removing the MediumMuon ID requirement, saved later
      nOfflineMuons+=1;
      idx.push_back(j);
    }
    j+=1;
  }
  
  // Requiring exactly two offline muons in the event
  // ------------------------------------------------
  if (!(idx.size()==2)) {/*cout<<"failed offline muons"<<endl;*/ return;}

  // Checking opposite charge                                                                                     
  // ------------------------
  int checkOSCharge = offlineMuonsH->at(idx[0]).charge() * offlineMuonsH->at(idx[1]).charge();
  //cout<<"Check OS requirement for offline muons = " << checkOSCharge << endl; 
  if (checkOSCharge > 0){/*cout<<"not OS pair"<<endl;*/return;}
  //cout<<"Check OS requirement for offline muons = PASSED!" << endl; 

  edm::Handle<edm::TriggerResults> triggerResultsH;
  iEvent.getByToken(triggerResultsToken, triggerResultsH);

  bool passDST=false;
  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
      if (triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) passDST=true;
  }

  // Scouting path true => scouting reconstruction enabled
  // -----------------------------------------------------
  if (!passDST) {/*cout<<"failed DST"<<endl;*/ return;}

  l1Result_.clear();

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
          //if (l1htbit) passScouting=true; // No selection on the Scouting seeds, selection stores
      }
      for( unsigned int iseed = 0; iseed < l1MonitorSeeds_.size(); iseed++ ) {
          bool l1htbit = 0;
          l1GtUtils_->getFinalDecisionByName(string(l1MonitorSeeds_[iseed]), l1htbit);
	  l1Result_mon_.push_back( l1htbit );
          //if (l1htbit) passMonitor=true; // No selection on the Monitoring seeds, selection stored
      }
  }
  //if (!passMonitor) {/*cout<<"failed L1 seed"<<endl;*/ return;} // No selection on the Monitoring seeds
  //if (!passScouting) {/*cout<<"failed L1 seed"<<endl;*/ return;} // No selection on the Scouting seeds
    
  Handle<double> rhoH;
  iEvent.getByToken(rhoToken, rhoH);


  // -------------- //
  // SCOUTING MUONS //
  // -------------- //

  Handle<vector<Run3ScoutingMuon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);
  nScoutingMuons=0;
  vector<int> idx_scout;
  int s=0;

  // Minimal selection on the scouting muons
  // ---------------------------------------
  for (auto muons_scout_iter = muonsH->begin(); muons_scout_iter != muonsH->end(); ++muons_scout_iter) {
    //cout<<"scouting muon pt: "<<muons_scout_iter->pt()<<" eta: "<<muons_scout_iter->eta()<<endl;                                 
    if (muons_scout_iter->pt()>3 and abs(muons_scout_iter->eta())<1.9) {
      nScoutingMuons+=1;
      idx_scout.push_back(s);
    }
    s+=1;
  }

  // Requiring exactly two scouting muons in the event
  // -------------------------------------------------
  if (!(idx_scout.size()==2)) {/*cout<<"failed scouting muons"<<endl;*/ return;}

  // Checking opposite charge                                                                                     
  // ------------------------
  int checkOSCharge_scout = muonsH->at(idx_scout[0]).charge() * muonsH->at(idx_scout[1]).charge();
  if (checkOSCharge_scout > 0){/*cout<<"not OS pair"<<endl;*/return;}


  // ---------------- //
  // VERTEX SELECTION //
  // ---------------- //

  // Checking on primary vertices to take the average
  // ------------------------------------------------
  Handle<vector<Run3ScoutingVertex> > primaryVerticesH;
  iEvent.getByToken(primaryVerticesToken, primaryVerticesH);

  vtxX.clear();
  vtxY.clear();
  vtxXError.clear();
  vtxYError.clear();

  nvtx = 0;
  for (auto vtx_iter = primaryVerticesH->begin(); vtx_iter != primaryVerticesH->end(); ++vtx_iter) {
    //std::cout<<"primary x: "<<vtx_iter->x() <<" y: "<<  vtx_iter->y()<<" ex: "<<vtx_iter->xError()<<" ey: "<<vtx_iter->yError()<<std::endl;
    nvtx++;
    vtxX.push_back(vtx_iter->x());   
    vtxY.push_back(vtx_iter->y());        
    vtxXError.push_back(vtx_iter->xError());        
    vtxYError.push_back(vtx_iter->yError());           
  }

  avgPrimaryX = 0;
  if (!vtxX.empty()) {
    avgPrimaryX = std::reduce(vtxX.begin(), vtxX.end(), 0.0) / vtxX.size();
  }

  avgPrimaryY = 0;
  if (!vtxY.empty()) {
    avgPrimaryY = std::reduce(vtxY.begin(), vtxY.end(), 0.0) / vtxY.size();
  }

  //std::cout<<"avgPrimaryX: "<<avgPrimaryX<<" avgPrimaryY: "<<avgPrimaryY<<std::endl;

  // Checking on displaced vertices
  // ------------------------------
  Handle<vector<Run3ScoutingVertex> > verticesH;
  iEvent.getByToken(verticesToken, verticesH);

  lxy=0;
  vtxErr=0;
  IP=0;

  std::vector<int> vtxIndx1 = (muonsH->at(idx_scout[0])).vtxIndx();
  std::vector<int> vtxIndx2 = (muonsH->at(idx_scout[1])).vtxIndx();

  //std::cout<<"vtxIndx1 size: "<<vtxIndx1.size()<<" vtxIndx2 size: "<<vtxIndx2.size()<<" num vtx: "<<verticesH->size()<<std::endl;
      
  ndvtx = verticesH->size();
  vtxMatch = vtxIndx1.size() > 0 && vtxIndx2.size() > 0 && vtxIndx1[0] == 0 && vtxIndx2[0] == 0 && ndvtx > 0;
      
  if (vtxMatch) {
    //std::cout<<"vtxIndx12: "<<(muonsH->at(idx_scout[0])).vtxIndx()[0]<<" "<<(muonsH->at(idx_scout[1])).vtxIndx()[0]<<endl;

    auto vtx = verticesH->begin();
    double dx = (vtx->x()) - avgPrimaryX;
    double dy = (vtx->y()) - avgPrimaryY;

    //std::cout<<"x: "<<vtx->x()<<" y: "<<vtx->y()<<std::endl;

    lxy = sqrt(dx*dx + dy*dy);
    vtxErr = sqrt(dx*dx*(vtx->xError())*(vtx->xError()) + dy*dy*(vtx->yError())*(vtx->yError())) / lxy;
    IP = lxy/vtxErr;
  }
  else{
    lxy = -999.;
    vtxErr = -999.;
    IP = -999.;
  }

  //std::cout<<"lxy: "<<lxy<<" vtxErr: "<<vtxErr<<" IP: "<<IP<<std::endl;


  //cout << "Filling offline muon quantities!" << endl;
  pt1=offlineMuonsH->at(idx[0]).pt();
  pt2=offlineMuonsH->at(idx[1]).pt();
  eta1=offlineMuonsH->at(idx[0]).eta();
  eta2=offlineMuonsH->at(idx[1]).eta();
  phi1=offlineMuonsH->at(idx[0]).phi();
  phi2=offlineMuonsH->at(idx[1]).phi();
  id1=offlineMuonsH->at(idx[0]).pdgId();
  id2=offlineMuonsH->at(idx[1]).pdgId();
  
  // Compute the pfIso for the muon. Note: PUCorr = 0.5*muons_iter->puChargedHadronIso()
  // -----------------------------------------------------------------------------------
  pfIso1 = (offlineMuonsH->at(idx[0]).chargedHadronIso() + std::max(offlineMuonsH->at(idx[0]).photonIso() + offlineMuonsH->at(idx[0]).neutralHadronIso() - 0.5*offlineMuonsH->at(idx[0]).puChargedHadronIso(), 0.0))/offlineMuonsH->at(idx[0]).pt();
  pfIso2 = (offlineMuonsH->at(idx[1]).chargedHadronIso() + std::max(offlineMuonsH->at(idx[1]).photonIso() + offlineMuonsH->at(idx[1]).neutralHadronIso() - 0.5*offlineMuonsH->at(idx[1]).puChargedHadronIso(), 0.0))/offlineMuonsH->at(idx[1]).pt();
  //std::cout << "----- pfIso: " << pfIso << endl;
  
  // Saving Muon ID info for various WP:
  // 0 isLooseMuon, 1 isMediumMuon, 2 isTightMuon
  // --------------------------------------------
  mu1_ID.push_back(offlineMuonsH->at(idx[0]).isLooseMuon());
  mu1_ID.push_back(offlineMuonsH->at(idx[0]).isMediumMuon());
  mu2_ID.push_back(offlineMuonsH->at(idx[1]).isLooseMuon());
  mu2_ID.push_back(offlineMuonsH->at(idx[1]).isMediumMuon());  
  
  TLorentzVector mu1;
  mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.105658);
  TLorentzVector mu2;
  mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.105658);
  TLorentzVector dimu = mu1+mu2;
  mass=dimu.M();
  ptmm=dimu.Pt();
  drmm=mu1.DeltaR(mu2);
  
  rho=*rhoH;

  //cout << "                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                      " <<endl;
  //cout << "%%% pt1_scout = " << muonsH->at(idx_scout[0]).pt() << endl;
  //cout << "%%% pt2_scout = " << muonsH->at(idx_scout[1]).pt() << endl;

    
  // Checking the muon matching based on dR                                                  
  // --------------------------------------
  float dr1 = 999.;
  int idx1_scout = -1;
  dr_matching_1 = 999.;

  float dr2 = 999.;
  int idx2_scout = -1;
  dr_matching_2 = 999.; 

  /*
  cout<<"############ BEFORE #############"<<endl;
  cout<<"dr_matching1 = "<<dr_matching_1<<"  and dr_matching2 = "<<dr_matching_2<<endl;
  cout<<"#################################"<<endl;
  */
  for (int mu_dr1=0; mu_dr1 < 2; mu_dr1++){
    /*
    cout << "---------- Loop on the scouting muons for mu0: checking mu " << mu_dr1 << endl;
    cout << "scoutMu pt: " << muonsH->at(idx_scout[mu_dr1]).pt()<<endl;
    cout << "scoutMu eta: " << muonsH->at(idx_scout[mu_dr1]).eta()<<endl;
    cout << "scoutMu phi: " << muonsH->at(idx_scout[mu_dr1]).phi()<<endl;
    */
    dr1 = sqrt( (offlineMuonsH->at(idx[0]).eta() - muonsH->at(idx_scout[mu_dr1]).eta())*(offlineMuonsH->at(idx[0]).eta() - muonsH->at(idx_scout[mu_dr1]).eta()) 
		   + (offlineMuonsH->at(idx[0]).phi() - muonsH->at(idx_scout[mu_dr1]).phi())*(offlineMuonsH->at(idx[0]).phi() - muonsH->at(idx_scout[mu_dr1]).phi()) );
    //cout<<"dR = "<<dr1<<endl;
    if (dr1 < dr_matching_1) {
      dr_matching_1 = dr1;
      idx1_scout = mu_dr1;
    }
    //cout<<" ** AFTER[1] dR = "<<dr_matching_1<<" and idx1_scout = "<< idx1_scout << endl;
  }
  
  for (int mu_dr2=0; mu_dr2 < 2; mu_dr2++){
    /*
    cout << "---------- Loop on the scouting muons for mu1: checking mu " << mu_dr2 << endl;
    cout << "scoutMu pt: " << muonsH->at(idx_scout[mu_dr2]).pt()<<endl;
    cout << "scoutMu eta: " << muonsH->at(idx_scout[mu_dr2]).eta()<<endl;
    cout << "scoutMu phi: " << muonsH->at(idx_scout[mu_dr2]).phi()<<endl;
    */
    dr2 = sqrt( (offlineMuonsH->at(idx[1]).eta() - muonsH->at(idx_scout[mu_dr2]).eta())*(offlineMuonsH->at(idx[1]).eta() - muonsH->at(idx_scout[mu_dr2]).eta()) 
	       + (offlineMuonsH->at(idx[1]).phi() - muonsH->at(idx_scout[mu_dr2]).phi())*(offlineMuonsH->at(idx[1]).phi() - muonsH->at(idx_scout[mu_dr2]).phi()) );
    //cout<<"dR = "<<dr2<<endl;
    if (dr2 < dr_matching_2) {
      dr_matching_2 = dr2;
      idx2_scout = mu_dr2;
    }
    //cout<<" ** AFTER[2] dR = "<<dr_matching_2<<" and idx2_scout = "<< idx2_scout << endl;
  }
  /*
  cout<<"############ AFTER #############"<<endl;
  cout<<"dr_matching1 = "<<dr_matching_1<<"  and dr_matching2 = "<<dr_matching_2<<endl;
  cout<<"idx1_scout = "<<idx1_scout<<"  and idx2_scout = "<<idx2_scout<<endl;
  cout<<"#################################"<<endl;
   
  cout<<">>>>>>>>>> scouting muon pt1: "<< muonsH->at(idx_scout[idx1_scout]).pt() << endl;                                 
  cout<<">>>>>>>>>> scouting muon pt2: "<< muonsH->at(idx_scout[idx2_scout]).pt() << endl;                                 

  cout<<">>>>>>>>>> offline muon pt1: "<< offlineMuonsH->at(idx[0]).pt() << endl;                                 
  cout<<">>>>>>>>>> offline muon pt2: "<< offlineMuonsH->at(idx[1]).pt() << endl;                                 
  */

  //if(dr_matching_1 > 0.2 || dr_matching_2 > 0.2){/*cout<<"no matching!"<<endl;*/return;};
  //cout << "-------------------------------PASSED dR Matching-------------------------"<<endl;
   
  if (idx1_scout == idx2_scout){
    if (dr_matching_1 < dr_matching_2) {idx1_scout = 0; idx2_scout = 1;}
    else if (dr_matching_1 > dr_matching_2) {idx1_scout = 1; idx2_scout = 0;}
  }

  // Filling Scouting quantities 
  // ---------------------------
  TLorentzVector mu1_scout;
  TLorentzVector mu2_scout;
  pt1_scout=muonsH->at(idx1_scout).pt();
  pt2_scout=muonsH->at(idx2_scout).pt();
  mu1_scout.SetPtEtaPhiM(pt1_scout,muonsH->at(idx1_scout).eta(),muonsH->at(idx1_scout).phi(),0.105658);
  mu2_scout.SetPtEtaPhiM(pt2_scout,muonsH->at(idx2_scout).eta(),muonsH->at(idx2_scout).phi(),0.105658);

  TLorentzVector dimu_scout = mu1_scout+mu2_scout;
  mass_scout=dimu_scout.M();
  ptmm_scout=dimu_scout.Pt();
  drmm_scout=mu1_scout.DeltaR(mu2_scout);

  /*cout<<"############ Checking the final angular separation #############"<<endl;
  cout<<"DIMUON dR = " << drmm << endl;
  cout<<"Scouting DIMUON dR = " << drmm_scout << endl;
  cout<<"################################################################"<<endl;
  */
  tree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void ScoutingTreeMakerRun3Monitor::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"      , "tree");
    tree->Branch("mass"                , &mass                         , "mass/F");
    tree->Branch("ptmm"                , &ptmm                         , "ptmm/F");
    tree->Branch("drmm"                , &drmm                         , "drmm/F");
    tree->Branch("pt1"                 , &pt1                          , "pt1/F");
    tree->Branch("pt2"                 , &pt2                          , "pt2/F");
    tree->Branch("eta1"                , &eta1                         , "eta1/F");
    tree->Branch("eta2"                , &eta2                         , "eta2/F");
    tree->Branch("id1"                 , &id1                          , "id1/I");
    tree->Branch("id2"                 , &id2                          , "id2/I");
    tree->Branch("rho"                 , &rho                          , "rho/F");
    tree->Branch("pfIso1"              , &pfIso1                       , "pfIso1/F");
    tree->Branch("pfIso2"              , &pfIso2                       , "pfIso2/F");
    tree->Branch("mu1_ID", "std::vector<bool>"                         , &mu1_ID, 32000, 0);
    tree->Branch("mu2_ID", "std::vector<bool>"                         , &mu2_ID, 32000, 0);
    tree->Branch("l1Result", "std::vector<bool>"                       , &l1Result_, 32000, 0);
    tree->Branch("l1Result_mon", "std::vector<bool>"                   , &l1Result_mon_, 32000, 0);
    tree->Branch("nOfflineMuons"       , &nOfflineMuons                , "nOfflineMuons/I");
    tree->Branch("nScoutingMuons"      , &nScoutingMuons               , "nScoutingMuons/I");
    tree->Branch("nvtx"                , &nvtx                         , "nvtx/I");
    tree->Branch("ndvtx"                , &ndvtx                         , "ndvtx/I");
    tree->Branch("pt1_scout"           , &pt1_scout                    , "pt1_scout/F");                                    
    tree->Branch("pt2_scout"           , &pt2_scout                    , "pt2_scout/F");                                      
    tree->Branch("mass_scout"          , &mass_scout                   , "mass_scout/F");                                          
    tree->Branch("ptmm_scout"          , &ptmm_scout                   , "ptmm_scout/F");                           
    tree->Branch("drmm_scout"          , &drmm_scout                   , "drmm_scout/F");
    tree->Branch("lxy"                 , &lxy                          , "lxy/F"     );
    tree->Branch("vtxErr"              , &vtxErr                       , "vtxErr/F"  );
    tree->Branch("vtxMatch"            , &vtxMatch                     , "vtxMatch/B");
    tree->Branch("IP"                  , &IP                           , "IP/F"      );
    tree->Branch("dr_matching_1"       , &dr_matching_1                , "dr_matching_1/F");
    tree->Branch("dr_matching_2"       , &dr_matching_2                , "dr_matching_2/F");
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
