// -*- C++ -*-
//
// Package:    Run3ScoutingAnalysisTools/ScoutingTreeMakerRun3
// Class:      ScoutingTreeMakerRun3
//
/**\class ScoutingTreeMakerRun3 ScoutingTreeMakerRun3.cc Run3ScoutingAnalysisTools/ScoutingTreeMakerRun3/plugins/ScoutingTreeMakerRun3.cc

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

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class ScoutingTreeMakerRun3 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingTreeMakerRun3(const edm::ParameterSet&);
  ~ScoutingTreeMakerRun3() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
    //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;   

  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>             triggerResultsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon> >      muonsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> >  electronsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >    primaryVerticesToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >    verticesToken;
  const edm::EDGetTokenT<double>                          rhoToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> >  photonsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle> >  pfcandsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPFJet> >  pfjetsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingTrack> >  tracksToken;

  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;

  bool doL1;
  triggerExpression::Data triggerCache_;

  edm::InputTag                algInputTag_;
  edm::InputTag                extInputTag_;
  edm::EDGetToken              algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;
  std::vector<bool>            l1Result_;

  TTree* tree;

  float trackIso1;
  float trackIso2;
  int nValidPixelHits1;
  int nValidPixelHits2;
  int nTrackerLayersWithMeasurement1;
  int nTrackerLayersWithMeasurement2;
  float trk_chi21;
  float trk_chi22;

  bool muonID1;
  bool muonID2;
  float mass;
  float pt;
  float dr;
  float pt1;
  float pt2;
  float eta1;
  float eta2;

  float rho;
  int nMuonsID;

  bool hasPvtx;

  int ndvtx;
  bool isValidVtx;
  float vtxChi2;
  int vtxNdof;
  bool vtxMatch;

  float vtxXError;
  float vtxYError;
  float vtxZError;

  float Lxy;
  float LxyErr;
  float LxySig;
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
ScoutingTreeMakerRun3::ScoutingTreeMakerRun3(const edm::ParameterSet& iConfig):
    triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
    triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),
    muonsToken               (consumes<std::vector<Run3ScoutingMuon> >             (iConfig.getParameter<edm::InputTag>("muons"))),
    electronsToken           (consumes<std::vector<Run3ScoutingElectron> >         (iConfig.getParameter<edm::InputTag>("electrons"))),
    primaryVerticesToken     (consumes<std::vector<Run3ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("primaryVertices"))),
    verticesToken            (consumes<std::vector<Run3ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("displacedVertices"))),
    rhoToken                 (consumes<double>                                 (iConfig.getParameter<edm::InputTag>("rho"))), 
    photonsToken             (consumes<std::vector<Run3ScoutingPhoton> >         (iConfig.getParameter<edm::InputTag>("photons"))),
    pfcandsToken             (consumes<std::vector<Run3ScoutingParticle> >         (iConfig.getParameter<edm::InputTag>("pfcands"))),
    pfjetsToken              (consumes<std::vector<Run3ScoutingPFJet> >            (iConfig.getParameter<edm::InputTag>("pfjets"))),
    tracksToken              (consumes<std::vector<Run3ScoutingTrack> >            (iConfig.getParameter<edm::InputTag>("tracks"))),
    doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1")            : false)
{
    usesResource("TFileService");
    if (doL1) {
        algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
        extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
        algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
        l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
        l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
    }
    else {
        l1Seeds_ = std::vector<std::string>();
        l1GtUtils_ = 0;
    }
}

ScoutingTreeMakerRun3::~ScoutingTreeMakerRun3() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void ScoutingTreeMakerRun3::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;

  Handle<vector<Run3ScoutingMuon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);

  if (muonsH->size()<2) return;

  int nMuons=0;
  nMuonsID=0;
  vector<int> idx;

  int j=0;
  for (auto muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
      //std::cout<<"pt: "<<muons_iter->pt()<<std::endl;                                             
      //std::cout<<"trkiso: "<<muons_iter->trackIso()<<" pix hits: "<< muons_iter->nValidPixelHits()<<" layers: "<<muons_iter->nTrackerLayersWithMeasurement()<<" trk chi2: "<<muons_iter->trk_chi2()<<std::endl; 

      /*
      if (muons_iter->pt()>4) {
          nMuons+=1;
          if ((muons_iter->trackIso()<0.15) &&
              (muons_iter->nValidPixelHits()>0) &&
              (muons_iter->nTrackerLayersWithMeasurement()>5)&&
              (muons_iter->trk_chi2()<10)) {
              nMuonsID+=1;
          }
      }
      */

      idx.push_back(j);
      j+=1;
  }

  //std::cout<<std::endl<<idx.size()<<std::endl;

  if (idx.size()>1) {
      //std::cout << "charge: " << (muonsH->at(idx[0]).charge()) << ", " << (muonsH->at(idx[1]).charge()) << std::endl;
      if ((muonsH->at(idx[0]).charge())*(muonsH->at(idx[1]).charge()) > 0) {
        return;
      }

      //muonID1 = (muonsH->at(idx[0]).pt()>4) && (muonsH->at(idx[0]).trackIso()<0.15) && (muonsH->at(idx[0]).nValidPixelHits()>0) && (muonsH->at(idx[0]).trk_chi2()<10);
      //muonID2 = (muonsH->at(idx[1]).pt()>4) && (muonsH->at(idx[1]).trackIso()<0.15) && (muonsH->at(idx[1]).nValidPixelHits()>0) && (muonsH->at(idx[1]).trk_chi2()<10);

      //std::cout << "ID" << muonID1 << ", " << muonID2 << std::endl;

      trackIso1 = muonsH->at(idx[0]).trackIso();
      trackIso2 = muonsH->at(idx[1]).trackIso();   
      nValidPixelHits1 = muonsH->at(idx[0]).nValidPixelHits();
      nValidPixelHits2 = muonsH->at(idx[1]).nValidPixelHits();
      nTrackerLayersWithMeasurement1 = muonsH->at(idx[0]).nTrackerLayersWithMeasurement();
      nTrackerLayersWithMeasurement2 = muonsH->at(idx[1]).nTrackerLayersWithMeasurement();
      trk_chi21 = muonsH->at(idx[0]).trk_chi2();
      trk_chi22 = muonsH->at(idx[1]).trk_chi2();

      pt1=muonsH->at(idx[0]).pt();
      pt2=muonsH->at(idx[1]).pt();

      eta1=muonsH->at(idx[0]).eta();
      eta2=muonsH->at(idx[1]).eta();      
      float phi1=muonsH->at(idx[0]).phi();
      float phi2=muonsH->at(idx[1]).phi();
      
      TLorentzVector mu1;
      mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.105658);

      TLorentzVector mu2;
      mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.105658);

      TLorentzVector dimu = mu1+mu2;
      mass=dimu.M();
      pt=dimu.Pt();
      dr=mu1.DeltaR(mu2);

      //std::cout<<"pt: "<<pt1<<", "<<pt2<<", nMuonsID: "<<nMuonsID<<std::endl;

      Handle<double> rhoH;
      iEvent.getByToken(rhoToken, rhoH);
      rho=*rhoH;

      Handle<vector<Run3ScoutingVertex> > primaryVerticesH;
      iEvent.getByToken(primaryVerticesToken, primaryVerticesH);

      std::vector<float> vtxX;
      std::vector<float> vtxY;

      int npvtx = 0;
      for (auto vtx_iter = primaryVerticesH->begin(); vtx_iter != primaryVerticesH->end(); ++vtx_iter) {
        //std::cout<<"primary x: "<<vtx_iter->x() <<" y: "<<  vtx_iter->y()<<" ex: "<<vtx_iter->xError()<<" ey: "<<vtx_iter->yError()<<std::endl;
        npvtx++;
        vtxX.push_back(vtx_iter->x());   
        vtxY.push_back(vtx_iter->y());                
      }

      hasPvtx = npvtx > 0;

      float avgPrimary[2];
      avgPrimary[0] = ( vtxX.empty() ) ? 0 : ( std::reduce(vtxX.begin(), vtxX.end(), 0.0) / vtxX.size() );
      avgPrimary[1] = ( vtxY.empty() ) ? 0 : ( std::reduce(vtxY.begin(), vtxY.end(), 0.0) / vtxY.size() );

      //std::cout << "npvtx: " << npvtx << " avgPrimaryX: " << avgPrimary[0] << " avgPrimaryY: " << avgPrimary[1] << std::endl;

      Handle<vector<Run3ScoutingVertex> > verticesH;
      iEvent.getByToken(verticesToken, verticesH);

      std::vector<int> vtxIndx1 = (muonsH->at(idx[0])).vtxIndx();
      std::vector<int> vtxIndx2 = (muonsH->at(idx[1])).vtxIndx();

      //std::cout<<"vtxIndx1 size: "<<vtxIndx1.size()<<" vtxIndx2 size: "<<vtxIndx2.size()<<" num vtx: "<<verticesH->size()<<std::endl;
      
      ndvtx = verticesH->size();
      vtxMatch = vtxIndx1.size() > 0 && vtxIndx2.size() > 0 && vtxIndx1[0] == 0 && vtxIndx2[0] == 0 && ndvtx > 0;

      Lxy=0;
      LxyErr=0;
      LxySig=0;

      if (vtxMatch) {
        auto vtx = verticesH->begin();
        double dx = (vtx->x()) - avgPrimary[0];
        double dy = (vtx->y()) - avgPrimary[1];

        //std::cout<<"x: "<<vtx->x()<<" y: "<<vtx->y()<<std::endl;

        vtxChi2 = vtx->chi2();
        vtxNdof = vtx->ndof();
        isValidVtx = vtx->isValidVtx();

        vtxXError = vtx->xError();
        vtxYError = vtx->yError();
        vtxZError = vtx->zError();

        Lxy = sqrt(dx*dx + dy*dy);
        LxyErr = sqrt(dx*dx*(vtx->xError())*(vtx->xError()) + dy*dy*(vtx->yError())*(vtx->yError())) / Lxy;
        LxySig = Lxy/LxyErr;
      }

      //std::cout<<"Lxy: "<<Lxy<<" LxyErr: "<<LxyErr<<" LxySig: "<<LxySig<<std::endl;
      
      l1Result_.clear();
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
          }
      }

      //std::cout<<"tree filling with mass: "<<mass<<", pt: "<<pt<<std::endl;
      tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void ScoutingTreeMakerRun3::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"      , "tree");
    //tree->Branch("muonID1"             , &muonID1                     , "muonID1/B" );
    //tree->Branch("muonID2"             , &muonID2                     , "muonID2/B" );

    tree->Branch("trackIso1", &trackIso1, "trackIso1/F");
    tree->Branch("trackIso2", &trackIso2, "trackIso2/F");
    tree->Branch("nValidPixelHits1", &nValidPixelHits1, "nValidPixelHits1/I");
    tree->Branch("nValidPixelHits2", &nValidPixelHits2, "nValidPixelHits2/I");
    tree->Branch("nTrackerLayersWithMeasurement1", &nTrackerLayersWithMeasurement1, "nTrackerLayersWithMeasurement1/I");
    tree->Branch("nTrackerLayersWithMeasurement2", &nTrackerLayersWithMeasurement2, "nTrackerLayersWithMeasurement2/I");
    tree->Branch("trk_chi21", &trk_chi21, "trk_chi21/F");
    tree->Branch("trk_chi22", &trk_chi22, "trk_chi22/F");

    tree->Branch("mass"                , &mass                        , "mass/F"    );
    tree->Branch("pt"                  , &pt                          , "pt/F"      );
    tree->Branch("dr"                  , &dr                          , "dr/F"      );
    tree->Branch("pt1"                 , &pt1                         , "pt1/F"     );
    tree->Branch("pt2"                 , &pt2                         , "pt2/F"     );
    tree->Branch("eta1"                , &eta1                        , "eta1/F"    );
    tree->Branch("eta2"                , &eta2                        , "eta2/F"    );
    tree->Branch("rho"                 , &rho                         , "rho/F"     );

    tree->Branch("vtxMatch"            , &vtxMatch                    , "vtxMatch/B");
    tree->Branch("vtxChi2"             , &vtxChi2                     , "vtxChi2/F" );
    tree->Branch("vtxNdof"             , &vtxNdof                     , "vtxNdof/I" );
    tree->Branch("Lxy"                 , &Lxy                         , "Lxy/F"     );
    tree->Branch("LxyErr"              , &LxyErr                      , "LxyErr/F"  );
    tree->Branch("LxySig"              , &LxySig                      , "LxySig/F"  );

    tree->Branch("vtxXError"           , &vtxXError                   , "vtxXError/F");
    tree->Branch("vtxYError"           , &vtxYError                   , "vtxYError/F");
    tree->Branch("vtxZError"           , &vtxZError                   , "vtxZError/F");

    tree->Branch("l1Result", "std::vector<bool>"             ,&l1Result_, 32000, 0  );
}

// ------------ method called once each job just after ending the event loop  ------------
void ScoutingTreeMakerRun3::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingTreeMakerRun3::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(ScoutingTreeMakerRun3);
