// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>

// ROOT includes
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// CMSSW data formats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// Other relevant CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include <DataFormats/TrackReco/interface/TrackBase.h>

#include "DataFormats/Math/interface/libminifloat.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/TrackReco/interface/fillCovariance.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"

// Root include files
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// User include files

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/XConePlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/contrib/SoftKiller.hh"

using namespace std;


class ScoutingNanoAOD : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit ScoutingNanoAOD(const edm::ParameterSet&);
  ~ScoutingNanoAOD();
		
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void clearVars();
  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>             	triggerResultsToken;

  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon> >      muonsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> >  	electronsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> >  	photonsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle> >  	pfcandsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPFJet> >  	pfjetsToken;
 const edm::EDGetTokenT<std::vector<Run3ScoutingTrack> >  	tracksToken;
  

  //const edm::EDGetTokenT<GenEventInfoProduct>             genEvtInfoToken;

  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;

        
	
  bool doL1;       
  triggerExpression::Data triggerCache_;
      
	
  // Generator-level information
  // Flags for the different types of triggers used in the analysis
  // For now we are interested in events passing either the single or double lepton triggers
  unsigned char                trig;
       
  edm::InputTag                algInputTag_;       
  edm::InputTag                extInputTag_;       
  edm::EDGetToken              algToken_;
  //l1t::L1TGlobalUtil          *l1GtUtils_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;
  std::vector<bool>            l1Result_;
       
        


  //Photon
  const static int 	max_pho = 1000;
  UInt_t n_pho;
  vector<Float16_t> 	    	Photon_pt;
  vector<Float16_t>        	Photon_eta;
  vector<Float16_t>        	Photon_phi;
  vector<Float16_t>	    	Photon_m;
  vector<Float16_t>	    	Photon_sigmaietaieta;
  vector<Float16_t>	    	Photon_HoE;
  vector<Float16_t>        	Photon_ecaliso;
  vector<Float16_t>	    	Photon_hcaliso;

  //Electron
  const static int 	max_ele = 1000;
  UInt_t n_ele;
  vector<Float16_t> 	    Electron_pt;
  vector<Float16_t>        Electron_eta;
  vector<Float16_t>        Electron_phi;
  vector<Float16_t>	    Electron_m;
  vector<Float16_t>        Electron_d0;
  vector<Float16_t>	    Electron_dz;
  vector<Float16_t>	    Electron_detain;
  vector<Float16_t>	    Electron_dphiin;
  vector<Float16_t>	    Electron_sigmaietaieta;
  vector<Float16_t>	    Electron_HoE;
  vector<Float16_t>	    Electron_ooEMOop;
  vector<Float16_t>	    Electron_mHits;
  vector<Float16_t>        Electron_charge;
  vector<Float16_t>        Electron_ecaliso;
  vector<Float16_t>	    Electron_hcaliso;
  vector<Float16_t>        Electron_tkiso;

  //Muon
  const static int 	max_mu = 1000;
  UInt_t n_mu;
  vector<Float16_t> 	Muon_pt;
  vector<Float16_t> 	Muon_eta;
  vector<Float16_t> 	Muon_phi;
  vector<Float16_t> 	Muon_m;
  vector<Float16_t> 	Muon_ecaliso;
  vector<Float16_t> 	Muon_hcaliso;
  vector<Float16_t> 	Muon_trkiso;
  vector<Float16_t> 	Muon_chi2;
  vector<Float16_t> 	Muon_ndof;
  vector<Float16_t> 	Muon_charge;
  vector<Float16_t> 	Muon_dxy;
  vector<Float16_t> 	Muon_dz;
  vector<Float16_t> 	Muon_nvalidmuon_hits;
  vector<Float16_t> 	Muon_nvalidpixelhits;
  
  vector<Float16_t> 	Muon_nmatchedstations;
  vector<Float16_t>       Muon_type;
  vector<Float16_t>       Muon_nvalidstriphits;
  vector<Float16_t>       Muon_trkqoverp;
  vector<Float16_t>       Muon_trklambda;
  vector<Float16_t>       Muon_trkpt;
  vector<Float16_t>       Muon_trkphi;
  vector<Float16_t>       Muon_trketa;
  vector<Float16_t>       Muon_trkqoverperror;
  vector<Float16_t>       Muon_trklambdaerror;
  vector<Float16_t>       Muon_trkpterror;
  vector<Float16_t>       Muon_trkphierror;
  vector<Float16_t>       Muon_trketaerror;
  vector<Float16_t>       Muon_trkdszerror;
  vector<Float16_t>       Muon_trkdsz;
  int muontvtxind[max_mu];

  //PFJets
  const static int 	max_jet = 1000;
  UInt_t n_jet;
  vector<Float16_t> 	    Jet_pt;
  vector<Float16_t>         Jet_eta;
  vector<Float16_t>         Jet_phi;
  vector<Float16_t>	    Jet_m;
  vector<Float16_t>	    Jet_area;
  vector<Float16_t>	    Jet_chargedHadronEnergy;
  vector<Float16_t>         Jet_neutralHadronEnergy;
  vector<Float16_t>	    Jet_photonEnergy;
  vector<Float16_t>	    Jet_electronEnergy;
  vector<Float16_t>	    Jet_muonEnergy;
  vector<Float16_t>	    Jet_HFHadronEnergy;
  vector<Float16_t>	    Jet_HFEMEnergy;
  vector<Float16_t>	    Jet_HOEnergy;
  vector<Float16_t>	    Jet_chargedHadronMultiplicity;
  vector<Float16_t>         Jet_neutralHadronMultiplicity;
  vector<Float16_t>	    Jet_photonMultiplicity;
  vector<Float16_t>	    Jet_electronMultiplicity;
  vector<Float16_t>	    Jet_muonMultiplicity;
  vector<Float16_t>	    Jet_HFHadronMultiplicity;
  vector<Float16_t>	    Jet_HFEMMultiplicity;
  vector<Float16_t> 	    Jet_csv;
  vector<Float16_t> 	    Jet_mvaDiscriminator;
  std::vector< std::vector<int16_t> >  	    Jet_constituents;

  //PFCand
  const static int 	max_pfcand = 10000;
  UInt_t n_pfcand;
  vector<Float16_t> 	    pfcandpt;
  vector<Float16_t>         pfcandeta;
  vector<Float16_t>         pfcandphi;
  vector<Float16_t>	    pdcandm;
  vector<Float16_t>	    pfcandpdgid;
  vector<Float16_t>	    pfcandvertex;

  UInt_t n_fatjet;
  vector<Float16_t> FatJet_area;
  vector<Float16_t> FatJet_eta;
  vector<Float16_t> FatJet_n2b1;
  vector<Float16_t> FatJet_n3b1;
  vector<Float16_t> FatJet_phi;
  vector<Float16_t> FatJet_pt;
  vector<Float16_t> FatJet_tau1;
  vector<Float16_t> FatJet_tau2;
  vector<Float16_t> FatJet_tau3;
  vector<Float16_t> FatJet_tau4;
  vector<Float16_t> FatJet_mass;
  vector<Float16_t> FatJet_msoftdrop;
  vector<Float16_t> FatJet_mtrim;
        
  // TTree carrying the event weight information
  TTree* tree;

  //Run and lumisection
  int run;
  int lumSec;

};

ScoutingNanoAOD::ScoutingNanoAOD(const edm::ParameterSet& iConfig): 
  triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
  triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),


  muonsToken               (consumes<std::vector<Run3ScoutingMuon> >             (iConfig.getParameter<edm::InputTag>("muons"))), 
  electronsToken           (consumes<std::vector<Run3ScoutingElectron> >         (iConfig.getParameter<edm::InputTag>("electrons"))), 
  photonsToken           (consumes<std::vector<Run3ScoutingPhoton> >         (iConfig.getParameter<edm::InputTag>("photons"))), 
  pfcandsToken             (consumes<std::vector<Run3ScoutingParticle> >         (iConfig.getParameter<edm::InputTag>("pfcands"))), 
  pfjetsToken              (consumes<std::vector<Run3ScoutingPFJet> >            (iConfig.getParameter<edm::InputTag>("pfjets"))),
 tracksToken              (consumes<std::vector<Run3ScoutingTrack> >            (iConfig.getParameter<edm::InputTag>("tracks"))), 
//  pileupInfoToken          (consumes<std::vector<PileupSummaryInfo> >        (iConfig.getParameter<edm::InputTag>("pileupinfo"))),
//  gensToken                (consumes<std::vector<reco::GenParticle> >        (iConfig.getParameter<edm::InputTag>("gens"))),
  //genEvtInfoToken          (consumes<GenEventInfoProduct>                    (iConfig.getParameter<edm::InputTag>("geneventinfo"))),    
  doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1")            : false)
{
  usesResource("TFileService");
  if (doL1) {
   algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
   extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
   algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
   l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
    /* l1GtUtils_ = new l1t::L1TGlobalUtil(iConfig,consumesCollector());*/	
   l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(
    iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
  }
  else {
    l1Seeds_ = std::vector<std::string>();
    l1GtUtils_ = 0;
  }

 // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree"       , "tree");

  // Event weights
    
  tree->Branch("lumSec"		, &lumSec			 , "lumSec/i" );
  tree->Branch("run"			, &run				 , "run/i" );
  //tree->Branch("nvtx"			, &nvtx				 , "nvtx/i" );
    
  // Triggers
  tree->Branch("trig"                 , &trig                          , "trig/b");
  tree->Branch("l1Result"		, "std::vector<bool>"             ,&l1Result_	, 32000, 0);		
  // Pileup info
  //tree->Branch("nvtx"                 , &nvtx                          , "nvtx/i"       );

  //Electrons
  tree->Branch("n_ele"            	   ,&n_ele 			, "n_ele/i"		);
  tree->Branch("Electron_pt"         ,&Electron_pt 		 		);
  tree->Branch("Electron_eta"               ,&Electron_eta 		  	);
  tree->Branch("Electron_phi"               ,&Electron_phi 		 	);
  tree->Branch("Electron_charge"            ,&Electron_charge 		 	);
  tree->Branch("Electron_m"            	   ,&Electron_m 			 );
tree->Branch("Electron_tkiso"               ,&Electron_tkiso 		 );
tree->Branch("Electron_HoE"            	   ,&Electron_HoE 		 );
tree->Branch("Electron_sigmaietaieta"       ,&Electron_sigmaietaieta 	 );
 tree->Branch("Electron_dphiin"              ,&Electron_dphiin 		 );
 tree->Branch("Electron_detain"              ,&Electron_detain 		 );
 tree->Branch("Electron_mHits"               ,&Electron_mHits 		 );
 tree->Branch("Electron_ooEMOop"             ,&Electron_ooEMOop  		 );

  //Photons
  tree->Branch("n_pho"            	   ,&n_pho 			, "n_pho/i"		);
  tree->Branch("Photon_pt"            	   ,&Photon_pt 			);
  tree->Branch("Photon_eta"            	   ,&Photon_eta 			);
  tree->Branch("Photon_phi"            	   ,&Photon_phi 			);	
  tree->Branch("Photon_m"            	   ,&Photon_m 			);
  tree->Branch("Photon_hcaliso"             ,&Photon_hcaliso 		);
  tree->Branch("Photon_ecaliso"             ,&Photon_ecaliso 		);
  tree->Branch("Photon_HoE"            	   ,&Photon_HoE 			);
  tree->Branch("Photon_sigmaietaieta"       ,&Photon_sigmaietaieta		 );

  tree->Branch("n_pfcand"            	   ,&n_pfcand 		,"n_pfcand/i"		);	
  tree->Branch("pfcandpt"        	   ,&pfcandpt 		 );
  tree->Branch("pfcandeta"            	   ,&pfcandeta 		 );
  tree->Branch("pfcandphi"            	   ,&pfcandphi		 );
  tree->Branch("pdcandm"            	   ,&pdcandm 		 );
  tree->Branch("pfcandpdgid"               ,&pfcandpdgid		 );
  tree->Branch("pfcandvertex"              ,&pfcandvertex 	 );

  tree->Branch("n_mu"            	   ,&n_mu 			, "n_mu/i"		);
  tree->Branch("Muon_pt", &Muon_pt	);
  tree->Branch("Muon_eta", &Muon_eta	);
  tree->Branch("Muon_phi", &Muon_phi	);
  tree->Branch("Muon_m", &Muon_m	);
  tree->Branch("Muon_ecaliso", &Muon_ecaliso	);
  tree->Branch("Muon_hcaliso", &Muon_hcaliso	);
  tree->Branch("Muon_trkiso", &Muon_trkiso	);
  tree->Branch("Muon_chi2", &Muon_chi2	);
  tree->Branch("Muon_ndof", &Muon_ndof	);
  tree->Branch("Muon_charge", &Muon_charge	);
  tree->Branch("Muon_dxy", &Muon_dxy	);
  tree->Branch("Muon_dz", &Muon_dz	);
  tree->Branch("Muon_nvalidmuon_hits", &Muon_nvalidmuon_hits	);
  tree->Branch("Muon_validpixelhits", &Muon_nvalidpixelhits  );
  
  tree->Branch("Muon_nmatchedstations", &Muon_nmatchedstations	);
  tree->Branch("Muon_type",   &Muon_type    );
  tree->Branch("Muon_nvalidstriphits",    &Muon_nvalidstriphits   );
  tree->Branch("Muon_trkqoverp",    &Muon_trkqoverp   );
  tree->Branch("Muon_trklambda",   &Muon_trklambda    );
  tree->Branch("Muon_trkpt",   &Muon_trkpt    );
  tree->Branch("Muon_trkphi",  &Muon_trkphi     );
  tree->Branch("Muon_trketa",   &Muon_trketa    );
  tree->Branch("Muon_trkqoverperror",   &Muon_trkqoverperror    );
  tree->Branch("Muon_trklambdaerror",   &Muon_trklambdaerror    );
  tree->Branch("Muon_trkpterror",   &Muon_trkpterror    );
  tree->Branch("Muon_trkphierror",   &Muon_trkphierror    );
  tree->Branch("Muon_trketaerror",   &Muon_trketaerror    );
  tree->Branch("Muon_trkdszerror",   &Muon_trkdszerror    );
  tree->Branch("Muon_trkdsz",    &Muon_trkdsz   );


  tree->Branch("n_jet"            	   	,&n_jet 			, "n_jet/i"		);
  tree->Branch("Jet_pt"            	   	,&Jet_pt 				 );
  tree->Branch("Jet_eta"            	   	,&Jet_eta 			 );
  tree->Branch("Jet_phi"            	   	,&Jet_phi 			 );
  tree->Branch("Jet_m"            	   	,&Jet_m 				 );
  tree->Branch("Jet_area"            	   	,&Jet_area			 );
  tree->Branch("Jet_chargedHadronEnergy"         ,&Jet_chargedHadronEnergy 	 );
  tree->Branch("Jet_neutralHadronEnergy"         ,&Jet_neutralHadronEnergy 	 );
  tree->Branch("Jet_photonEnergy"            	,&Jet_photonEnergy 		 );
  tree->Branch("Jet_electronEnergy"              ,&Jet_electronEnergy 		 );
  tree->Branch("Jet_muonEnergy"    		   ,&Jet_muonEnergy 		 );
  tree->Branch("Jet_HFHadronEnergy"            	   ,&Jet_HFHadronEnergy 		 );
  tree->Branch("Jet_HFEMEnergy"            	   ,&Jet_HFEMEnergy 		 );
  tree->Branch("Jet_HOEnergy"            	   ,&Jet_HOEnergy 		 );
  tree->Branch("Jet_chargedHadronMultiplicity"      ,&Jet_chargedHadronMultiplicity 		 );
  tree->Branch("Jet_neutralHadronMultiplicity"      ,&Jet_neutralHadronMultiplicity 		 );
  tree->Branch("Jet_photonMultiplicity"            	   ,&Jet_photonMultiplicity 		 );
  tree->Branch("Jet_electronMultiplicity"            	   ,&Jet_electronMultiplicity 		 );
  tree->Branch("Jet_muonMultiplicity"            	   ,&Jet_muonMultiplicity 		 );
  tree->Branch("Jet_HFHadronMultiplicity"            	   ,&Jet_HFHadronMultiplicity 		 );
  tree->Branch("Jet_HFEMMultiplicity"            	   ,&Jet_HFEMMultiplicity 		 );
  tree->Branch("Jet_csv"            	   	,&Jet_csv 		 );
  tree->Branch("Jet_mvaDiscriminator"            	   ,&Jet_mvaDiscriminator 		 );
  tree->Branch("Jet_constituents"            	, "std::vector< vector<int16_t> >"   , &Jet_constituents 		, 32000, 0);
  
  tree->Branch("FatJet_area",&FatJet_area);
  tree->Branch("FatJet_eta",&FatJet_eta);
  tree->Branch("FatJet_n2b1",&FatJet_n2b1);
  tree->Branch("FatJet_n3b1",&FatJet_n3b1);
  tree->Branch("FatJet_phi",&FatJet_phi);
  tree->Branch("FatJet_pt",&FatJet_pt);
  tree->Branch("FatJet_tau1",&FatJet_tau1);
  tree->Branch("FatJet_tau2",&FatJet_tau2);
  tree->Branch("FatJet_tau3",&FatJet_tau3);
  tree->Branch("FatJet_tau4",&FatJet_tau4);
  tree->Branch("FatJet_mass",&FatJet_mass);
  tree->Branch("FatJet_msoftdrop",&FatJet_msoftdrop);
  tree->Branch("FatJet_mtrim",&FatJet_mtrim);
  

 
}


ScoutingNanoAOD::~ScoutingNanoAOD() {
}

void ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace fastjet;
  using namespace fastjet::contrib;
    
  // Handles to the EDM content
  edm::Handle<edm::TriggerResults> triggerResultsH;
  iEvent.getByToken(triggerResultsToken, triggerResultsH);
    
  Handle<vector<Run3ScoutingElectron> > electronsH;
  iEvent.getByToken(electronsToken, electronsH);

  Handle<vector<Run3ScoutingMuon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);

  Handle<vector<Run3ScoutingPhoton> > photonsH;
  iEvent.getByToken(photonsToken, photonsH);

  Handle<vector<Run3ScoutingPFJet> > pfjetsH;
  iEvent.getByToken(pfjetsToken, pfjetsH);
    
  Handle<vector<Run3ScoutingParticle> > pfcandsH;
  iEvent.getByToken(pfcandsToken, pfcandsH);

  Handle<vector<Run3ScoutingTrack> > tracksH;
  iEvent.getByToken(tracksToken, tracksH);

  run = iEvent.eventAuxiliary().run();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();


  // Which triggers fired
  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
    if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;
    if (i == 0  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   1; // DST_L1DoubleMu_CaloScouting_PFScouting
    if (i == 1  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   2; // DST_DoubleMu3_Mass10_CaloScouting_PFScouting
    if (i == 2  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   4; // DST_ZeroBias_CaloScouting_PFScouting
    if (i == 3  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=   8; // DST_L1HTT_CaloScouting_PFScouting
    if (i == 4  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=  16; // DST_CaloJet40_CaloScouting_PFScouting
    if (i == 5  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=  32; // DST_HT250_CaloScouting
    if (i == 6  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig +=  64; // DST_HT410_PFScouting
    if (i == 7  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) trig += 128; // DST_HT450_PFScouting
  }
  
  Jet_constituents.clear();
  //0.1396
  //built transient tracks, check Vertex fitting
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);

  vector<TransientTrack> t_tks;
  std::unique_ptr<VertexCollection> vertexCollection(new VertexCollection());

  /* for (auto tracks_iter = tracksH->begin(); tracks_iter != tracksH->end(); ++tracks_iter) 
    {
      Particle::PolarLorentzVector pp;
      Particle::LorentzVector p;
      pp=Particle::PolarLorentzVector(tracks_iter->tk_pt(), tracks_iter->tk_eta(), tracks_iter->tk_phi(), 0.1396);
      p =Particle::LorentzVector(pp);
      Track::Point  tk1_refPoint(tracks_iter->tk_vx(), tracks_iter->tk_vy(), tracks_iter->tk_vz());
      Track::Vector tk1_momentum(p.px(), p.py(), p.pz());

      float tk1_covdata[]={
      tracks_iter->tk_qoverp_Error() * tracks_iter->tk_qoverp_Error(),
      tracks_iter->tk_qoverp_lambda_cov(),
      tracks_iter->tk_lambda_Error() * tracks_iter->tk_lambda_Error(),
      tracks_iter->tk_qoverp_phi_cov(),
      tracks_iter->tk_lambda_phi_cov(),
      tracks_iter->tk_phi_Error() * tracks_iter->tk_phi_Error(),
      tracks_iter->tk_qoverp_dxy_cov(),
      tracks_iter->tk_lambda_dxy_cov(),
      tracks_iter->tk_phi_dxy_cov(),
      tracks_iter->tk_dxy_Error() * tracks_iter->tk_dxy_Error(),
      tracks_iter->tk_qoverp_dsz_cov(),
      tracks_iter->tk_lambda_dsz_cov(),
      tracks_iter->tk_phi_dsz_cov(),
      tracks_iter->tk_dxy_dsz_cov(),
      tracks_iter->tk_dsz_Error() * tracks_iter->tk_dsz_Error()
    };
    reco::TrackBase::CovarianceMatrix tk1_cov;
    fillCovariance(tk1_cov, tk1_covdata);

 	Track tk1(tracks_iter->tk_chi2(), 
	      tracks_iter->tk_ndof(), 
	      tk1_refPoint,
	      tk1_momentum,
	      tracks_iter->tk_charge(),
	      tk1_cov);

	TransientTrack ttkp1 = (*theB).build(tk1);
	t_tks.push_back(ttkp1);
    }

	KalmanVertexFitter kvf;

      TransientVertex tv = kvf.vertex(t_tks);


  std::cout<<"vertex position: "<<tv.position().x()<<"  "<<tv.position().y()<<std::endl;
  */
  n_ele = 0;
  for (auto electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) 
    {
      Electron_pt.push_back(electrons_iter->pt());
      Electron_eta.push_back(electrons_iter->eta());
      Electron_phi.push_back(electrons_iter->phi());	
      Electron_m.push_back(electrons_iter->m());
      Electron_detain.push_back(electrons_iter->dEtaIn());
      Electron_dphiin.push_back(electrons_iter->dPhiIn());
      Electron_sigmaietaieta.push_back(electrons_iter->sigmaIetaIeta());
      Electron_HoE.push_back(electrons_iter->hOverE());	
      Electron_ooEMOop.push_back(electrons_iter->ooEMOop());
      Electron_mHits.push_back(electrons_iter->missingHits());
      Electron_charge.push_back(electrons_iter->charge());
      Electron_tkiso.push_back(electrons_iter->trackIso());
      Electron_ecaliso.push_back(electrons_iter->ecalIso());
      Electron_hcaliso.push_back(electrons_iter->hcalIso());
      n_ele++;
    }

  n_pho = 0;

  for (auto photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
    Photon_pt.push_back(photons_iter->pt());
    Photon_eta.push_back(photons_iter->eta());
    Photon_phi.push_back(photons_iter->phi());
    Photon_m.push_back(photons_iter->m());
    Photon_sigmaietaieta.push_back(photons_iter->sigmaIetaIeta());
    Photon_HoE.push_back(photons_iter->hOverE());
    Photon_ecaliso.push_back(photons_iter->ecalIso());
    Photon_hcaliso.push_back(photons_iter->hcalIso());
    
    n_pho++;
  }

  vector<PseudoJet> fj_part;
  n_pfcand = 0;
    for (auto pfcands_iter = pfcandsH->begin(); pfcands_iter != pfcandsH->end(); ++pfcands_iter) {
      pfcandpt.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->pt())));
      pfcandeta.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->eta())));
      pfcandphi.push_back(MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->phi())));
    pdcandm.push_back(pfcands_iter->m());
    pfcandpdgid.push_back(pfcands_iter->pdgId());
    pfcandvertex.push_back(pfcands_iter->vertex());
    PseudoJet temp_jet = PseudoJet(0, 0, 0, 0);
    temp_jet.reset_PtYPhiM(pfcands_iter->pt(), pfcands_iter->eta(), pfcands_iter->phi(), pfcands_iter->m());
    temp_jet.set_user_index(pfcands_iter->pdgId());
    fj_part.push_back(temp_jet);

    n_pfcand++;
  } 

     n_mu=0;
for (auto muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
 	  Muon_pt.push_back(muons_iter->pt());
   	Muon_eta.push_back(muons_iter->eta());
   	Muon_phi.push_back(muons_iter->phi());
   	Muon_m.push_back(muons_iter->m());
   	Muon_ecaliso.push_back(muons_iter->ecalIso());
   	Muon_hcaliso.push_back(muons_iter->hcalIso());
   	Muon_trkiso.push_back(muons_iter->trk_chi2());
   	Muon_chi2.push_back(muons_iter->trk_ndof());
   	Muon_ndof.push_back(muons_iter->charge());
   	Muon_charge.push_back(muons_iter->trk_dxy());
   	Muon_dxy.push_back(muons_iter->trk_dz());
   	Muon_dz.push_back(muons_iter->nValidRecoMuonHits());
   	Muon_nvalidmuon_hits.push_back(muons_iter->nValidPixelHits());
   	Muon_nvalidpixelhits.push_back(muons_iter->nRecoMuonMatchedStations());
   	Muon_nmatchedstations.push_back(muons_iter->nTrackerLayersWithMeasurement());
    Muon_type.push_back(muons_iter->type());
    Muon_nvalidstriphits.push_back(muons_iter->nValidStripHits());
    Muon_trkqoverp.push_back(muons_iter->trk_qoverp());
    Muon_trklambda.push_back(muons_iter->trk_lambda());
    Muon_trkpt.push_back(muons_iter->trk_pt());
    Muon_trkphi.push_back(muons_iter->trk_phi());
    Muon_trketa.push_back(muons_iter->trk_eta());
    Muon_trkqoverperror.push_back(muons_iter->trk_dxyError());
    Muon_trklambdaerror.push_back(muons_iter->trk_dzError());
    Muon_trkpterror.push_back(muons_iter->trk_qoverpError());
    Muon_trkphierror.push_back(muons_iter->trk_lambdaError());
    Muon_trketaerror.push_back(muons_iter->trk_phiError());
    Muon_trkdszerror.push_back(muons_iter->trk_dsz());
    Muon_trkdsz.push_back(muons_iter->trk_dszError());
    n_mu++;
 }


  n_jet = 0;
   for (auto pfjets_iter = pfjetsH->begin(); pfjets_iter != pfjetsH->end(); ++pfjets_iter) {
    Jet_pt.push_back(pfjets_iter->pt());
    Jet_eta.push_back(pfjets_iter->eta());
    Jet_phi.push_back(pfjets_iter->phi());
    Jet_m.push_back(pfjets_iter->m());
    Jet_area.push_back(pfjets_iter->jetArea());
    Jet_chargedHadronEnergy.push_back(pfjets_iter->chargedHadronEnergy());
    Jet_neutralHadronEnergy.push_back(pfjets_iter->neutralHadronEnergy());
    Jet_photonEnergy.push_back(pfjets_iter->photonEnergy());
    Jet_electronEnergy.push_back(pfjets_iter->electronEnergy());
    Jet_muonEnergy.push_back(pfjets_iter->muonEnergy());
    Jet_HFHadronEnergy.push_back(pfjets_iter->HFHadronEnergy());
    Jet_HFEMEnergy.push_back(pfjets_iter->HFEMEnergy());
    Jet_HOEnergy.push_back(pfjets_iter->HOEnergy());
    
    Jet_chargedHadronMultiplicity.push_back(pfjets_iter->chargedHadronMultiplicity());
    Jet_neutralHadronMultiplicity.push_back(pfjets_iter->neutralHadronMultiplicity());
    Jet_photonMultiplicity.push_back(pfjets_iter->photonMultiplicity());
    Jet_electronMultiplicity.push_back(pfjets_iter->electronMultiplicity());
    Jet_muonMultiplicity.push_back(pfjets_iter->muonMultiplicity());
    Jet_HFHadronMultiplicity.push_back(pfjets_iter->HFHadronMultiplicity());
    Jet_HFEMMultiplicity.push_back(pfjets_iter->HFEMMultiplicity());
    Jet_csv.push_back(pfjets_iter->csv());
    Jet_mvaDiscriminator.push_back(pfjets_iter->mvaDiscriminator());
    // Jet_constituents.push_back(vector<int16_t>(pfjets_iter->constituents()));
    n_jet++;
  }

//FatJet stuff
  JetDefinition ak8_def = JetDefinition(antikt_algorithm, 0.8);
  double sd_z_cut = 0.10;
  double sd_beta = 0;
  SoftDrop sd_groomer = SoftDrop(sd_z_cut, sd_beta, 1.0);
  Filter trimmer = Filter(JetDefinition(kt_algorithm, 0.2), SelectorPtFractionMin(0.03));

  double beta = 1.0;
  Nsubjettiness nSub1 = Nsubjettiness(1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub2 = Nsubjettiness(2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub3 = Nsubjettiness(3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub4 = Nsubjettiness(4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
  Nsubjettiness nSub5 = Nsubjettiness(5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));

  EnergyCorrelatorN2 N2=EnergyCorrelatorN2(1.0);
  EnergyCorrelatorN3 N3=EnergyCorrelatorN3(1.0);

  fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  ClusterSequenceArea ak8_cs(fj_part, ak8_def, area_def);
  vector<PseudoJet> ak8_jets = sorted_by_pt(ak8_cs.inclusive_jets(100.0));

  for(auto &j: ak8_jets) {
    FatJet_area.push_back(j.area());
    FatJet_eta.push_back(j.pseudorapidity());
    FatJet_mass.push_back(j.m());
    
    PseudoJet sd_ak8 = sd_groomer(j);
    FatJet_msoftdrop.push_back(sd_ak8.m());
    
    PseudoJet trimmed_ak8 = trimmer(j);
    FatJet_mtrim.push_back(trimmed_ak8.m());
    
    FatJet_n2b1.push_back(N2(sd_ak8));
    FatJet_n3b1.push_back(N3(sd_ak8));
    FatJet_phi.push_back(j.phi_std());
    FatJet_pt.push_back(j.pt());
    FatJet_tau1.push_back(nSub1.result(j));
    FatJet_tau2.push_back(nSub2.result(j));
    FatJet_tau3.push_back(nSub3.result(j));
    FatJet_tau4.push_back(nSub4.result(j));
  }
  
 if (doL1) {
    l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
    	for( int r = 0; r<280; r++){
	string name ("empty");
	bool algoName_ = false;
	algoName_ = l1GtUtils_->getAlgNameFromBit(r,name);
	cout << "getAlgNameFromBit = " << algoName_  << endl;
	cout << "L1 bit number = " << r << " ; L1 bit name = " << name << endl;
	}
    for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
      bool l1htbit = 0;	
			
      l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
      //cout<<string(l1Seeds_[iseed])<<"  "<<l1htbit<<endl;
      l1Result_.push_back( l1htbit );
      }
 }


 tree->Fill();	
	
}

void ScoutingNanoAOD::clearVars(){
  Photon_pt.clear();
  Photon_eta.clear();
  Photon_phi.clear();
  Photon_m.clear();
  Photon_sigmaietaieta.clear();
  Photon_HoE.clear();
  Photon_ecaliso.clear();
  Photon_hcaliso.clear();
  Electron_pt.clear();
  Electron_eta.clear();
  Electron_phi.clear();
  Electron_m.clear();
  Electron_d0.clear();
  Electron_dz.clear();
  Electron_detain.clear();
  Electron_dphiin.clear();
  Electron_sigmaietaieta.clear();
  Electron_HoE.clear();
  Electron_ooEMOop.clear();
  Electron_mHits.clear();
  Electron_charge.clear();
  Electron_ecaliso.clear();
  Electron_hcaliso.clear();
  Electron_tkiso.clear();
  Muon_pt.clear();
  Muon_eta.clear();
  Muon_phi.clear();
  Muon_m.clear();
  Muon_ecaliso.clear();
  Muon_hcaliso.clear();
  Muon_trkiso.clear();
  Muon_chi2.clear();
  Muon_ndof.clear();
  Muon_charge.clear();
  Muon_dxy.clear();
  Muon_dz.clear();
  Muon_nvalidmuon_hits.clear();
  Muon_nvalidpixelhits.clear();
  Muon_nmatchedstations.clear();
  Muon_type.clear();
  Muon_nvalidstriphits.clear();
  Muon_trkqoverp.clear();
  Muon_trklambda.clear();
  Muon_trkpt.clear();
  Muon_trkphi.clear();
  Muon_trketa.clear();
  Muon_trkqoverperror.clear();
  Muon_trklambdaerror.clear();
  Muon_trkpterror.clear();
  Muon_trkphierror.clear();
  Muon_trketaerror.clear();
  Muon_trkdszerror.clear();
  Muon_trkdsz.clear();
  Jet_pt.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_m.clear();
  Jet_area.clear();
  Jet_chargedHadronEnergy.clear();
  Jet_neutralHadronEnergy.clear();
  Jet_photonEnergy.clear();
  Jet_electronEnergy.clear();
  Jet_muonEnergy.clear();
  Jet_HFHadronEnergy.clear();
  Jet_HFEMEnergy.clear();
  Jet_HOEnergy.clear();
  Jet_chargedHadronMultiplicity.clear();
  Jet_neutralHadronMultiplicity.clear();
  Jet_photonMultiplicity.clear();
  Jet_electronMultiplicity.clear();
  Jet_muonMultiplicity.clear();
  Jet_HFHadronMultiplicity.clear();
  Jet_HFEMMultiplicity.clear();
  Jet_csv.clear();
  Jet_mvaDiscriminator.clear();
  pfcandpt.clear();
  pfcandeta.clear();
  pfcandphi.clear();
  pdcandm.clear();
  pfcandpdgid.clear();
  pfcandvertex.clear();
  FatJet_area.clear();
  FatJet_eta.clear();
  FatJet_n2b1.clear();
  FatJet_n3b1.clear();
  FatJet_phi.clear();
  FatJet_pt.clear();
  FatJet_tau1.clear();
  FatJet_tau2.clear();
  FatJet_tau3.clear();
  FatJet_tau4.clear();
  FatJet_mass.clear();
  FatJet_msoftdrop.clear();
  FatJet_mtrim.clear();
}

void ScoutingNanoAOD::beginJob() {
  
}

void ScoutingNanoAOD::endJob() {
}

void ScoutingNanoAOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  // HLT paths

  triggerPathsVector.push_back("DST_DoubleMu1_noVtx_CaloScouting_v*");
  triggerPathsVector.push_back("DST_DoubleMu3_noVtx_CaloScouting_v*");
  triggerPathsVector.push_back("DST_DoubleMu3_noVtx_Mass10_PFScouting_v*");
  triggerPathsVector.push_back("DST_L1HTT_CaloScouting_PFScouting_v*");
  triggerPathsVector.push_back("DST_CaloJet40_CaloScouting_PFScouting_v*");
  triggerPathsVector.push_back("DST_HT250_CaloScouting_v*");
  triggerPathsVector.push_back("DST_HT410_PFScouting_v*");
  triggerPathsVector.push_back("DST_HT450_PFScouting_v*");

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

void ScoutingNanoAOD::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ScoutingNanoAOD::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void ScoutingNanoAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void ScoutingNanoAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScoutingNanoAOD);
