// -*- C++ -*-
//
// Package:    Analysis/Presel
// Class:      Presel
// 
/**\class Presel Presel.cc Analysis/Presel/plugins/Presel.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Robert Stringer
//         Created:  Thu, 01 Oct 2015 17:58:07 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisDataFormats/BoostedObjects/interface/GenParticleWithDaughters.h"
#include "AnalysisDataFormats/BoostedObjects/interface/Jet.h"
#include "AnalysisDataFormats/BoostedObjects/interface/ResolvedVjj.h"
#include "Analysis/VLQAna/interface/JetSelector.h"
#include "Analysis/VLQAna/interface/HT.h"
#include "Analysis/VLQAna/interface/VCandProducer.h"
#include "Analysis/VLQAna/interface/Utilities.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>

#include <sstream>
#include <bitset>
#include <boost/regex.hpp>

//bool sortByPt (const vlq::Jet& jet1, const vlq::Jet& jet2) {
//  return jet1.getPt() > jet2.getPt() ;  
//}
//
//bool sortByMass (const vlq::Jet& jet1, const vlq::Jet& jet2) {
//  return jet1.getMass() > jet2.getMass() ;  
//}
//
//bool sortByTrimmedMass (const vlq::Jet& jet1, const vlq::Jet& jet2) {
//  return jet1.getTrimmedMass() > jet2.getTrimmedMass() ;  
//}
//
//bool sortBySoftDropMass (const vlq::Jet& jet1, const vlq::Jet& jet2) {
//  return jet1.getSoftDropMass() > jet2.getSoftDropMass() ;  
//}
//
//bool sortByCSV (const vlq::Jet& jet1, const vlq::Jet& jet2) {
//  return jet1.getCSV() > jet2.getCSV() ;  
//}

//
// class declaration
//

class Presel : public edm::EDFilter {
  public:
    explicit Presel(const edm::ParameterSet&);
    ~Presel();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    std::string pn_;
    edm::InputTag l_trigName                     ; 
    edm::InputTag l_trigBit                      ; 
    edm::InputTag l_jetAK8Pt                     ; 
    edm::InputTag l_jetAK8Eta                    ; 
    edm::InputTag l_jetAK8Phi                    ; 
    edm::InputTag l_jetAK8Mass                   ; 
    edm::InputTag l_jetAK8FilteredMass           ; 
    edm::InputTag l_jetAK8TrimmedMass            ; 
    edm::InputTag l_jetAK8PrunedMass             ; 
    edm::InputTag l_jetAK8SoftDropMass           ; 
    edm::InputTag l_jetAK8Energy                 ; 
    edm::InputTag l_jetAK8Flavour                ; 
    edm::InputTag l_jetAK8CSV                    ; 
    edm::InputTag l_jetAK8JEC                    ; 
    edm::InputTag l_jetAK8Area                   ; 
    edm::InputTag l_jetAK8Tau1                   ;  
    edm::InputTag l_jetAK8Tau2                   ;  
    edm::InputTag l_jetAK8Tau3                   ;  
    edm::InputTag l_jetAK8nSubJets               ;  
    edm::InputTag l_jetAK8minmass                ;  
    edm::InputTag l_jetAK8VSubjetIndex0          ;  
    edm::InputTag l_jetAK8VSubjetIndex1          ;  
    edm::InputTag l_jetAK8TopSubjetIndex0        ; 
    edm::InputTag l_jetAK8TopSubjetIndex1        ; 
    edm::InputTag l_jetAK8TopSubjetIndex2        ; 
    edm::InputTag l_jetAK8TopSubjetIndex3        ; 
    edm::InputTag l_subjetAK8BDisc               ; 
    edm::InputTag l_subjetAK8Pt                  ; 
    edm::InputTag l_subjetAK8Eta                 ; 
    edm::InputTag l_subjetAK8Phi                 ; 
    edm::InputTag l_subjetAK8Mass                ; 
    edm::InputTag l_subjetCmsTopTagBDisc         ; 
    edm::InputTag l_subjetCmsTopTagPt            ; 
    edm::InputTag l_subjetCmsTopTagEta           ; 
    edm::InputTag l_subjetCmsTopTagPhi           ; 
    edm::InputTag l_subjetCmsTopTagMass          ; 
    edm::InputTag l_jetAK4Pt                     ; 
    edm::InputTag l_jetAK4Eta                    ; 
    edm::InputTag l_jetAK4Phi                    ; 
    edm::InputTag l_jetAK4Mass                   ; 
    edm::InputTag l_jetAK4Energy                 ; 
    edm::InputTag l_jetAK4Flavour                ; 
    edm::InputTag l_jetAK4CSV                    ; 
    edm::InputTag l_jetAK4JEC                    ; 
    edm::InputTag l_jetAK4nHadEnergy             ;
    edm::InputTag l_jetAK4nEMEnergy              ;
    edm::InputTag l_jetAK4HFHadronEnergy         ;
    edm::InputTag l_jetAK4cHadEnergy             ;
    edm::InputTag l_jetAK4cEMEnergy              ;
    edm::InputTag l_jetAK4numDaughters           ;
    edm::InputTag l_jetAK4cMultip                ;
    edm::InputTag l_jetAK4Y                      ;
    edm::InputTag l_jetAK4Area                   ; 
    edm::InputTag l_HbbCands                     ; 
    std::vector<std::string> hltPaths_           ; 
    std::vector<std::string> metFilters_         ; 
    edm::ParameterSet GenHSelParams_             ; 
    edm::ParameterSet AK4JetSelParams_           ; 
    edm::ParameterSet AK4TrigJetSelParams_       ; 
    edm::ParameterSet BTaggedLooseAK4SelParams_  ; 
    edm::ParameterSet BTaggedMediumAK4SelParams_  ; 
    edm::ParameterSet BTaggedTightAK4SelParams_  ; 
    edm::ParameterSet AK8JetSelParams_           ; 
    edm::ParameterSet TJetSelParams_             ; 
    edm::ParameterSet HJetSelParams_             ; 
    edm::ParameterSet WJetSelParams_             ; 
    double ak8jetsPtMin_                         ;
    double ak8jetsEtaMax_                        ; 
    double ak4jetsPtMin_                         ;
    double ak4jetsEtaMax_                        ; 
    double HTMin_                                ; 
    edm::InputTag l_metFiltersName               ; 
    edm::InputTag l_metFiltersBit                ; 
    edm::InputTag l_hbheNoiseFilter              ; 
    edm::InputTag l_vtxRho                       ; 
    edm::InputTag l_vtxZ                         ; 
    edm::InputTag l_vtxNdf                       ; 
    bool _isData				 ;
    string skipEvents				 ; 
    string weightFileName_			 ;   
    std::vector<std::string> weights_		 ;

    edm::Service<TFileService> fs                ; 
    std::map<std::string, TH1D*> h1_             ; 
    std::map<std::string, TH2D*> h2_             ; 

    TFile * weightFile;
    TGraphAsymmErrors * wGraph[10]; 
    std::vector<string> weightnames;

      // ----------member data ---------------------------
};

using namespace std; 

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Presel::Presel(const edm::ParameterSet& iConfig) :
  pn_                     (iConfig.getParameter<string>            ("processName")),
  l_trigName              (iConfig.getParameter<edm::InputTag>     ("trigNameLabel")),
  l_trigBit               (iConfig.getParameter<edm::InputTag>     ("trigBitLabel")),
  l_jetAK8Pt              (iConfig.getParameter<edm::InputTag>     ("jetAK8PtLabel")),  
  l_jetAK8Eta             (iConfig.getParameter<edm::InputTag>     ("jetAK8EtaLabel")),  
  l_jetAK8Phi             (iConfig.getParameter<edm::InputTag>     ("jetAK8PhiLabel")),  
  l_jetAK8Mass            (iConfig.getParameter<edm::InputTag>     ("jetAK8MassLabel")),  
  l_jetAK8FilteredMass    (iConfig.getParameter<edm::InputTag>     ("jetAK8FilteredMassLabel")),  
  l_jetAK8TrimmedMass     (iConfig.getParameter<edm::InputTag>     ("jetAK8TrimmedMassLabel")),  
  l_jetAK8PrunedMass      (iConfig.getParameter<edm::InputTag>     ("jetAK8PrunedMassLabel")),  
  l_jetAK8SoftDropMass    (iConfig.getParameter<edm::InputTag>     ("jetAK8SoftDropMassLabel")),  
  l_jetAK8Energy          (iConfig.getParameter<edm::InputTag>     ("jetAK8EnergyLabel")),  
  l_jetAK8Flavour         (iConfig.getParameter<edm::InputTag>     ("jetAK8FlavourLabel")),  
  l_jetAK8CSV             (iConfig.getParameter<edm::InputTag>     ("jetAK8CSVLabel")),  
  l_jetAK8JEC             (iConfig.getParameter<edm::InputTag>     ("jetAK8JECLabel")), 
  l_jetAK8Area            (iConfig.getParameter<edm::InputTag>     ("jetAK8AreaLabel")),
  l_jetAK8Tau1            (iConfig.getParameter<edm::InputTag>     ("jetAK8Tau1Label")),  
  l_jetAK8Tau2            (iConfig.getParameter<edm::InputTag>     ("jetAK8Tau2Label")),  
  l_jetAK8Tau3            (iConfig.getParameter<edm::InputTag>     ("jetAK8Tau3Label")),  
  l_jetAK8nSubJets        (iConfig.getParameter<edm::InputTag>     ("jetAK8nSubJetsLabel")),  
  l_jetAK8minmass         (iConfig.getParameter<edm::InputTag>     ("jetAK8minmassLabel")),  
  l_jetAK8VSubjetIndex0   (iConfig.getParameter<edm::InputTag>     ("jetAK8VSubjetIndex0Label")),  
  l_jetAK8VSubjetIndex1   (iConfig.getParameter<edm::InputTag>     ("jetAK8VSubjetIndex1Label")),  
  l_jetAK8TopSubjetIndex0 (iConfig.getParameter<edm::InputTag>     ("jetAK8TopSubjetIndex0Label")), 
  l_jetAK8TopSubjetIndex1 (iConfig.getParameter<edm::InputTag>     ("jetAK8TopSubjetIndex1Label")), 
  l_jetAK8TopSubjetIndex2 (iConfig.getParameter<edm::InputTag>     ("jetAK8TopSubjetIndex2Label")), 
  l_jetAK8TopSubjetIndex3 (iConfig.getParameter<edm::InputTag>     ("jetAK8TopSubjetIndex3Label")), 
  l_subjetAK8BDisc        (iConfig.getParameter<edm::InputTag>     ("subjetAK8BDiscLabel")), 
  l_subjetAK8Pt           (iConfig.getParameter<edm::InputTag>     ("subjetAK8PtLabel")), 
  l_subjetAK8Eta          (iConfig.getParameter<edm::InputTag>     ("subjetAK8EtaLabel")), 
  l_subjetAK8Phi          (iConfig.getParameter<edm::InputTag>     ("subjetAK8PhiLabel")), 
  l_subjetAK8Mass         (iConfig.getParameter<edm::InputTag>     ("subjetAK8MassLabel")), 
  l_subjetCmsTopTagBDisc  (iConfig.getParameter<edm::InputTag>     ("subjetCmsTopTagBDiscLabel")), 
  l_subjetCmsTopTagPt     (iConfig.getParameter<edm::InputTag>     ("subjetCmsTopTagPtLabel")), 
  l_subjetCmsTopTagEta    (iConfig.getParameter<edm::InputTag>     ("subjetCmsTopTagEtaLabel")), 
  l_subjetCmsTopTagPhi    (iConfig.getParameter<edm::InputTag>     ("subjetCmsTopTagPhiLabel")), 
  l_subjetCmsTopTagMass   (iConfig.getParameter<edm::InputTag>     ("subjetCmsTopTagMassLabel")), 
  l_jetAK4Pt              (iConfig.getParameter<edm::InputTag>     ("jetAK4PtLabel")),  
  l_jetAK4Eta             (iConfig.getParameter<edm::InputTag>     ("jetAK4EtaLabel")),  
  l_jetAK4Phi             (iConfig.getParameter<edm::InputTag>     ("jetAK4PhiLabel")),  
  l_jetAK4Mass            (iConfig.getParameter<edm::InputTag>     ("jetAK4MassLabel")),  
  l_jetAK4Energy          (iConfig.getParameter<edm::InputTag>     ("jetAK4EnergyLabel")),  
  l_jetAK4Flavour         (iConfig.getParameter<edm::InputTag>     ("jetAK4FlavourLabel")),  
  l_jetAK4CSV             (iConfig.getParameter<edm::InputTag>     ("jetAK4CSVLabel")),  
  l_jetAK4JEC             (iConfig.getParameter<edm::InputTag>     ("jetAK4JECLabel")),
  l_jetAK4nHadEnergy      (iConfig.getParameter<edm::InputTag>     ("jetAK4nHadEnergyLabel")),
  l_jetAK4nEMEnergy       (iConfig.getParameter<edm::InputTag>     ("jetAK4nEMEnergyLabel")),
  l_jetAK4HFHadronEnergy  (iConfig.getParameter<edm::InputTag>     ("jetAK4HFHadronEnergyLabel")),
  l_jetAK4cHadEnergy      (iConfig.getParameter<edm::InputTag>     ("jetAK4cHadEnergyLabel")),
  l_jetAK4cEMEnergy       (iConfig.getParameter<edm::InputTag>     ("jetAK4cEMEnergyLabel")),
  l_jetAK4numDaughters    (iConfig.getParameter<edm::InputTag>     ("jetAK4numDaughtersLabel")),
  l_jetAK4cMultip         (iConfig.getParameter<edm::InputTag>     ("jetAK4cMultipLabel")),
  l_jetAK4Y               (iConfig.getParameter<edm::InputTag>     ("jetAK4YLabel")),
  l_jetAK4Area            (iConfig.getParameter<edm::InputTag>     ("jetAK4AreaLabel")),
  l_HbbCands              (iConfig.getParameter<edm::InputTag>     ("HbbCandsLabel")),
  hltPaths_               (iConfig.getParameter<vector<string>>    ("hltPaths")), 
  metFilters_             (iConfig.getParameter<vector<string>>    ("metFilters")), 
  GenHSelParams_          (iConfig.getParameter<edm::ParameterSet> ("GenHSelParams")),
  AK4JetSelParams_        (iConfig.getParameter<edm::ParameterSet> ("AK4JetSelParams")),
  AK4TrigJetSelParams_    (iConfig.getParameter<edm::ParameterSet> ("AK4TrigJetSelParams")),
  BTaggedLooseAK4SelParams_ (iConfig.getParameter<edm::ParameterSet> ("BTaggedLooseAK4SelParams")),
  BTaggedMediumAK4SelParams_ (iConfig.getParameter<edm::ParameterSet> ("BTaggedMediumAK4SelParams")),
  BTaggedTightAK4SelParams_ (iConfig.getParameter<edm::ParameterSet> ("BTaggedTightAK4SelParams")),
  AK8JetSelParams_        (iConfig.getParameter<edm::ParameterSet> ("AK8JetSelParams")),
  TJetSelParams_          (iConfig.getParameter<edm::ParameterSet>  ("TJetSelParams")),
  HJetSelParams_          (iConfig.getParameter<edm::ParameterSet>  ("HJetSelParams")),
  WJetSelParams_          (iConfig.getParameter<edm::ParameterSet>  ("WJetSelParams")),
  ak8jetsPtMin_           (iConfig.getParameter<double>            ("ak8jetsPtMin")),
  ak8jetsEtaMax_          (iConfig.getParameter<double>            ("ak8jetsEtaMax")), 
  ak4jetsPtMin_           (iConfig.getParameter<double>            ("ak4jetsPtMin")),
  ak4jetsEtaMax_          (iConfig.getParameter<double>            ("ak4jetsEtaMax")), 
  HTMin_                  (iConfig.getParameter<double>            ("HTMin")),
  l_metFiltersName        (iConfig.getParameter<edm::InputTag>     ("metFiltersNameLabel")),
  l_metFiltersBit         (iConfig.getParameter<edm::InputTag>     ("metFiltersBitLabel")),
  l_hbheNoiseFilter       (iConfig.getParameter<edm::InputTag>     ("hbheNoiseFilterLabel")),
  l_vtxRho                (iConfig.getParameter<edm::InputTag>     ("vtxRhoLabel")),  
  l_vtxZ                  (iConfig.getParameter<edm::InputTag>     ("vtxZLabel")),  
  l_vtxNdf                (iConfig.getParameter<edm::InputTag>     ("vtxNdfLabel")),  
  _isData		  (iConfig.getParameter<bool>		   ("isData")),
  skipEvents		  (iConfig.getParameter<string>		   ("skipEvents")),
  weightFileName_ 	  (iConfig.getParameter<string>	   	   ("weightFileName")),                      
  weights_ 		  (iConfig.getParameter<vector<string>>    ("weights"))

{
  produces<unsigned>("ngoodAK4Jets");
  produces<unsigned>("ngoodAK8Jets");
  produces<unsigned>("nbtaggedlooseAK4");
  produces<unsigned>("nbtaggedmediumAK4");
  produces<unsigned>("nbtaggedtightAK4");
  produces<unsigned>("nTJets");
  produces<unsigned>("nHJets");
  produces<unsigned>("nWJets");
  produces<double>("htak4jets");
  produces<double>("htak8jets");
  produces<double>("htak4trigjets");
  produces<double>("maxetaak4");
  produces<double>("MassLeading2AK8");
  produces<double>("DeltaEtaLeading2AK8");
  produces<double>("pt1stAK8");
  produces<double>("pt2ndAK8");
  produces<double>("mass1stAK8");
  produces<double>("mass2ndAK8");
  produces<vector<unsigned> >("ak4goodjets");
  produces<vector<unsigned> >("ak8goodjets");
  produces<vector<unsigned>>("bjetIdxs");
  produces<vector<unsigned>>("tjetIdxs");
  produces<vector<unsigned>>("hjetIdxs");
  produces<vector<unsigned>>("wjetIdxs");

  produces<vector<string> >("WeightName");
  produces<vector<double> >("EvtWeight"); 

 

  //register your products
  /* Examples
     produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");

  //if you want to put into the Run
  produces<ExampleData2,InRun>();
  */
  //now do what ever other initialization is needed

}


Presel::~Presel()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
bool Presel::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;

  typedef Handle <vector<string>> hstring ; 
  typedef Handle <vector<float>> hfloat ; 
  typedef Handle <vector<int>> hint ; 

  hstring h_trigName             ; evt.getByLabel (l_trigName               , h_trigName             );
  hfloat  h_trigBit              ; evt.getByLabel (l_trigBit                , h_trigBit              ); 
  hstring h_metFiltersName       ; evt.getByLabel (l_metFiltersName         , h_metFiltersName       );
  hfloat  h_metFiltersBit        ; evt.getByLabel (l_metFiltersBit          , h_metFiltersBit        ); 
 
  Handle <bool> h_hbheNoiseFilter ; evt.getByLabel (l_hbheNoiseFilter, h_hbheNoiseFilter);

  hfloat  h_vtxRho               ; evt.getByLabel (l_vtxRho                 , h_vtxRho               );
  hfloat  h_vtxZ                 ; evt.getByLabel (l_vtxZ                   , h_vtxZ                 );
  hint    h_vtxNdf               ; evt.getByLabel (l_vtxNdf                 , h_vtxNdf               );
  hfloat  h_jetAK8Pt             ; evt.getByLabel (l_jetAK8Pt               , h_jetAK8Pt             );
  hfloat  h_jetAK8Eta            ; evt.getByLabel (l_jetAK8Eta              , h_jetAK8Eta            );
  hfloat  h_jetAK8Phi            ; evt.getByLabel (l_jetAK8Phi              , h_jetAK8Phi            );
  hfloat  h_jetAK8Mass           ; evt.getByLabel (l_jetAK8Mass             , h_jetAK8Mass           );
  hfloat  h_jetAK8FilteredMass   ; evt.getByLabel (l_jetAK8FilteredMass     , h_jetAK8FilteredMass   );
  hfloat  h_jetAK8TrimmedMass    ; evt.getByLabel (l_jetAK8TrimmedMass      , h_jetAK8TrimmedMass    );
  hfloat  h_jetAK8PrunedMass     ; evt.getByLabel (l_jetAK8PrunedMass       , h_jetAK8PrunedMass     );
  hfloat  h_jetAK8SoftDropMass   ; evt.getByLabel (l_jetAK8SoftDropMass     , h_jetAK8SoftDropMass   );
  hfloat  h_jetAK8Energy         ; evt.getByLabel (l_jetAK8Energy           , h_jetAK8Energy         );
  hfloat  h_jetAK8Flavour        ; evt.getByLabel (l_jetAK8Flavour          , h_jetAK8Flavour        );
  hfloat  h_jetAK8CSV            ; evt.getByLabel (l_jetAK8CSV              , h_jetAK8CSV            );
  hfloat  h_jetAK8JEC            ; evt.getByLabel (l_jetAK8JEC              , h_jetAK8JEC            );
  hfloat  h_jetAK8Area           ; evt.getByLabel (l_jetAK8Area             , h_jetAK8Area           );
  hfloat  h_jetAK8Tau1           ; evt.getByLabel (l_jetAK8Tau1             , h_jetAK8Tau1           ); 
  hfloat  h_jetAK8Tau2           ; evt.getByLabel (l_jetAK8Tau2             , h_jetAK8Tau2           ); 
  hfloat  h_jetAK8Tau3           ; evt.getByLabel (l_jetAK8Tau3             , h_jetAK8Tau3           ); 
  hfloat  h_jetAK8nSubJets       ; evt.getByLabel (l_jetAK8nSubJets         , h_jetAK8nSubJets       ); 
  hfloat  h_jetAK8minmass        ; evt.getByLabel (l_jetAK8minmass          , h_jetAK8minmass        ); 
  hfloat  h_jetAK8VSubjetIndex0  ; evt.getByLabel (l_jetAK8VSubjetIndex0    , h_jetAK8VSubjetIndex0  );  
  hfloat  h_jetAK8VSubjetIndex1  ; evt.getByLabel (l_jetAK8VSubjetIndex1    , h_jetAK8VSubjetIndex1  );  
  hfloat  h_jetAK8TopSubjetIndex0; evt.getByLabel (l_jetAK8TopSubjetIndex0  , h_jetAK8TopSubjetIndex0); 
  hfloat  h_jetAK8TopSubjetIndex1; evt.getByLabel (l_jetAK8TopSubjetIndex1  , h_jetAK8TopSubjetIndex1); 
  hfloat  h_jetAK8TopSubjetIndex2; evt.getByLabel (l_jetAK8TopSubjetIndex2  , h_jetAK8TopSubjetIndex2); 
  hfloat  h_jetAK8TopSubjetIndex3; evt.getByLabel (l_jetAK8TopSubjetIndex3  , h_jetAK8TopSubjetIndex3); 
  hfloat  h_subjetAK8BDisc       ; evt.getByLabel (l_subjetAK8BDisc         , h_subjetAK8BDisc       ); 
  hfloat  h_subjetAK8Pt          ; evt.getByLabel (l_subjetAK8Pt            , h_subjetAK8Pt          ); 
  hfloat  h_subjetAK8Eta         ; evt.getByLabel (l_subjetAK8Eta           , h_subjetAK8Eta         ); 
  hfloat  h_subjetAK8Phi         ; evt.getByLabel (l_subjetAK8Phi           , h_subjetAK8Phi         ); 
  hfloat  h_subjetAK8Mass        ; evt.getByLabel (l_subjetAK8Mass          , h_subjetAK8Mass        ); 
  hfloat  h_subjetCmsTopTagBDisc ; evt.getByLabel (l_subjetCmsTopTagBDisc   , h_subjetCmsTopTagBDisc ); 
  hfloat  h_subjetCmsTopTagPt    ; evt.getByLabel (l_subjetCmsTopTagPt      , h_subjetCmsTopTagPt    ); 
  hfloat  h_subjetCmsTopTagEta   ; evt.getByLabel (l_subjetCmsTopTagEta     , h_subjetCmsTopTagEta   ); 
  hfloat  h_subjetCmsTopTagPhi   ; evt.getByLabel (l_subjetCmsTopTagPhi     , h_subjetCmsTopTagPhi   ); 
  hfloat  h_subjetCmsTopTagMass  ; evt.getByLabel (l_subjetCmsTopTagMass    , h_subjetCmsTopTagMass  ); 
  hfloat  h_jetAK4Pt             ; evt.getByLabel (l_jetAK4Pt               , h_jetAK4Pt             );
  hfloat  h_jetAK4Eta            ; evt.getByLabel (l_jetAK4Eta              , h_jetAK4Eta            );
  hfloat  h_jetAK4Phi            ; evt.getByLabel (l_jetAK4Phi              , h_jetAK4Phi            );
  hfloat  h_jetAK4Mass           ; evt.getByLabel (l_jetAK4Mass             , h_jetAK4Mass           );
  hfloat  h_jetAK4Energy         ; evt.getByLabel (l_jetAK4Energy           , h_jetAK4Energy         );
  hfloat  h_jetAK4Flavour        ; evt.getByLabel (l_jetAK4Flavour          , h_jetAK4Flavour        );
  hfloat  h_jetAK4CSV            ; evt.getByLabel (l_jetAK4CSV              , h_jetAK4CSV            );
  hfloat  h_jetAK4JEC            ; evt.getByLabel (l_jetAK4JEC              , h_jetAK4JEC            );
  hfloat  h_jetAK4nHadEnergy     ; evt.getByLabel (l_jetAK4nHadEnergy       , h_jetAK4nHadEnergy     );
  hfloat  h_jetAK4nEMEnergy      ; evt.getByLabel (l_jetAK4nEMEnergy        , h_jetAK4nEMEnergy      );
  hfloat  h_jetAK4HFHadronEnergy ; evt.getByLabel (l_jetAK4HFHadronEnergy   , h_jetAK4HFHadronEnergy );
  hfloat  h_jetAK4cHadEnergy     ; evt.getByLabel (l_jetAK4cHadEnergy       , h_jetAK4cHadEnergy     );
  hfloat  h_jetAK4cEMEnergy      ; evt.getByLabel (l_jetAK4cEMEnergy        , h_jetAK4cEMEnergy      );
  hfloat  h_jetAK4numDaughters   ; evt.getByLabel (l_jetAK4numDaughters     , h_jetAK4numDaughters   );
  hfloat  h_jetAK4cMultip        ; evt.getByLabel (l_jetAK4cMultip          , h_jetAK4cMultip        );
  hfloat  h_jetAK4Y              ; evt.getByLabel (l_jetAK4Y                , h_jetAK4Y              );
  hfloat  h_jetAK4Area           ; evt.getByLabel (l_jetAK4Area             , h_jetAK4Area           );

 // if(!bData) {
 //     Handle<vlq::GenParticleWithDaughtersCollection> h_HbbCands ; evt.getByLabel (l_HbbCands , h_HbbCands );
 //     vlq::GenParticleWithDaughtersCollection::const_iterator ihbb ;
 // }
  
  //// Preselection HLT
  // Get all trig names
  //for ( vector<string>::const_iterator it = (h_trigName.product())->begin(); it != (h_trigName.product())->end(); ++it) {
  //  cout << *it << endl ; 
  //}

  //return false ; 
  //Skip event

  unsigned long run, event, lumi;

  run = evt.id().run();
  event = evt.id().event();
  lumi = evt.luminosityBlock();

 // cout << run << " " << lumi << " " << event << endl;  

  vector<string> events;
  istringstream f(skipEvents);
  string s;    
  while (getline(f, s, ';')) {
      events.push_back(s);
  }
  
  

  for( auto& it : events ) {
    vector<string> theevt;
    istringstream f(it);
    string s;    
    while (getline(f, s, ',')) {
      theevt.push_back(s);
    }
//    cout << theevt[0] << " " << theevt[1] << " " << theevt[2] << endl;
//    cout << "Checking " << run << " " << lumi << " " << event << endl;

    if (stoul(theevt[0].c_str(),nullptr,0) == run && stoul(theevt[1].c_str(),nullptr,0) == lumi && stoul(theevt[2].c_str(),nullptr,0) == event ) {
	cout << "Skipping " << run << " " << lumi << " " << event << endl;
	return false;
    }
  }

  unsigned int hltdecisions(0) ; 
  for ( unsigned int i = 0 ; i < h_trigName.product()->size() ; i++ ) {
//    vector<string>::const_iterator it = find( (h_trigName.product())->begin(), (h_trigName.product())->end(), myhltpath) ; 

    //Support regex pattern matching: wildcards . and .*
    //
     bool bPass = 0; 
     for( const string &myhltpath : hltPaths_ ) {

	string hltname = h_trigName.product()->at(i);
	boost::regex pattern(myhltpath);


//	cout << hltname << "  " << pattern << "  " << result << endl;

	bPass |= boost::regex_match(hltname , pattern);
        hltdecisions += (pow(2,i) && bPass);
     }

//    if ( it != (h_trigName.product())->end() ) {
//      std::string hltname = (h_trigName.product()) -> at (it - (h_trigName.product())->begin()) ; 
//	cout << "hlt" << hltname << endl;

      //cout << " HLT path " << *it << " decision = " << int((h_trigBit.product())->at(it - (h_trigName.product())->begin())) << endl ;
  }
  
  if (hltdecisions == 0) return false ; 
 
  
  vlq::JetCollection goodAK8Jets, goodAK4Jets, btaggedlooseAK4, btaggedmediumAK4, btaggedtightAK4, goodAK4TrigJets ;
  vector<unsigned> ak4selIdxs, ak8selIdxs, bjetIdxs;

  //// Store good AK8 jets
  JetSelector ak8jetsel(AK8JetSelParams_) ;
  for ( unsigned ijet = 0; ijet < (h_jetAK8Pt.product())->size(); ++ijet) {
    bool retak8jetsel = false ; 
    if (ak8jetsel(evt, ijet,retak8jetsel) == 0) { 
      LogDebug("Presel") << " ak8 jet with pt = " << (h_jetAK8Pt.product())->at(ijet) << " fail jet sel\n" ; 
      continue ;
    }
    TLorentzVector jetP4 ;
    jetP4.SetPtEtaPhiM( (h_jetAK8Pt.product())->at(ijet), 
        (h_jetAK8Eta.product())->at(ijet), 
        (h_jetAK8Phi.product())->at(ijet), 
        (h_jetAK8Mass.product())->at(ijet) ) ;
    vlq::Jet jet;
    jet.setP4(jetP4) ; 
    jet.setFilteredMass((h_jetAK8FilteredMass.product())->at(ijet)) ; 
    jet.setTrimmedMass((h_jetAK8TrimmedMass.product())->at(ijet)) ; 
    jet.setPrunedMass((h_jetAK8PrunedMass.product())->at(ijet)) ; 
    jet.setSoftDropMass((h_jetAK8SoftDropMass.product())->at(ijet)) ; 
    jet.setCSV((h_jetAK8CSV.product())->at(ijet)) ;
    goodAK8Jets.push_back(jet) ;
    ak8selIdxs.push_back(ijet);
  }


  //// Store good AK4 jets 
  JetSelector ak4jetsel(AK4JetSelParams_) ;
  JetSelector btaggedlooseak4sel(BTaggedLooseAK4SelParams_) ;
  JetSelector btaggedmediumak4sel(BTaggedMediumAK4SelParams_) ;
  JetSelector btaggedtightak4sel(BTaggedTightAK4SelParams_) ;
  bool retak4jetsel = false ; 
  bool retbtaggedlooseak4sel = false ; 
  bool retbtaggedmediumak4sel = false ; 
  bool retbtaggedtightak4sel = false ; 
  for ( unsigned ijet = 0; ijet < (h_jetAK4Pt.product())->size(); ++ijet) {
    retak4jetsel = false ;
    retbtaggedlooseak4sel = false ;
    retbtaggedmediumak4sel = false ;
    retbtaggedtightak4sel = false;
    if (ak4jetsel(evt, ijet,retak4jetsel) == 0) { 
      LogDebug("Presel") << " ak4 jet with pt = " << (h_jetAK4Pt.product())->at(ijet) << " fail jet sel\n" ; 
      continue ;
    }
    TLorentzVector jetP4 ;
    jetP4.SetPtEtaPhiM( (h_jetAK4Pt.product())->at(ijet), 
        (h_jetAK4Eta.product())->at(ijet), 
        (h_jetAK4Phi.product())->at(ijet), 
        (h_jetAK4Mass.product())->at(ijet) ) ;
    vlq::Jet jet;  
    jet.setP4(jetP4) ;
    jet.setCSV((h_jetAK4CSV.product())->at(ijet)) ;
    goodAK4Jets.push_back(jet) ;
    ak4selIdxs.push_back(ijet);
    if ( btaggedlooseak4sel(evt, ijet,retbtaggedlooseak4sel) != 0 ) {
      bjetIdxs.push_back(ijet) ; 
      btaggedlooseAK4.push_back(jet) ; 
  //    cout << "index: "<< ijet << "   CSV: "<< (h_jetAK4CSV.product())->at(ijet) << endl;
    }
    if ( btaggedmediumak4sel(evt, ijet,retbtaggedmediumak4sel) != 0 ) {
      btaggedmediumAK4.push_back(jet) ; 
    }
    if ( btaggedtightak4sel(evt, ijet,retbtaggedtightak4sel) != 0 ) {
      btaggedtightAK4.push_back(jet) ; 
    }

  }

  // Store good AK4 Jets (Trigger Definition)
  JetSelector ak4trigjetsel(AK4TrigJetSelParams_) ;
  bool retak4trigjetsel = false;
  for ( unsigned ijet = 0; ijet < (h_jetAK4Pt.product())->size(); ++ijet) {
    retak4trigjetsel = false;
    if (ak4trigjetsel(evt, ijet, retak4trigjetsel ) == 0)
       continue;
    TLorentzVector jetP4 ;
    jetP4.SetPtEtaPhiM( (h_jetAK4Pt.product())->at(ijet),
        (h_jetAK4Eta.product())->at(ijet),
        (h_jetAK4Phi.product())->at(ijet),
        (h_jetAK4Mass.product())->at(ijet) ) ;
    vlq::Jet jet;
    jet.setP4(jetP4) ;
    jet.setCSV((h_jetAK4CSV.product())->at(ijet)) ;
    goodAK4TrigJets.push_back(jet) ;
  }
  

  HT htak4trig(goodAK4TrigJets);
  HT htak4(goodAK4Jets) ; 
  HT htak8(goodAK8Jets) ; 


  std::vector<double> weight;

  for (unsigned i = 0 ; i < weights_.size() && i < 10 ; i++) {
      if (wGraph[i]) {
//	cout << "HT: " << htak4.getHT() << endl;
//	cout << "weight: " << wGraph[i]->Eval(htak4.getHT()) << endl;
	
	weight.push_back(wGraph[i]->Eval(htak4.getHT()));
      }
      else 
	cout << "Can't find hist" << endl;
  }
  

  //Fill no-cut histos

  h1_["nak8_nocuts"] -> Fill(goodAK8Jets.size()) ; 
  h1_["nak4_nocuts"] -> Fill(goodAK4Jets.size()) ; 
  h1_["nbloose_nocuts"] -> Fill(btaggedlooseAK4.size()) ; 
  h1_["nbmedium_nocuts"] -> Fill(btaggedmediumAK4.size()) ; 
  h1_["nbtight_nocuts"] -> Fill(btaggedtightAK4.size()) ; 
  h1_["htak4_nocuts"] -> Fill(htak4.getHT());
  h1_["htak4trig_nocuts"] -> Fill(htak4trig.getHT());


  for ( vlq::Jet& jet : goodAK4Jets) {
    h1_["ak4_pt_nocuts"] -> Fill(jet.getPt());
    h1_["ak4_eta_nocuts"] -> Fill(jet.getEta());
  }


  //// Preselection at least four b-tagged AK4 jet 

   if ( goodAK4Jets.size() < 4) return false;

  //// Preselection one AK8 jet 
  if ( goodAK8Jets.size() < 1 ) return false ; 
  //if ( goodAK8Jets.size() > 1 && ( goodAK8Jets.at(0).getPt() < 300 || goodAK8Jets.at(1).getPt() < 220.) )  return false ; 


  //// Make W, top and H jets 
  vector<unsigned> seltjets, selhjets, selwjets;
  vlq::JetCollection tjets, hjets, wjets ; 
  JetSelector tjetsel(TJetSelParams_) ;
  JetSelector hjetsel(HJetSelParams_) ;
  JetSelector wjetsel(WJetSelParams_) ;
  bool rettjetsel = false ;
  bool rethjetsel = false ;
  bool retwjetsel = false ;
  for ( unsigned  &ijet :  ak8selIdxs) {
    TLorentzVector jetP4 ;
    vector<float>drhjet_hpart, drhjet_bpart, drhjet_bbarpart ; 
    jetP4.SetPtEtaPhiM( (h_jetAK8Pt.product())->at(ijet), 
        (h_jetAK8Eta.product())->at(ijet), 
        (h_jetAK8Phi.product())->at(ijet), 
        (h_jetAK8Mass.product())->at(ijet) ) ;

/*    for ( ihbb = h_HbbCands.product()->begin(); ihbb != h_HbbCands.product()->end(); ++ihbb ) {
      TLorentzVector hp4 = (ihbb->getMom()).getP4() ; 
      drhjet_hpart.push_back(hp4.DeltaR(jetP4)) ; 
      vlq::GenParticleCollection bs =  ( ihbb->getDaughters() ) ; 
      for ( vlq::GenParticle& b : bs ) {
        TLorentzVector bp4 = b.getP4() ; 
        if ( b.getPdgID() == 5 ) {
          drhjet_bpart.push_back(bp4.DeltaR(jetP4)) ; 
        }
        if ( b.getPdgID() == -5 ) {
          drhjet_bbarpart.push_back(bp4.DeltaR(jetP4)) ; 
        }
      }
    }
*
*/

  try {
    std::vector<float>::iterator iminh    = std::min_element(drhjet_hpart.begin(), drhjet_hpart.end());
    std::vector<float>::iterator iminb    = std::min_element(drhjet_bpart.begin(), drhjet_bpart.end());
    std::vector<float>::iterator iminbbar = std::min_element(drhjet_bbarpart.begin(), drhjet_bbarpart.end());
    if ( iminh != drhjet_hpart.end() && iminb != drhjet_bpart.end() && iminbbar != drhjet_bbarpart.end() )
      if ( *iminh < 0.8 && *iminb < 0.8 && *iminbbar < 0.8 ) 
        h1_["ptak8"]->Fill(jetP4.Pt()) ; 
    vlq::Jet jet; 
    jet.setP4(jetP4) ;
    jet.setCSV((h_jetAK8CSV.product())->at(ijet)) ;
    rettjetsel = false ;
    if (tjetsel(evt, ijet,rettjetsel) == true ) { 
      tjets.push_back(jet) ; 
      seltjets.push_back(ijet) ; 
    }
    rethjetsel = false ;
    if (hjetsel(evt, ijet,rethjetsel) == true ) { 
     hjets.push_back(jet) ; 
     selhjets.push_back(ijet) ; 
     h1_["csvhjets"] -> Fill((h_jetAK8CSV.product())->at(ijet)) ; 
     if ( iminh != drhjet_hpart.end() && iminb != drhjet_bpart.end() && iminbbar != drhjet_bbarpart.end() )
       if ( *iminh < 0.8  && *iminb < 0.8 && *iminbbar < 0.8 ) 
         h1_["pthjets"]->Fill(jetP4.Pt()) ;
    } 
    retwjetsel = false ;
    if (wjetsel(evt, ijet,retwjetsel) == true ) { 
      wjets.push_back(jet) ; 
      selwjets.push_back(ijet) ; 
    }
   } catch( ... )
	{
 	   cout << "Exception..." << run << " " << lumi << " " << event << endl; 
	}
 
  }



  //// Preselection HT
  if ( htak4.getHT() < HTMin_ ) return false; 


 // Pre-sel: Good primary vertices
 //
   if ( _isData ) {
     unsigned npv(0) ; 
     for ( unsigned ipv = 0; ipv < (h_vtxRho.product())->size(); ++ipv) {
        double vtxRho = (h_vtxRho.product())->at(ipv) ; 
        double vtxZ = (h_vtxZ.product())->at(ipv) ; 
        double vtxNdf = (h_vtxNdf.product())->at(ipv) ; 
        if ( abs(vtxRho) < 2. && abs(vtxZ) <= 24. && vtxNdf > 4 ) ++npv ; 
     }
     if ( npv < 1 ) return false ; 

   //h1_["cutflow"] -> AddBinContent(3, evtwt) ; 
 //  h1_["cutflow_nowt"] -> AddBinContent(3) ; 
 //
 // Pre-sel: MET filters: CSC beam halo and HBHE noise filters
     bool hbheNoiseFilter = h_hbheNoiseFilter.product() ; 
     if ( !hbheNoiseFilter ) return false ; 

     bool metfilterdecision(1) ; 
       for ( const string& metfilter : metFilters_ ) {
         vector<string>::const_iterator it ; 
         for (it = (h_metFiltersName.product())->begin(); it != (h_metFiltersName.product())->end(); ++it) {
           if ( it->find(metfilter) < std::string::npos) {
             metfilterdecision *= (h_metFiltersBit.product())->at( it - (h_metFiltersName.product())->begin() ) ; 
           }
         }
       }
     if ( !metfilterdecision ) return false ; 

   }

  if (goodAK4Jets.size() > 0) {
    std::sort(goodAK4Jets.begin(), goodAK4Jets.end(), Utilities::sortByCSV) ; 
    h1_["ak4highestcsv_nocuts"] -> Fill((goodAK4Jets.at(0)).getCSV()) ;
    std::sort(goodAK4Jets.begin(), goodAK4Jets.end(), Utilities::sortByPt<vlq::Jet>) ; 
  }



  //// Pick forwardmost AK4 jet
  double maxeta(0) ;
  vlq::Jet forwardestjet ; 
  for ( auto& jet : goodAK4Jets ) {
    if ( abs(jet.getEta()) > abs(maxeta) ) { 
      forwardestjet = jet ; 
      maxeta = jet.getEta() ; 
    }
  }

  double ptak8_1 = goodAK8Jets.at(0).getP4().Pt() ;
  double ptak8_2(0) ; 
  goodAK8Jets.size() > 1 ?  ptak8_2 = goodAK8Jets.at(1).getP4().Pt() : 0 ; 
  double mak8_1 = goodAK8Jets.at(0).getP4().Mag() ;
  double mak8_2(0) ; 
  goodAK8Jets.size() > 1 ?  mak8_2 = goodAK8Jets.at(1).getP4().Mag() : 0 ; 
  double mak8_12(0) ; 
  double detaLeading2AK8(-1) ; 
  if (goodAK8Jets.size() > 1) {
    TLorentzVector p4_ak8_12(goodAK8Jets.at(0).getP4() + goodAK8Jets.at(1).getP4()) ;
    mak8_12 = p4_ak8_12.Mag() ; 
    detaLeading2AK8 = abs(goodAK8Jets.at(0).getEta() - goodAK8Jets.at(1).getEta() ) ;
  }

  double ptak8leading ((goodAK8Jets.at(0)).getPt()) ; 
  double mak8leading ((goodAK8Jets.at(0)).getMass()) ; 
  double csvak8leading ((goodAK8Jets.at(0)).getCSV()) ; 


  h1_["ptak8leading"] -> Fill(ptak8leading,weight[0]) ; 
  h1_["etaak8leading"] -> Fill((goodAK8Jets.at(0)).getEta(),weight[0]) ;
  h1_["mak8leading"] -> Fill(mak8leading,weight[0]) ; 
  h1_["trimmedmak8leading"] -> Fill((goodAK8Jets.at(0)).getTrimmedMass(),weight[0]) ;
  h1_["softdropmak8leading"] -> Fill((goodAK8Jets.at(0)).getSoftDropMass(),weight[0]) ;
  h1_["csvak8leading"] -> Fill(csvak8leading,weight[0]) ;

  double ptak82nd (0) ;
  double mak82nd (0) ;
  double csvak82nd (0) ; 
  if (goodAK8Jets.size() > 1) {
    ptak82nd = (goodAK8Jets.at(1)).getPt() ; 
    mak82nd = (goodAK8Jets.at(1)).getMass() ; 
    csvak82nd = (goodAK8Jets.at(1)).getCSV() ; 
  }

  h1_["ptak82nd"] -> Fill(ptak82nd,weight[0]) ; 
  h1_["mak82nd"] -> Fill(mak82nd,weight[0]) ; 
  h1_["csvak82nd"] -> Fill(csvak82nd,weight[0]) ;

  h1_["ptak8leadingPlus2nd"] -> Fill(ptak8leading+ptak82nd,weight[0]) ; 

  std::sort(goodAK8Jets.begin(), goodAK8Jets.end(), Utilities::sortByMass<vlq::Jet>) ; 
  h1_["mak8highestm"] -> Fill((goodAK8Jets.at(0)).getMass(),weight[0]) ; 

  std::sort(goodAK8Jets.begin(), goodAK8Jets.end(), Utilities::sortByTrimmedMass) ; 
  h1_["trimmedmak8highesttrimmedm"] -> Fill((goodAK8Jets.at(0)).getTrimmedMass(),weight[0]) ; 

  std::sort(goodAK8Jets.begin(), goodAK8Jets.end(), Utilities::sortBySoftDropMass) ; 
  h1_["softdropmak8highestsoftdropm"] -> Fill((goodAK8Jets.at(0)).getSoftDropMass(),weight[0]) ; 

  std::sort(goodAK8Jets.begin(), goodAK8Jets.end(), Utilities::sortByPt<vlq::Jet>) ; 

  double ptak4leading ((h_jetAK4Pt.product())->at(ak4selIdxs.at(0))) ;   
  h1_["ptak4leading"] -> Fill(ptak4leading,weight[0]) ; 
  h1_["etaak4leading"] -> Fill((h_jetAK4Eta.product())->at(ak4selIdxs.at(0)),weight[0]) ; 

  double ptbjetleading (-1);
  double csvbjetleading (-1);
  double csvbjethighestcsv(-1); 
  if ( bjetIdxs.size() > 0 ) {
    ptbjetleading = (h_jetAK4Pt.product())->at(bjetIdxs.at(0)) ; 
    csvbjetleading = (h_jetAK4CSV.product())->at(bjetIdxs.at(0)) ; 
    h1_["ptbjetleading"] -> Fill(ptbjetleading,weight[0]) ; 
    h1_["etabjetleading"] -> Fill((h_jetAK4Eta.product())->at(bjetIdxs.at(0)),weight[0]) ; 
    h1_["csvbjetleading"] -> Fill(csvbjetleading,weight[0]) ; 

    std::sort(btaggedlooseAK4.begin(), btaggedlooseAK4.end(), Utilities::sortByCSV) ; 
    csvbjethighestcsv = (btaggedlooseAK4.at(0)).getCSV() ; 
    h1_["csvbjethighestcsv"] -> Fill(csvbjethighestcsv,weight[0]) ; 
    h1_["ptak4highestcsv"] -> Fill((btaggedlooseAK4.at(0)).getPt(),weight[0]) ;
    h1_["etaak4highestcsv"] -> Fill((btaggedlooseAK4.at(0)).getEta(),weight[0]) ;
  }

  h1_["ptak4forwardmost"] -> Fill(forwardestjet.getPt(),weight[0]) ; 
  h1_["etaak4forwardmost"] -> Fill(forwardestjet.getEta(),weight[0]) ; 

  h2_["pt_ak8_leading_2nd"] -> Fill(ptak8leading, ptak82nd,weight[0]) ; 
  h2_["m_ak8_leading_2nd"] -> Fill(mak8leading, mak82nd,weight[0]) ; 
  h2_["csv_ak8_leading_2nd"] -> Fill(csvak8leading, csvak82nd,weight[0]) ; 

  h1_["nak8_presel"] -> Fill(goodAK8Jets.size(),weight[0]) ; 
  h1_["nak4_presel"] -> Fill(goodAK4Jets.size(),weight[0]) ; 
  h1_["nbloose_presel"] -> Fill(btaggedlooseAK4.size(),weight[0]) ; 
  h1_["nbmedium_presel"] -> Fill(btaggedmediumAK4.size(),weight[0]) ; 
  h1_["nbtight_presel"] -> Fill(btaggedtightAK4.size(),weight[0]) ; 

  h1_["nwjet_presel"] -> Fill(wjets.size(),weight[0]) ; 
  h1_["nhjet_presel"] -> Fill(hjets.size(),weight[0]) ; 
  h1_["ntjet_presel"] -> Fill(tjets.size(),weight[0]) ; 

  h1_["ht_presel"] ->Fill(htak4.getHT(),weight[0]) ; 

  // Make H cands
  std::vector<vlq::ResolvedVjj> wcands, hcands ; 
  if (goodAK4Jets.size() > 1 && wjets.size() == 0) {
    double mmin (60), mmax(100), drmax(1.2), smdmin(0.0), smdmax(0.5) ;  
    VCandProducer WCandProducer(goodAK4Jets, mmin, mmax,drmax, smdmin, smdmax) ;  
    WCandProducer.getCands(wcands) ; 
  }
  if (btaggedmediumAK4.size() > 1 && hjets.size() == 0) {
    double mmin (100), mmax(140), drmax(1.2), smdmin(0.0), smdmax(0.5) ;  
    VCandProducer HCandProducer(btaggedmediumAK4, mmin, mmax,drmax, smdmin, smdmax) ;  
    HCandProducer.getCands(hcands) ; 
  }

  h1_["nwcand_presel"] -> Fill(hcands.size(),weight[0]) ; 
  h1_["nhcand_presel"] -> Fill(hcands.size(),weight[0]) ; 


  std::auto_ptr<unsigned> ngoodAK4Jets ( new unsigned(goodAK4Jets.size()) );
  std::auto_ptr<unsigned> ngoodAK8Jets ( new unsigned(goodAK8Jets.size()) );
  std::auto_ptr<unsigned> nTJets ( new unsigned(tjets.size()) );
  std::auto_ptr<unsigned> nHJets ( new unsigned(hjets.size()) );
  std::auto_ptr<unsigned> nWJets ( new unsigned(wjets.size()) );
  std::auto_ptr<unsigned> nbtaggedlooseAK4 ( new unsigned(btaggedlooseAK4.size()) );
  std::auto_ptr<unsigned> nbtaggedmediumAK4 ( new unsigned(btaggedmediumAK4.size()) );
  std::auto_ptr<unsigned> nbtaggedtightAK4 ( new unsigned(btaggedtightAK4.size()) );
  std::auto_ptr<double> htak4jets ( new double(htak4.getHT()) );
  std::auto_ptr<double> htak4trigjets ( new double(htak4trig.getHT()) );
  std::auto_ptr<double> htak8jets ( new double(htak8.getHT()) );
  std::auto_ptr<double> maxetaak4 ( new double(maxeta) );
  std::auto_ptr<double> MassLeading2AK8 ( new double(mak8_12) );
  std::auto_ptr<double> DeltaEtaLeading2AK8 ( new double(detaLeading2AK8) );
  std::auto_ptr<double> pt1stAK8   ( new double(ptak8_1) );
  std::auto_ptr<double> pt2ndAK8   ( new double(ptak8_2) );
  std::auto_ptr<double> mass1stAK8 ( new double(mak8_1) );
  std::auto_ptr<double> mass2ndAK8 ( new double(mak8_2) );
  std::auto_ptr<vector<unsigned> > ak4goodjets ( new vector<unsigned>(ak4selIdxs));
  std::auto_ptr<vector<unsigned> > ak8goodjets ( new vector<unsigned>(ak8selIdxs));
  std::auto_ptr<vector<unsigned> > bjetIdxsptr ( new vector<unsigned>(bjetIdxs));
  std::auto_ptr<vector<unsigned> > tjetIdxs ( new vector<unsigned>(seltjets));
  std::auto_ptr<vector<unsigned> > hjetIdxs ( new vector<unsigned>(selhjets));
  std::auto_ptr<vector<unsigned> > wjetIdxs ( new vector<unsigned>(selwjets));
  std::auto_ptr<vector<double> > EvtWeight ( new vector<double>(weight));
  std::auto_ptr<vector<string> > WeightName ( new vector<string>(weightnames));


  evt.put(ngoodAK4Jets, "ngoodAK4Jets") ; 
  evt.put(ngoodAK8Jets, "ngoodAK8Jets") ; 
  evt.put(nTJets, "nTJets") ; 
  evt.put(nHJets, "nHJets") ; 
  evt.put(nWJets, "nWJets") ; 
  evt.put(nbtaggedlooseAK4, "nbtaggedlooseAK4") ; 
  evt.put(nbtaggedmediumAK4, "nbtaggedmediumAK4") ; 
  evt.put(nbtaggedtightAK4, "nbtaggedtightAK4") ; 
  evt.put(maxetaak4, "maxetaak4") ; 
  evt.put(MassLeading2AK8, "MassLeading2AK8") ; 
  evt.put(DeltaEtaLeading2AK8, "DeltaEtaLeading2AK8") ; 
  evt.put(pt1stAK8  , "pt1stAK8") ; 
  evt.put(pt2ndAK8  , "pt2ndAK8") ; 
  evt.put(mass1stAK8, "mass1stAK8") ; 
  evt.put(mass2ndAK8, "mass2ndAK8") ; 
  evt.put(htak4jets, "htak4jets") ; 
  evt.put(htak8jets, "htak8jets") ; 
  evt.put(htak4trigjets, "htak4trigjets");
  evt.put(ak4goodjets, "ak4goodjets");
  evt.put(ak8goodjets, "ak8goodjets");
  evt.put(bjetIdxsptr, "bjetIdxs");
  evt.put(tjetIdxs, "tjetIdxs");
  evt.put(hjetIdxs, "hjetIdxs");
  evt.put(wjetIdxs, "wjetIdxs");
  evt.put(EvtWeight, "EvtWeight");
  evt.put(WeightName, "WeightName");

  return true ; 

}

// ------------ method called once each job just before starting event loop  ------------
void Presel::beginJob() {

  weightFile = new TFile(weightFileName_.c_str(),"READ");


  for (unsigned i = 0 ; i < weights_.size() ; i++) {
      wGraph[i] = (TGraphAsymmErrors * ) weightFile->Get(weights_[i].c_str());
      if (wGraph[i] == 0) cout << "Can't find hist: " << weights_[i].c_str() <<endl; 
      weightnames.push_back(weights_[i]);      
//	cout << "Added " << weights_[i]<<endl;
  }

  h1_["ptak8leading"]  = fs->make<TH1D>("ptak8leading"  ,";p_T(leading AK8 jet) [GeV];;"      , 60, 0., 3000.) ; 
  h1_["ptak4leading"]  = fs->make<TH1D>("ptak4leading"  ,";p_T(leading AK4 jet) [GeV];;"      , 60, 0., 3000.) ; 
  h1_["ptbjetleading"]  = fs->make<TH1D>("ptbjetleading"  ,";p_T(leading b jet) [GeV];;"      , 60, 0., 3000.) ; 

  h1_["etaak8leading"] = fs->make<TH1D>("etaak8leading", ";#eta(leading AK8 jet);;" , 100 ,-5. ,5.) ; 
  h1_["etaak4leading"] = fs->make<TH1D>("etaak4leading", ";#eta(leading AK4 jet);;" , 100 ,-5. ,5.) ; 
  h1_["etabjetleading"] = fs->make<TH1D>("etabjetleading", ";#eta(leading b jet);;" , 100 ,-5. ,5.) ; 

  h1_["ptak82nd"]  = fs->make<TH1D>("ptak82nd"  ,";p_T(2nd AK8 jet) [GeV];;"      , 60, 0., 3000.) ; 
  h1_["ptak8leadingPlus2nd"]  = fs->make<TH1D>("ptak8leadingPlus2nd"  ,";p_T(leading AK8 jet)+p_T (2nd AK8 jet) [GeV];;"      , 60, 0., 3000.) ; 

  h1_["csvbjetleading"] = fs->make<TH1D>("csvbjetleading", ";CSV (leading b jet);;" ,50 ,0. ,1.) ; 
  h1_["csvbjethighestcsv"] = fs->make<TH1D>("csvbjethighestcsv", ";max. CSV b jet;;" ,50 ,0. ,1.) ; 

  h1_["ak4highestcsv_nocuts"] = fs->make<TH1D>("ak4highestcsv_nocuts", ";max. CSV of AK4 jets;;" , 50, 0., 1.) ; 

  h1_["ptak4highestcsv"] = fs->make<TH1D>("ptak4highestcsv", ";p_T(highest CSV AK4 jet);;" , 150, 0., 3000.) ; 
  h1_["etaak4highestcsv"] = fs->make<TH1D>("etaak4highestcsv", ";#eta(highest CSV AK4 jet);;" , 100 ,-5. ,5.) ; 

  h1_["ptak4forwardmost"] = fs->make<TH1D>("ptak4forwardmost", ";p_T(forwardmost AK4 jet);;" , 150, 0., 3000.) ; 
  h1_["etaak4forwardmost"] = fs->make<TH1D>("etaak4forwardmost", ";#eta(forwardmost AK4 jet);;" , 100 ,-5. ,5.) ; 

  h1_["mak8leading"] = fs->make<TH1D>("mak8leading", ";M(leading AK8 jet) [GeV];;" ,500 ,0., 1000.) ; 
  h1_["mak8highestm"] = fs->make<TH1D>("mak8highestm", ";M(most massive AK8 jet) [GeV];;" ,500 ,0., 1000.) ; 

  h1_["prunedmak8leading"] = fs->make<TH1D>("prunedmak8leading", ";M(pruned leading AK8 jet) [GeV];;" ,500 ,0., 1000.) ; 

  h1_["trimmedmak8leading"] = fs->make<TH1D>("trimmedmak8leading", ";M(trimmed leading AK8 jet) [GeV];;" ,500 ,0., 1000.) ; 
  h1_["trimmedmak8highesttrimmedm"] = fs->make<TH1D>("trimmedmak8highesttrimmedm", ";M(highest trimmed mass AK8 jet) [GeV];;" ,500 ,0., 1000.) ; 

  h1_["softdropmak8leading"] = fs->make<TH1D>("softdropmak8leading", ";M(leading AK8 jet) [GeV];;" ,500 ,0., 1000.) ; 
  h1_["softdropmak8highestsoftdropm"] = fs->make<TH1D>("softdropmak8highestsoftdropm", ";M(highest soft drop mass AK8 jet) [GeV];;" ,500 ,0., 1000.) ; 

  h1_["mak82nd"] = fs->make<TH1D>("mak82nd", ";M(2nd AK8 jet) [GeV];;" ,500 ,0., 1000.) ; 

  h1_["csvak8leading"] = fs->make<TH1D>("csvak8leading", ";CSV (leading AK8 jet);;" ,50 ,0. ,1.) ;
  h1_["csvak82nd"] = fs->make<TH1D>("csvak82nd", ";CSV (2nd AK8 jet);;" ,50 ,0. ,1.) ;

  h2_["pt_ak8_leading_2nd"] = fs->make<TH2D>("pt_ak8_leading_2nd", ";p_T (leading AK8 jet) [GeV];p_T (2nd AK8 jet) [GeV];" ,60, 0., 3000. ,60, 0., 3000.) ; 

  h2_["m_ak8_leading_2nd"] = fs->make<TH2D>("m_ak8_leading_2nd", ";M(leading AK8 jet) [GeV];M(2nd AK8 jet) [GeV];" ,500, 0., 1000. ,500, 0., 1000.) ; 

  h2_["csv_ak8_leading_2nd"] = fs->make<TH2D>("csv_ak8_leading_2nd", ";CSV(leading AK8 jet) ;CSV(2nd AK8 jet) ;" ,50 ,0. ,1. ,50 ,0. ,1.) ;  

  h1_["nak8_nocuts"] = fs->make<TH1D>("nak8_nocuts", ";AK8 jet multiplicity before cuts;;" , 11, -0.5,10.5) ; 
  h1_["nak4_nocuts"] = fs->make<TH1D>("nak4_nocuts", ";AK4 jet multiplicity before cuts;;" , 21, -0.5,20.5) ; 
  h1_["nbloose_nocuts"] = fs->make<TH1D>("nbloose_nocuts", ";b jet (loose) multiplicity before cuts;;" , 11, -0.5,10.5) ; 
  h1_["nbmedium_nocuts"] = fs->make<TH1D>("nbmedium_nocuts", ";b jet (medium) multiplicity before cuts;;" , 11, -0.5,10.5) ; 
  h1_["nbtight_nocuts"] = fs->make<TH1D>("nbtight_nocuts", ";b jet (tight) multiplicity before cuts;;" , 11, -0.5,10.5) ; 

  h1_["ak4_pt_nocuts"] = fs->make<TH1D>("ak4_pt_nocuts", "p_T (all AK4 Jets - no cuts);;" , 150, 0 , 3000);
  h1_["ak4_eta_nocuts"] = fs->make<TH1D>("ak4_eta_nocuts", "#eta (all AK4 Jets - no cuts);;" , 100, -5, 5);

  h1_["htak4_nocuts"] = fs->make<TH1D>("htak4_nocuts", " H_T (AK4 Jets);;", 300 ,0 ,6000);
  h1_["htak4trig_nocuts"] = fs->make<TH1D>("htak4trig_nocuts", " H_T (AK4 Jets - Trig. Def);;", 300 ,0 ,6000);

  h1_["nak8_presel"] = fs->make<TH1D>("nak8_presel", ";AK8 jet multiplicity;;" , 11, -0.5,10.5) ; 
  h1_["nak4_presel"] = fs->make<TH1D>("nak4_presel", ";AK4 jet multiplicity;;" , 21, -0.5,20.5) ; 
  h1_["nbloose_presel"] = fs->make<TH1D>("nbloose_presel", ";b jet multiplicity (Loose OP);;" , 11, -0.5,10.5) ; 
  h1_["nbmedium_presel"] = fs->make<TH1D>("nbmedium_presel", ";b jet multiplicity (Medium OP);;" , 11, -0.5,10.5) ; 
  h1_["nbtight_presel"] = fs->make<TH1D>("nbtight_presel", ";b jet multiplicity (Tight OP);;" , 11, -0.5,10.5) ; 

  h1_["nwjet_presel"] = fs->make<TH1D>("nwjet_presel", ";W jet multiplicity;;" , 11, -0.5,10.5) ; 
  h1_["nhjet_presel"] = fs->make<TH1D>("nhjet_presel", ";H jet multiplicity;;" , 11, -0.5,10.5) ; 
  h1_["ntjet_presel"] = fs->make<TH1D>("ntjet_presel", ";top jet multiplicity;;" , 11, -0.5,10.5) ; 

  h1_["nwcand_presel"] = fs->make<TH1D>("nwcand_presel", ";W candidate multiplicity;;" , 11, -0.5,10.5) ; 
  h1_["nhcand_presel"] = fs->make<TH1D>("nhcand_presel", ";H candidate multiplicity;;" , 11, -0.5,10.5) ; 

  h1_["ht_presel"] = fs->make<TH1D>("ht_presel" ,";H_T (AK4 jets) [GeV]", 300, 0., 6000.) ; 

  h1_["nbloose"] = fs->make<TH1D>("nbloose", ";b jet multiplicity (Loose OP);;" , 11, -0.5,10.5) ; 
  h1_["nbmedium"] = fs->make<TH1D>("nbmedium", ";b jet multiplicity (Medium OP);;" , 11, -0.5,10.5) ; 

  h1_["nwjet"] = fs->make<TH1D>("nwjet", ";W jet multiplicity;;" , 11, -0.5,10.5) ; 
  h1_["nhjet"] = fs->make<TH1D>("nhjet", ";H jet multiplicity;;" , 11, -0.5,10.5) ; 
  h1_["ntjet"] = fs->make<TH1D>("ntjet", ";top jet multiplicity;;" , 11, -0.5,10.5) ; 

  h1_["nwcand"] = fs->make<TH1D>("nwcand", ";W candidate multiplicity;;" , 11, -0.5,10.5) ; 
  h1_["nhcand"] = fs->make<TH1D>("nhcand", ";H candidate multiplicity;;" , 11, -0.5,10.5) ; 

//  h1_["ht"] = fs->make<TH1D>("ht" ,";H_T (AK4 jets) [GeV]", 200, 0., 4000.) ; 

  h1_["ptleadinghcand"] = fs->make<TH1D>("ptleadinghcand" ,";p_T (leading H cands) [GeV]", 150, 0., 3000.) ; 
  h1_["ptleadingwcand"] = fs->make<TH1D>("ptleadingwcand" ,";p_T (leading W cands) [GeV]", 150, 0., 3000.) ; 

  h1_["ptleadinghjet"] = fs->make<TH1D>("ptleadinghjet" ,";p_T (leading H-tagged jets) [GeV]", 150, 0., 3000.) ; 
  h1_["ptleadingwjet"] = fs->make<TH1D>("ptleadingwjet" ,";p_T (leading W-tagged jets) [GeV]", 150, 0., 3000.) ; 
  h1_["ptleadingtjet"] = fs->make<TH1D>("ptleadingtjet" ,";p_T (leading top-tagged jets) [GeV]", 150, 0., 3000.) ; 
/*
  h1_["ht_nwjetGt0"] = fs->make<TH1D>("ht_nwjetGt0" ,";H_T (AK4 jets) [GeV]", 200, 0., 4000.) ; 
  h1_["ht_nhjetGt0"] = fs->make<TH1D>("ht_nhjetGt0" ,";H_T (AK4 jets) [GeV]", 200, 0., 4000.) ; 
  h1_["ht_ntjetGt0"] = fs->make<TH1D>("ht_ntjetGt0" ,";H_T (AK4 jets) [GeV]", 200, 0., 4000.) ; 
  h1_["ht_nwhjetGt0"] = fs->make<TH1D>("ht_nwhjetGt0" ,";H_T (AK4 jets) [GeV]", 200, 0., 4000.) ; 
  h1_["ht_nhtjetGt0"] = fs->make<TH1D>("ht_nhtjetGt0" ,";H_T (AK4 jets) [GeV]", 200, 0., 4000.) ; 
  h1_["ht_nwcandGt0"] = fs->make<TH1D>("ht_nwcandGt0" ,";H_T (AK4 cands) [GeV]", 200, 0., 4000.) ; 
  h1_["ht_nhcandGt0"] = fs->make<TH1D>("ht_nhcandGt0" ,";H_T (AK4 cands) [GeV]", 200, 0., 4000.) ; 
*/
  h1_["ptak8"]  = fs->make<TH1D>("ptak8"  ,";p_T(AK8 jet) [GeV]"         , 150, 0., 3000.) ; 
  h1_["pthjets"] = fs->make<TH1D>("pthjets" ,";p_T (H-tagged jets) [GeV]", 150, 0., 3000.) ; 
  h1_["csvhjets"] = fs->make<TH1D>("csvhjets", ";CSV (H-tagged jets);;" ,50 ,0. ,1.) ;

  h1_["drwh"] = fs->make<TH1D>("drwh", ";#DeltaR(H,W);;", 50, 0, 5.) ;  


}

// ------------ method called once each job just after ending the event loop  ------------
void 
Presel::endJob() {

  weightFile->Close();

}

// ------------ method called when starting to processes a run  ------------
/*
void
Presel::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
Presel::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
Presel::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
Presel::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Presel::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(Presel);
