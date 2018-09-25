/**
 * @file   SSNetTest_module.cc
 * @brief  Takes an SSNetHits.root file 
 * @author ksutton
 * 
*/


// LArSoft includes
 #include "lardataobj/Simulation/SimChannel.h"
 #include "larsim/Simulation/LArG4Parameters.h"
 #include "lardataobj/RecoBase/Hit.h"
 #include "lardataobj/RecoBase/Cluster.h"
 #include "larcore/Geometry/Geometry.h"
 #include "larcore/Geometry/GeometryCore.h"
 #include "nusimdata/SimulationBase/MCParticle.h"
 #include "nusimdata/SimulationBase/MCTruth.h"
 #include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
 #include "lardataobj/RecoBase/Vertex.h"
 #include "lardataobj/RecoBase/PFParticle.h"
 #include "lardataobj/RecoBase/Track.h"
 #include "lardataobj/RecoBase/Shower.h"
 #include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
 
// Framework includes
 #include "canvas/Utilities/Exception.h"
 #include "art/Framework/Core/EDAnalyzer.h"
 #include "art/Framework/Principal/Event.h"
 #include "art/Framework/Principal/Handle.h"
 #include "art/Framework/Services/Registry/ServiceHandle.h"
 #include "art/Framework/Services/Optional/TFileService.h"
 #include "art/Framework/Core/ModuleMacros.h"
 #include "canvas/Persistency/Common/FindManyP.h"

 // utility libraries
 #include "messagefacility/MessageLogger/MessageLogger.h"
 #include "fhiclcpp/ParameterSet.h"
 #include "fhiclcpp/types/Table.h"
 #include "fhiclcpp/types/Atom.h"
 #include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT includes. Note: To look up the properties of the ROOT classes,
 // use the ROOT web site; e.g.,// <http://root.cern.ch/root/html532/ClassIndex.html>
 #include "TH1.h"
 #include "TH2.h"
 #include "TTree.h"
 #include "TLorentzVector.h"
 #include "TVector3.h"

 // C++ Includes
 #include <map>
 #include <vector>
 #include <string>
 #include <cmath>
 #include <memory>

namespace{
double DetectorDiagonal(geo::GeometryCore const& geom);
bool matches(const recob::Hit* Hit, const recob::Hit* ssHit);
double distToVtx(double fTimeToCMConstant, double fWireToCMConstant, double radius, double vertex_time, double vertex_wire, const recob::Hit*  hit);
bool inROI(double radius, double dist);
double calcWire(double Y, double Z, int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo );
TVector3 findEnd(TVector3* shower_start, double* shower_length, TVector3* shower_dir);
double calcTime(double X,int plane,int fTPC,int fCryostat, detinfo::DetectorProperties const& detprop);
TVector2 calcNormVec(TVector2 shower_start_plane,TVector2 shower_end_plane);
double angleFromShower(TVector2 shower_direction_plane, const recob::Hit*  hit, double vertex_time, double vertex_wire);
//double getRadius(double rad_in_cm);
}

namespace lar {
namespace example {
  class SSNetTest : public art::EDAnalyzer
  {
  public:
	struct Config {
      
      // save some typing:
      	using Name = fhicl::Name;
     	using Comment = fhicl::Comment;
      // one Atom for each parameter
        fhicl::Atom<art::InputTag> SimulationLabel {Name("SimulationLabel"),
      		Comment("tag of the input data product with the detector simulation information")
        };
      
      fhicl::Atom<int> PxType {
      Name("PxType"),
      Comment("")
      };
                                                                 
     fhicl::Atom<art::InputTag> SSHitLabel {
      Name("SSHitLabel"),
      Comment("tag of the input SSnet data product with reconstructed hits")
      };
                                          
     fhicl::Atom<art::InputTag> HitLabel {
      Name("HitLabel"),
      Comment("tag of the input data product with reconstructed hits")
      };
   
      
      fhicl::Atom<art::InputTag> ClusterLabel {
      Name("ClusterLabel"),
      Comment("tag of the input data product with reconstructed clusters")
      };
                                                                                                                                          
      fhicl::Atom<int> PDGcode {
      Name("PDGcode"),
      Comment("particle type (PDG ID) of the primary particle to be selected")
      };     //                                                                                                                                                                         
     fhicl::Atom<double> BinSize {
	Name("BinSize"),                                                                                                                                                                    Comment("dx [cm] used for the dE/dx calculation")                                                                                                                                                                                       };
     fhicl::Atom<art::InputTag> VertexLabel {
	Name("Vertex"),
	Comment("currently using Pandora vertex")                                                                                                                                                                                       };
     fhicl::Atom<art::InputTag> PFPLabel {
	Name("PFP"),
	Comment("tag of the PFParticle producer")                                                                                                                                                                                       };
     fhicl::Atom<art::InputTag> TrackLabel {
	Name("TrackLabel"),
	Comment("tag of the track producer")                                                                                                                                                                                       };
     fhicl::Atom<art::InputTag> ShowerLabel {
	Name("ShowerLabel"),
	Comment("tag of the shower producer")                                                                                                                                                                                       };
    fhicl::Atom<int> Planes {
	Name("Planes"),
	Comment("number of planes to be considered (should be 3)")
 };    

   fhicl::Atom<double> Radius {
	Name("Radius"),
	Comment("radius in cm of the ROI from the vertex")
 }; 
   
   fhicl::Atom<double> WireToCMConstant {
	Name("WireToCMConstant"),
	Comment("conversion factor for wire # to cm")
 }; 
  
   fhicl::Atom<double> TimeToCMConstant {
	Name("TimeToCMConstant"),
	Comment("conversion factor for time # to cm")
 };  
};//struct
    
     using Parameters = art::EDAnalyzer::Table<Config>;

    explicit SSNetTest(Parameters const& config); 
  
    virtual void beginJob() override;
    virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;

  private:

    // The parameters we'll read from the .fcl file.
    art::InputTag fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector
    art::InputTag fHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fSSHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fClusterProducerLabel;    ///< The name of the producer that created clusters
    int fSelectedPDG;                     ///< PDG code of particle we'll focus on
    double fBinSize;                      ///< For dE/dx work: the value of dx. 
    int fPxType;
    art::InputTag fVertexProducerLabel; ///<name of the producer that created the vertices
    art::InputTag fPFPLabel; ///name of PFParticle producer
    art::InputTag fTrackProducerLabel; ///name of track producer
    art::InputTag fShowerProducerLabel; ///name of track producer
    int fPlanes; ///the number of planes
    double fRadius; //radius of ROI in cm
    double fWireToCMConstant; //converstion factor wire to cm
    double fTimeToCMConstant; //converstion factor time to cm

// vector of shower-like hit indices
//   std::vector<size_t> _shrhits;

 // Pointers to the histograms we'll create. 
    TH1D* fPDGCodeHist;     ///< PDG code of all particles
    TH1D* fMomentumHist;    ///< momentum [GeV] of all selected particles
    TH1D* fTrackLengthHist; ///< true length [cm] of all selected particles

    // The n-tuples we'll create.
    TTree* fmytree;     ///< tuple with information about all events
    TTree* fselectTree; ///< contains information about 1 track 1 shower events
    TTree* fROITree; ///<contains information about sshits inside the ROI

    // The variables that will go into the n-tuple.
    int fEvent;     ///< number of the event being processed
    int fRun;       ///< number of the run being processed
    int fSubRun;    ///< number of the sub-run being processed
    int fSimPDG;       ///< PDG ID of the particle begin processed
    int fSimTrackID;   ///< GEANT ID of the particle begin processed
    

    int fTPC;        ///the TPC ID
    int fCryostat;   ///the cryostat ID
    int total_num_vertices; ///equivalent to number of ROI's in file
    
    Float_t _xpos, _ypos, _zpos; // xyz of vertex

    double fStartXYZT[4]; ///< (x,y,z,t) of the true start of the particle
    double fEndXYZT[4];   ///< (x,y,z,t) of the true end of the particle
    double fStartPE[4];   ///< (Px,Py,Pz,E) at the true start of the particle
    double fEndPE[4];     ///< (Px,Py,Pz,E) at the true end of the particle

    int fnum_total_hits; //the total number of gaushits 
    int fnum_total_sshits; //the total number of sshits
    int fnum_total_matched_hits; //the total number of matched sshits
    double fratio_total_sshits_hits; //the total number of sshits over total number of hits
    
    int fnum_sshits_ROI_no_shower; //total number of sshits within each ROI once the shower is removed
    
    int fnum_sshits_shower; //number of sshits in a shower
    int fnum_hits_shower; //number of hits in a shower
    double fratio_hits_shower; //ratio of sshits to hits in a shower

    int fnum_sshits_track; //number of sshits in a track
    int fnum_hits_track; //number of hits in a track
    double fratio_hits_track; //ratio of sshits to hits in a track

    double fradial_dist_sshit_vtx; //the distance of each sshit in ROI to the vertex (shower removed)
    double fopening_angle_shower_sshit; //the angle on each plane between the shower direction and the vtx-sshit
	 
    //int fnum_sshits_ROI_with_shower; 
 
    geo::GeometryCore const* fGeometry; 

    //stuff to do the vertex to plane mapping
   // art::ServiceHandle<geo::Geometry>  geo;
  //  geo::Geometry const* fGeo;
    detinfo::DetectorProperties const* fDetprop;    

 }; // class SSNetTest

SSNetTest::SSNetTest(Parameters const& config) // Initialize member data here.
	: EDAnalyzer(config)
 	, fSimulationProducerLabel(config().SimulationLabel())
    	, fHitProducerLabel       (config().HitLabel())
    	, fSSHitProducerLabel       (config().SSHitLabel())
    	, fClusterProducerLabel   (config().ClusterLabel())
    	, fSelectedPDG            (config().PDGcode())
    	, fBinSize                (config().BinSize())
//	, fVertexProducerLabel (config().Vertex())
       // , fLArCVLocation  	  (config().LArCVLocation())
        , fPxType          	  (config().PxType())
        , fVertexProducerLabel    (config().VertexLabel()) 
	, fPFPLabel		  (config().PFPLabel())  
	, fTrackProducerLabel		  (config().TrackLabel())  	
	, fShowerProducerLabel		  (config().ShowerLabel())  	
	, fPlanes			  (config().Planes())
	, fRadius			  (config().Radius())
	, fWireToCMConstant		  (config().WireToCMConstant())
	, fTimeToCMConstant		  (config().TimeToCMConstant())
	// fInHitProducer   = p.get<std::string>("InHitProducer","gaushit");
        //fPxThresholdHigh = p.get<double>     ("PxThresholdHigh"        );
        //fPxThresholdLow  = p.get<double>     ("PxThresholdLow"         );
       // produces< std::vector< recob::Hit > >(Form("shrhit%zu",(size_t)(fPxThresholdHigh*100.)));
        //produces< std::vector< recob::Hit > >(Form("shrhit%zu",(size_t)(fPxThresholdLow*100.) ));

  {
    // get a pointer to the geometry service provider
     fGeometry = lar::providerFrom<geo::Geometry>();
     //fGeo =  lar::providerFrom<geo::Geometry>();

    //get a pointer to the detector properties service provider
     fDetprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
 }


  
  //-----------------------------------------------------------------------
  void SSNetTest::beginJob()
  {
    // Get the detector length, to determine the maximum bin edge of one
    // of the histograms.
    const double detectorLength = DetectorDiagonal(*fGeometry);
  
   //get the tpc and cryostat ID's  
    auto const TPC = (*fGeometry).begin_TPC();
    auto ID = TPC.ID();
    fCryostat = ID.Cryostat;
    fTPC = ID.TPC;
    std::cout<<TPC.ID()<<"= the beginning TPC ID" <<std::endl;
    std::cout<<"the cryostat id = "<<fCryostat<<std::endl;  
    std::cout<<"the tpc id = "<<fTPC<<std::endl;  

    total_num_vertices = 0;
//check that there's only one definition
/* geo::GeometryCore::TPC_iterator iTPC = (*fGeometry).begin_TPC();  
 auto const end_TPC = (*fGeometry).end_TPC();
    //std::cout<<end_TPC.ID()<<"= the ending TPC ID" <<std::endl;
    while (iTPC != end_TPC) {
  	std::cout << "current TPC: " << iTPC.ID() << std::endl;
  // the TPC descriptor object
 // 	geo::TPCGeo const& TPC = *iTPC;
       ++iTPC;
         } // while
*/

    
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;
  
    fPDGCodeHist     = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
    fMomentumHist    = tfs->make<TH1D>("mom",     ";particle Momentum (GeV);",    100, 0.,    10.);
    fTrackLengthHist = tfs->make<TH1D>("length",  ";particle track length (cm);", 200, 0, detectorLength);
 
    // Define our n-tuples, which are limited forms of ROOT
    // TTrees. Start with the TTree itself.
    fmytree     = tfs->make<TTree>("SSNetTestSimulation",    "SSNetTestSimulation");
    fselectTree = tfs->make<TTree>("SSNetTestSelectTopologies",    "SSNetTestSelectTopologies");
    fROITree = tfs->make<TTree>("SSNetTestROI",    "SSNetTestROI");

    // Define the branches (columns) of our simulation n-tuple. When
    // we write a variable, we give the address of the variable to
    // TTree::Branch.
    fmytree->Branch("Event",       &fEvent,          "Event/I");
    fmytree->Branch("SubRun",      &fSubRun,         "SubRun/I");
    fmytree->Branch("Run",         &fRun,            "Run/I");
    fmytree->Branch("TrackID",     &fSimTrackID,     "TrackID/I");
    fmytree->Branch("PDG",         &fSimPDG,         "PDG/I");
    // When we write arrays, we give the address of the array to
    // TTree::Branch; in C++ this is simply the array name.
    fmytree->Branch("StartXYZT",   fStartXYZT,       "StartXYZT[4]/D");
    fmytree->Branch("EndXYZT",     fEndXYZT,         "EndXYZT[4]/D");
    fmytree->Branch("StartPE",     fStartPE,         "StartPE[4]/D");
    fmytree->Branch("EndPE",       fEndPE,           "EndPE[4]/D");
    // For a variable-length array: include the number of bins.
    //fmytree->Branch("NdEdx",       &fSimNdEdxBins,   "NdEdx/I");
    // ROOT can understand fairly well vectors of numbers (and little more)
    //fmytree->Branch("dEdx",        &fSimdEdxBins);
    fmytree->Branch("X_pos_vertex", &_xpos,   "_xpos/F");
    fmytree->Branch("Y_pos_vertex", &_ypos,   "_ypos/F");
    fmytree->Branch("Z_pos_vertex", &_zpos,   "_zpos/F");
    fmytree->Branch("total_hits", &fnum_total_hits,   "total_hits/I");
    fmytree->Branch("total_sshits", &fnum_total_sshits,   "total_sshits/I");
    fmytree->Branch("total_matched_hits", &fnum_total_matched_hits,   "total_matched_hits/I");
    fmytree->Branch("ratio_total_sshits_hits", &fratio_total_sshits_hits,   "ratio_total_sshits_hits/D");
   
/*
     fmytree->Branch("num_sshits_ROI_no_shower", &fnum_sshits_ROI_no_shower,   "num_sshits_ROI_no_shower/I");
    fmytree->Branch("num_sshits_shower", &fnum_sshits_shower,   "num_sshits_shower/I");
    fmytree->Branch("num_hits_shower", &fnum_hits_shower,   "num_hits_shower/I");
    fmytree->Branch("ratio_hits_shower", &fratio_hits_shower,   "ratio_hits_shower/D");
    
    fmytree->Branch("num_sshits_track", &fnum_sshits_track,   "num_sshits_track/I");
    fmytree->Branch("num_hits_track", &fnum_hits_track,   "num_hits_track/I");
    fmytree->Branch("ratio_hits_track", &fratio_hits_track,   "ratio_hits_track/D");
*/
    fROITree->Branch("num_sshits_ROI_no_shower", &fnum_sshits_ROI_no_shower,   "num_sshits_ROI_no_shower/I");
    fselectTree->Branch("num_sshits_shower", &fnum_sshits_shower,   "num_sshits_shower/I");
    fselectTree->Branch("num_hits_shower", &fnum_hits_shower,   "num_hits_shower/I");
    fselectTree->Branch("ratio_hits_shower", &fratio_hits_shower,   "ratio_hits_shower/D");
    
    fselectTree->Branch("num_sshits_track", &fnum_sshits_track,   "num_sshits_track/I");
    fselectTree->Branch("num_hits_track", &fnum_hits_track,   "num_hits_track/I");
    fselectTree->Branch("ratio_hits_track", &fratio_hits_track,   "ratio_hits_track/D");

    fROITree->Branch("radial_dist_sshit_vtx", &fradial_dist_sshit_vtx, "radial_dist_sshit_vtx/D");
    fROITree->Branch("fopening_angle_shower_sshit", &fopening_angle_shower_sshit, "opening_angle_shower_sshit/D");


}   
  //-----------------------------------------------------------------------
  void SSNetTest::beginRun(const art::Run& /*run*/)
  {
 //   std::cout<<"flag 1"<<std::endl;
    art::ServiceHandle<sim::LArG4Parameters> larParameters;
   // fElectronsToGeV = 1./larParameters->GeVToElectrons();
  //std::cout<<"flag 2"<<std::endl;
}

  //-----------------------------------------------------------------------
  void SSNetTest::analyze(const art::Event& event) 
  {
    std::cout<<"------- starting event --------"<<std::endl;
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

   //read in PFParticles
   art::ValidHandle< std::vector<recob::PFParticle> > PFPHandle
   = event.getValidHandle<std::vector<recob::PFParticle>>
        (fPFPLabel);

  //read in clusters
  art::ValidHandle< std::vector<recob::Cluster> > clusterHandle
      = event.getValidHandle<std::vector<recob::Cluster>>
        (fClusterProducerLabel);

  //read in tracks
  art::ValidHandle< std::vector<recob::Track> > trackHandle
      = event.getValidHandle<std::vector<recob::Track>>
        (fTrackProducerLabel);

 //read in showers
  art::ValidHandle< std::vector<recob::Shower> > showerHandle
      = event.getValidHandle<std::vector<recob::Shower>>
        (fShowerProducerLabel);
/*
 //read in vertices
  art::ValidHandle< std::vector<recob::Vertex> > vertexHandle
      = event.getValidHandle<std::vector<recob::Vertex>>
        (fVertexProducerLabel);
*/
 // grab vertices, tracks, and showers associated to PFParticle
  const art::FindManyP<recob::Vertex> pfp_vtx_assn_v(PFPHandle, event, fPFPLabel);  
  const art::FindManyP<recob::Track  > pfp_trk_assn_v(PFPHandle, event, fPFPLabel);
  const art::FindManyP<recob::Shower> pfp_shr_assn_v(PFPHandle, event, fPFPLabel);
/* 
// grab showers  associated to vertices
   art::FindManyP<recob::Shower> vtx_shr_assn_v(vertexHandle, event, fVertexProducerLabel);
*/ 
// grab hits associated to track
  art::FindManyP<recob::Hit> trk_hit_assn_v(trackHandle, event, fTrackProducerLabel);
 
  // grab hits associated to shower
  art::FindManyP<recob::Hit> shr_hit_assn_v(showerHandle, event, fShowerProducerLabel);


/*
 *Perform hit-matching between hits and sshits
 */

//load hits
 art::Handle< std::vector<recob::Hit> > hitHandle;
 if (!event.getByLabel(fHitProducerLabel, hitHandle)) return;
 fnum_total_hits = hitHandle->size();

//load sshits  
 art::Handle< std::vector<recob::Hit> > sshitHandle;
 if (!event.getByLabel(fSSHitProducerLabel, sshitHandle)) return;

//list storing pairs of associated hits and sshits
std::list<std::pair<const recob::Hit*, const recob::Hit* >> _hitlist; //a list of hits
 
int n = 0; //the total number of sshits
fnum_total_sshits = sshitHandle->size();

// For every ssHit:
    for ( auto const& sshit : (*sshitHandle) )
      {
	n++;
	//for every Hit
	for ( auto const& hit : (*hitHandle) ){

		//if the wire, plane, and time match
		if (matches(&hit, &sshit)== true){	
			//make a pair of pointers to hits and add them to the list
			_hitlist.emplace(_hitlist.begin(), &hit, &sshit);
		}//if they match	
	}//for each Hit
 } // for each SSHit

fnum_total_matched_hits = _hitlist.size();
fratio_total_sshits_hits = (double)fnum_total_sshits/fnum_total_hits;
std::cout<<"number of hits = "<<fnum_total_hits<<", number of sshits = "<<n<<"/"<<fnum_total_sshits<<", number of matches = "<<fnum_total_matched_hits<<std::endl;
fmytree->Fill();
 
/*
 * loop over PFParticles and find 1 shower 1 track topologies
 *
 */

  std::vector<art::Ptr<recob::Track>> my_trks; //vector of pointers to tracks in case of 1 shower 1 track events
  std::vector<art::Ptr<recob::Shower>> my_shrs; //vector of pointers to tracks in case of 1 shower 1 track events
  std::vector<double*> my_vtxs; //vector of pointers to the first element in a double for the vertex XYZ position (3D)
 
 // std::cout<<"number of PFParticles: "<<PFPHandle->size()<<std::endl;
  //for each PFP
  for ( size_t pfp_index = 0; pfp_index != PFPHandle->size(); ++pfp_index ){
	auto const& pfp = PFPHandle->at(pfp_index);
	int pdg = pfp.PdgCode();
	
	//only want primary particles
	if (pfp.IsPrimary() == false){continue;}	
	
	//check for proton
	if (pdg == fSelectedPDG){
		//std::cout<<"found proton"<<std::endl;
		std::cout<<"PDG  = "<<pdg<<std::endl;
	}
	
	//check parent of particle
	//std::cout<<"The parent particle index is : "<<pfp.Parent()<<std::endl;

	//check number of daughters
	auto const& daughters = pfp.Daughters(); //get collection of daughter particles
	if (daughters.size() != 2){continue;} //want precisely 1 track 1 shower topology

	bool has_shower = false; //is true if at least one daughter is a shower
	bool has_track = false;	//is true if at least one daughter is a track
	//size_t pfp_daughter_trk_ind = -1; //index of PFP in PFPHandle corresponding to daughter track
	//size_t pfp_daughter_shr_ind = -1; //index of PFP in PFPHandle corresponding to daughter shower
	art::Ptr<recob::Track> this_track; //pointer to track
	art::Ptr<recob::Shower> this_shower; //pointer to shower

	std::cout<<"number of daughters = "<<daughters.size()<<std::endl;
	for (size_t d = 0; d != daughters.size(); ++d ){//for each daughter
		size_t daughter_index = daughters.at(d);
		//std::cout<<"the daughter particle(s) ID :"<<daughter_index<<std::endl;
		auto const& daughterPFP = PFPHandle->at(daughter_index);
		
		//check if the daughter is a track
		auto const& pfp_trk_v = pfp_trk_assn_v.at(daughter_index);
		std::cout<<"the number of tracks associated to this daughter = "<<pfp_trk_v.size()<<std::endl;
		if(pfp_trk_v.size() == 1){
			has_track = true;
			//get the track
			this_track = pfp_trk_v.at(0);
		}
		
		//check if the daughter is a shower
		auto const& pfp_shr_v = pfp_shr_assn_v.at(daughter_index);
		std::cout<<"the number of showers associated to this daughter = "<<pfp_shr_v.size()<<std::endl;
		if(pfp_shr_v.size() == 1){
			has_shower = true;
			this_shower = pfp_shr_v.at(0);
			//pfp_daughter_shr_ind = daughter_index;
	}
	
		std::cout<<"the daughter PDG code is :"<<daughterPFP.PdgCode()<<std::endl;
	}//for each daughter

	//skip if there isn't exactly one track and one shower
	if (has_shower == false || has_track == false){
		//std::cout<<"this event doesn't have one track and one shower"<<std::endl;
		continue;
	}else{
		my_trks.push_back(this_track);//add track pointer to vector
		my_shrs.push_back(this_shower);  //add shower pointer to vector
	}
	

	//if(pfp_daughter_shr_ind <= (size_t)-1 || pfp_daughter_trk_ind <= (size_t)-1){continue;}//check that track+shower vars were filled
 
	//get PFP corresponding to shower and track
	// auto const& shrPFP = PFPHandle->at(pfp_daughter_shr_ind);
	//auto const& trkPFP = PFPHandle->at(pfp_daughter_trk_ind);


	//for this event, get the vertex
	auto const& pfp_vtx_v = pfp_vtx_assn_v.at(pfp_index);

	//for each vertex
	for (size_t vtx_index = 0; vtx_index != pfp_vtx_v.size(); ++vtx_index ){
		auto vtx = *(pfp_vtx_v.at(vtx_index));
		std::cout<<"vertex ID  = "<<vtx.ID()<<std::endl;

		//get the vertex position
		double* xyz = new double[3];
        	vtx.XYZ(&xyz[0]);
        	_xpos = xyz[0];
       		_ypos = xyz[1];
        	_zpos = xyz[2];
        	
		//store the vertex postion
		my_vtxs.push_back(xyz);
		fmytree->Fill();
        	std::cout<<"vertex pos xyz = "<<_xpos<<", "<<_ypos<<", "<< _zpos<<std::endl;
	//	std::cout<<"the address to the vertex is "<<xyz<<std::endl;
		//delete[] xyz;
	} //for each vertex
}//for each PFP

std::cout<<"the number of stored vertices, tracks, and showers = "<<my_vtxs.size()<<", "<<my_trks.size()<<", "<<my_shrs.size()<<std::endl;

/* 
 * Remove shrhits matched with hits in the shower and look at shrhits matched to the track
 */


//std::cout<<"checking shower"<<std::endl;
//if(showerHandle->size() != 0) { //skip events with no shower
std::vector<TVector3> ShowerDir;
std::vector<double> ShowerLen;
std::vector<TVector3> ShowerStart;
 
  for ( size_t shr_index = 0; shr_index != showerHandle->size(); ++shr_index ){
	auto const& shr = showerHandle->at(shr_index);
	auto const& start  = shr.ShowerStart();
	for (size_t s = 0; s != my_shrs.size(); ++s ){
		auto const& current_shr = *(my_shrs.at(s));
		auto const& current_start = current_shr.ShowerStart();
		if (start == current_start){
			std::cout<<"matched showers"<<std::endl;
			std::cout<<"the shower length is  = "<<shr.Length()<<", and the opening angle is = "<<shr.OpenAngle()<<std::endl;
			
			 //get the associated hits for the shower
			 //const std::vector<art::Ptr<recob::Hit> > shr_hit_v = shr_hit_assn_v.at(shr_index);
			 auto shr_hit_v = shr_hit_assn_v.at(shr_index);
			fnum_hits_shower = shr_hit_v.size();
			std::cout<<"the number of hits in the shower = "<<fnum_hits_shower<<std::endl;		 
			int num_sshit_in_shower = 0;
			
			auto this_start = current_shr.ShowerStart();
			ShowerStart.push_back(this_start);

			auto this_len = current_shr.Length();
			ShowerLen.push_back(this_len);

			auto this_dir = current_shr.Direction();
			ShowerDir.push_back(this_dir);

			std::cout<<"the shower direction X component = "<<this_dir.X()<<std::endl;
		
			//for each hit
			 for (size_t h=0; h < shr_hit_v.size(); h++){
			 	auto const hit = *(shr_hit_v.at(h)); //get the hit from the pointer in the vector
				const recob::Hit* this_hit = &hit; //get the address to the hit that was stored in the vector
			
				//in the hit list, remove corresponding shr hits
				for(auto const& item : _hitlist){
					auto const stored_hit = (std::get<0>(item));
					if(matches(this_hit, stored_hit)==true){
						//std::cout<<"removing shower hit"<<std::endl;	
						_hitlist.remove(item);
						num_sshit_in_shower++;
						break;
					}	
				}

			}//for each hit
			fnum_sshits_shower = num_sshit_in_shower;

			if (num_sshit_in_shower > 300){
				std::cout<<"warning, very large number of sshits"<<std::endl;
			}
			std::cout<<"number of remaining matched shr hits = "<<_hitlist.size()<<std::endl;
			std::cout<<"the number of sshits in the shower = "<<fnum_sshits_shower<<std::endl;
			
			fratio_hits_shower = (double)fnum_sshits_shower/fnum_hits_shower;
			//std::cout<<"the ratio of sshits to hits is "<<fratio_hits_shower<<std::endl;
	//		fselectTree->Fill();
		}//if the showers match
	
	}//for each shower from a 1 shower 1 track topology 

	//auto const& length = shr.Length();
	//std::cout<<"shower length = "<<length<<std::endl;
	
//	auto const& start  = shr.ShowerStart();
//	std::cout<<"shower start xyz = "<<start.X()<<", "<<start.Y()<<", "<< start.Z()<<std::endl;
}//for each shower

//for each track
  for ( size_t trk_index = 0; trk_index != trackHandle->size(); ++trk_index ){
	auto const& trk = trackHandle->at(trk_index);
	auto const& start  = trk.Start();
	for (size_t tr = 0; tr != my_trks.size(); ++tr ){
		auto const& current_trk = *(my_trks.at(tr));
		auto const& current_start = current_trk.Start();
		int number_matched_trk_hits = 0;
		if (start == current_start){
			std::cout<<"matched tracks"<<std::endl;
			//std::cout<<"the shower length is  = "<<shr.Length()<<", and the opening angle is = "<<shr.OpenAngle()<<std::endl;
			
			 //get the associated hits for the track
			 const std::vector<art::Ptr<recob::Hit> > trk_hit_v = trk_hit_assn_v.at(trk_index);
			fnum_hits_track = trk_hit_v.size();
			std::cout<<"the number of hits in the track = "<<fnum_hits_track<<std::endl;		 
			
			//for each hit
			 for (size_t h=0; h < trk_hit_v.size(); h++){
			 	auto const hit = *(trk_hit_v.at(h));
				const recob::Hit* this_hit = &hit;
				//float this_channel = this_hit.Channel();
				
				//if the peak time is in the hit map, remove the shrhits from the map
				for (auto const& item : _hitlist){
					//_hitmap.erase(this_time);
					auto const stored_hit = (std::get<0>(item));
					if(matches(this_hit, stored_hit)== true){
						number_matched_trk_hits++;
					}
				}

			}//for each hit
		fnum_sshits_track = number_matched_trk_hits;
		fratio_hits_track = (double)fnum_sshits_track/fnum_hits_track;
		fselectTree->Fill();	
	
		std::cout<<"the number of matched shr hits in the track = "<<number_matched_trk_hits<<std::endl;
		}//if the tracks match
		
	}//for each track from a 1 shower 1 track topology 
}//for each track

/*
 * Look at remaining shrhits in ROI
 *
 */

//calc the radius for the ROI
//double radius = getRadius(fRadius);
	
//for each vertex position in the vector
for ( size_t vtx_index = 0; vtx_index != my_vtxs.size(); ++vtx_index ){
	//std::cout<<"looking at vertex # "<<vtx_index<<std::endl;
	double* xyz_pos = my_vtxs.at(vtx_index);
	//std::cout<<"the vertex pointer is "<<xyz_pos<<std::endl;
	
	//double X, Y, Z; 
	double X = *(xyz_pos);
	double Y = *(++xyz_pos);
	double Z = *(++xyz_pos);
	std::cout<<"the vertex XYZ = "<<X<<", "<<Y<<", "<<Z<<std::endl;
	
	// subset of hit <-> sshit pairs within ROI of the vertex
	std::list<std::pair<const recob::Hit*, const recob::Hit* >>  _ROIhitlist;
	
	//get the associated shower direction
	TVector3 shower_start = ShowerStart.at(vtx_index);
	double shower_length = ShowerLen.at(vtx_index);
	TVector3 shower_dir = ShowerDir.at(vtx_index);

	//create shower end point
	TVector3 shower_end = findEnd(&shower_start, &shower_length, &shower_dir);

	std::cout<<"The shower start is = "<<shower_start.X()<<", "<<shower_start.Y()<<", "<<shower_start.Z()<<std::endl;
	std::cout<<"The shower length is = "<<shower_length<<std::endl;
	std::cout<<"The calculated shower end is = "<<shower_end.X()<<", "<<shower_end.Y()<<", "<<shower_end.Z()<<std::endl;

	//for each plane
	for (int plane = 0; plane <fPlanes; ++plane){	
		//std::cout<<"starting plane "<<plane<<std::endl;
	
		//get the wire and time for each plane
		//grabbed this from lareventdisplay/EventDisplay/RecoBaseDrawer.cxx
		//int plane = 0;
		//double wire = geo->WireCoordinate(Y, Z, plane, fTPC, fCryostat);
		//double time = detprop->ConvertXToTicks(X, plane, fTPC,fCryostat);
		//std::cout<<"the time and wire on plane "<<plane<<" is ("<<time<<", "<<wire<<")"<<std::endl;

		double wire = calcWire(Y, Z, plane, fTPC, fCryostat, *fGeometry);
		double time = calcTime(X, plane, fTPC,fCryostat, *fDetprop);
		std::cout<<"the time and wire on plane "<<plane<<" is ("<<time<<", "<<wire<<")"<<std::endl;

		TVector2 shower_start_plane = TVector2(calcTime(shower_start.X(), plane, fTPC,fCryostat, *fDetprop), calcWire(shower_start.Y(), shower_start.Z(), plane, fTPC, fCryostat, *fGeometry));	
		TVector2 shower_end_plane =  TVector2(calcTime(shower_end.X(), plane, fTPC,fCryostat, *fDetprop), calcWire(shower_end.Y(), shower_end.Z(), plane, fTPC, fCryostat, *fGeometry));

		TVector2 shower_direction_plane = calcNormVec(shower_start_plane, shower_end_plane); 
		
		//std::cout<<"the direction of the shower on this plane is "<< shower_direction_plane.X() <<" starting and ending from "<< shower_start_plane.X()<<shower_end_plane.X() <<std::endl;
		//for each remaining sshit
		for(auto const& item : _hitlist){
			//if the sshit falls inside the ROI
			auto const this_sshit = std::get<1>(item);
			double dist = distToVtx(fTimeToCMConstant, fWireToCMConstant, fRadius, time, wire, this_sshit);
			if (inROI(fRadius, dist) == true){
			//	if(_ROIhitlist.size() < 3){
			//		std::cout<<"corresponding shower dir x = "<<shower_dir.X()<<std::endl;
			//	}
				_ROIhitlist.push_back(item);
                                fradial_dist_sshit_vtx = dist;
				fopening_angle_shower_sshit = angleFromShower(shower_direction_plane, this_sshit, time,wire);
				//std::cout<<"the angle from the shower is "<<angle<<std::endl;
				
				fROITree->Fill();	
			} //if in ROI
		}//each sshit

		//std::cout<<"ending plane "<<plane<<std::endl;

	}//for each plane
fnum_sshits_ROI_no_shower = _ROIhitlist.size();
std::cout<<"the number of shower hits within the ROI is after removing matched =  "<<fnum_sshits_ROI_no_shower<<std::endl;
fROITree->Fill();

}//loop over vertices

total_num_vertices+= my_vtxs.size();
std::cout<<"the number processed ROI's in file = "<<total_num_vertices<<std::endl;

 auto particleHandle
      = event.getValidHandle<std::vector<simb::MCParticle>>
      (fSimulationProducerLabel);
  
 if (!particleHandle){std::cout<<"missing particle handle"<<std::endl;}
    
 
std::map< int, const simb::MCParticle* > particleMap;
//std::cout<<"flag 3.1"<<std::endl;
    for ( auto const& particle : (*particleHandle) )
      {
	// For the methods you can call to get particle information,
	// see ${NUTOOLS_INC}/SimulationBase/MCParticle.h.
	fSimTrackID = particle.TrackId();
	//std::cout<<"TrackID: "<<fSimTrackID<<std::endl;
//std::cout<<"flag 3.1.1"<<std::endl;
	// Add the address of the MCParticle to the map, with the track ID as the key.
	particleMap[fSimTrackID] = &particle; 
//std::cout<<"flag 3.1.2"<<std::endl;
	// Histogram the PDG code of every particle in the event.
	fSimPDG = particle.PdgCode();
	//if(fSimPDG == fSelectedPDG){std::cout<<"Found particle with PDG: "<<fSimPDG<<std::endl;}
	fPDGCodeHist->Fill( fSimPDG );
//std::cout<<"flag 3.2"<<std::endl;


	//if ( particle.Process() == "primary"  &&  fSimPDG == fSelectedPDG )
	if (fSimPDG == fSelectedPDG ) { //for all events containing a proton
	  //std::cout<<"flag 3.2.1"<<std::endl;
	    // A particle has a trajectory, consisting of a set of
	    // 4-positions and 4-mommenta.
	    const size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

	    // For trajectories, as for vectors and arrays, the
	    // first point is #0, not #1.
	    const int last = numberTrajectoryPoints - 1;
	    const TLorentzVector& positionStart = particle.Position(0);
	    const TLorentzVector& positionEnd   = particle.Position(last);
	    const TLorentzVector& momentumStart = particle.Momentum(0);
	    const TLorentzVector& momentumEnd   = particle.Momentum(last);
//std::cout<<"flag 3.2.2"<<std::endl;

	    // Make a histogram of the starting momentum.
	   // std::cout<<"P: "<<momentumStart.P()<<std::endl;
	  //  fMomentumHist->Fill( momentumStart.P() );
//std::cout<<"flag 3.2.2.1"<<std::endl;

	    // Fill arrays with the 4-values. (Don't be fooled by
	    // the name of the method; it just puts the numbers from
	    // the 4-vector into the array.)
	    positionStart.GetXYZT( fStartXYZT );
//std::cout<<"flag 3.2.2.2"<<std::endl;

	    positionEnd.GetXYZT( fEndXYZT );
	//std::cout<<"flag 3.2.2.3"<<std::endl;

	    momentumStart.GetXYZT( fStartPE );
	//std::cout<<"flag 3.2.2.4"<<std::endl;

	    momentumEnd.GetXYZT( fEndPE );
//std::cout<<"flag 3.2.3"<<std::endl;

	    // Use a polar-coordinate view of the 4-vectors to
	    // get the track length.
	    const double trackLength = ( positionEnd - positionStart ).Rho();
	    LOG_DEBUG("SSNetTest")
	      << "Track length: " << trackLength << " cm";
	    
	    // Fill a histogram of the track length.
	    fTrackLengthHist->Fill( trackLength ); 
//std::cout<<"flag 3.2.4"<<std::endl

            LOG_DEBUG("SSNetTest")
	      << "track ID=" << fSimTrackID << " (PDG ID: " << fSimPDG << ") "
	      << trackLength << " cm long, momentum " << momentumStart.P();
	   //std::cout<<"filling tree"<<std::endl;
	    fmytree->Fill();

	  } // if primary and PDG selected by user
      } // loop over all particles in the event. 
 //   std::cout<<"flag 4"<<std::endl;
    // std::cout<<"flag 6"<<std::endl;
    // We have a map of dE/dx vectors. Write each one of them to the
    // reconstruction n-tuple. 
    //for ( const auto& dEdxEntry : dEdxMap )
      //{
	// At this point, we've filled in all the reconstruction
	// n-tuple's variables. Write it.
	//fReconstructionNtuple->Fill();
      //}
 
    // First, read in the clusters.
    // Again, here we require clusters to be available by getting a ValidHandle.
    // Otherwise, the following code would not work.
      // const art::FindManyP<recob::PFParticle> findManyPFP(clusterHandle, event, fClusterProducerLabel);

    const art::FindManyP<recob::Hit> findManyHits(clusterHandle, event, fClusterProducerLabel);

    if ( findManyHits.isValid() )
      {
        for ( size_t cluster_index = 0; cluster_index != clusterHandle->size(); ++cluster_index )
	  {
            auto const& hits = findManyHits.at( cluster_index );
	    mf::LogInfo("SSNetTest")  
	      << "Cluster ID=" << clusterHandle->at( cluster_index ).ID()
	      << " has " << hits.size() << " hits";
	  }
      } // findManyHits valid
    else
      {
	mf::LogError("SSNetTest")  
	  << "findManyHits recob::Hit for recob::Cluster failed;"
	  << " cluster label='" << fClusterProducerLabel << "'";
      }
 } // SSNetTest::analyze()
 
  
  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see SSNetTest.fcl for more information.
  DEFINE_ART_MODULE(SSNetTest)

} // namespace example
} // namespace lar 

namespace {
 // time to define that function...
 double DetectorDiagonal(geo::GeometryCore const& geom) {
	const double length = geom.DetLength();
	const double width = 2. * geom.DetHalfWidth();
	const double height = 2. * geom.DetHalfHeight();
	return std::sqrt(cet::sum_of_squares(length, width, height));
	} // DetectorDiagonal()

//takes pointers to any two hits (hit/sshit or hit/hit etc.) and returns true if they match
 bool matches(const recob::Hit* Hit, const recob::Hit* ssHit){
	auto hit = *Hit;
	auto sshit = *ssHit;
	
	auto t = hit.PeakTime();
        auto w = hit.WireID().Wire;
        auto plane  = hit.View();

	auto sst = sshit.PeakTime();
        auto ssw = sshit.WireID().Wire;
        auto ssplane  = sshit.View();

	if(t == sst && w == ssw && plane == ssplane){
		return true;
	} else{
		return false;
	}
}

//takes the plane and 2D vertex position and pointer to a 3D hit
//returns the radial distance in cm of the hit from the vertex in cm
double distToVtx(double fTimeToCMConstant, double fWireToCMConstant, double radius, double vertex_time, double vertex_wire, const recob::Hit*  hit){
	double hit_time = (*hit).PeakTime();
	double hit_wire = (*hit).WireID().Wire;
	double diff_t = hit_time - vertex_time;
	double diff_w = hit_wire - vertex_wire;

	//convert to cm
	diff_t = diff_t * fTimeToCMConstant;
	diff_w = diff_w * fWireToCMConstant;
	
	//calc dist from vertex
	double dist_from_vtx = (diff_t*diff_t) + (diff_w*diff_w);
	return dist_from_vtx;
}//distToVtx

//takes the distance to the vertex of a hit on a given plane and the radius
//returns true if hit is in ROI, false otherwise
bool inROI(double radius, double dist ){
	if(dist <= radius){
	//	std::cout<<"matched sshit at ("<<hit_time<<", "<<hit_wire<<")"<<std::endl;
//		std::cout<<"the distance from the vertex is "<<dist_from_vtx<<std::endl;
		return true;
	}else{
		return false;
	}
}//inROI

//given the shower start, direction, and length, calculates an approximate end point in 3D
TVector3 findEnd(TVector3* shower_start, double* shower_length, TVector3* shower_dir){
	//scale the direction by the length
	TVector3 scaled_dir = (*shower_length) * (*shower_dir);

	//add it to the start
	TVector3 end = (*shower_start)+scaled_dir;
	return end;
}


//takes Y and Z  coord in cm for a 3D spacepoint and returns the wire coord for a given plane
double calcWire(double Y, double Z, int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo ){
	double wire = geo.WireCoordinate(Y, Z, plane, fTPC, fCryostat);
	return wire;
}

//takes X coord for a 3D spacepoint and returns the time coord for a given plane
double calcTime(double X,int plane,int fTPC,int fCryostat, detinfo::DetectorProperties const& detprop){
	double time = detprop.ConvertXToTicks(X, plane, fTPC,fCryostat);
	return time;
}

//takes 2 vectors and returns normalized vector between them pointing from start to end
TVector2 calcNormVec(TVector2 shower_start_plane,TVector2 shower_end_plane){
	TVector2 diff = shower_end_plane - shower_start_plane;
	return diff.Unit();
}

//for a given plane, takes shower end points (2D projected) and then calculates angle
//between vector from vertex to hit and the shower
double angleFromShower(TVector2 shower_direction_plane, const recob::Hit*  hit, double vertex_time, double vertex_wire){
 	double hit_time = (*hit).PeakTime();
        double hit_wire = (*hit).WireID().Wire;	

	TVector2 this_hit = TVector2(hit_time, hit_wire);
	TVector2 this_vertex = TVector2(vertex_time, vertex_wire);

	TVector2 hit_vtx_direction_plane = calcNormVec(this_vertex, this_hit);

	double opening_angle = hit_vtx_direction_plane.Phi() - shower_direction_plane.Phi();
	return opening_angle;
	
}


//converts radius for ROI in units of time and wire from sm
//double getRadius(double rad_in_cm){
	//convert cm to time
	//convert cm to wire
	//r = sqrt(t^2 + w^2)
//	return 500;
//}//getRadius


} // local namespace
