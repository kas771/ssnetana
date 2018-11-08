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
 #include "lardata/Utilities/GeometryUtilities.h"
 
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
 #include "TF1.h"
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
 #include <iostream>
 #include <fstream> 

namespace{
double DetectorDiagonal(geo::GeometryCore const& geom);
bool matches(const recob::Hit* Hit, const recob::Hit* ssHit);
double distToVtx(double fTimeToCMConstant, double fWireToCMConstant, double radius, double vertex_time, double vertex_wire, const recob::Hit*  hit);
//double distToVtx2d(TVector2 my_point, TVector2 vertex);
bool inROI(double radius, double dist);
std::vector<TVector2> drawROI(double fTimeToCMConstant, double fWireToCMConstant, double radius, double vertex_time, double vertex_wire);
std::vector<TVector2> drawEllipse(double start_x, double start_y, double ax_x, double ax_y);
double calcWire(double Y, double Z, int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo );
TVector3 findEnd(TVector3* shower_start, double* shower_length, TVector3* shower_dir);
double calcTime(double X,int plane,int fTPC,int fCryostat, detinfo::DetectorProperties const& detprop);
TVector2 calcNormVec(TVector2 shower_start_plane,TVector2 shower_end_plane);
double angleFromShower(TVector2 shower_direction_plane, const recob::Hit*  hit, double vertex_time, double vertex_wire);
//std::list<std::pair<const recob::Hit*, const recob::Hit* >> removeHitsROIList(std::list<std::pair<const recob::Hit*, const recob::Hit* >> _ROIlist, std::list<std::pair<const recob::Hit*, const recob::Hit* >> _showerlist);

std::vector<TVector2> getNumHitsSectors(int n_inc,TVector2 shower_end, double radius, double vertex_time, double vertex_wire, std::list<std::pair<const recob::Hit*, const recob::Hit* >> _ROIlist, double fTimetoCMConstant, double fWiretoCMConstant);

bool inSector(TVector2 sector_start_cm, TVector2 sector_end_cm,  const recob::Hit*  hit, double vertex_time, double vertex_wire, double fTimetoCMConstant, double fWiretoCMConstant);

TVector2 calcCM(TVector2 vec_time_wire, double fTimeToCMConstant, double fWireToCMConstant);

bool areClockwise(TVector2 v1, TVector2 v2);

TVector2 getSectorEnd(TVector2 sector_start_cm, double theta);

std::vector<int> getIndHitsROIList(std::list<std::pair<const recob::Hit*, const recob::Hit* >> _ROIlist, std::list<std::pair<const recob::Hit*, const recob::Hit* >> _showerlist);

std::vector<std::array<double, 3>> Circle3D(const TVector3& centerPos, const TVector3& axisDir, const double& radius);

//std::vector<TVector2> fillROITree(std::list<std::pair<const recob::Hit*, const recob::Hit* >> _listnotrack, std::vector<TVector2>vertex_wire_time_plane, std::vector<TVector2> shower_dir_plane, int fPlanes, double fTimeToCMConstant, double fWireToCMConstant, double fRadius);

int contains_at_ind(int n, std::vector<int> list);

bool contains(int n, std::vector<int> _listnotrack);

//double fitFunction(double *x, double *par);
double mygaus(double x, double A, double mean, double sigma);
double fitf(Double_t *x,Double_t *par);

std::vector<TVector2> getMinMaxShowerPlane(std::vector<std::array<double, 3>> coneRim,int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo, detinfo::DetectorProperties const& detprop, TVector2 shower_start_plane, std::vector<TVector2> min_max);
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
  
   fhicl::Atom<std::string> InputVertexFile {
        Name("InputVertexFile"),
        Comment(".txt file containing information about the single photon vertices")
	};
   
    fhicl::Atom<int> Run {
        Name("Run"),
        Comment("The run for the event to for the output histograms")
	};

    fhicl::Atom<int> Subrun {
        Name("Subrun"),
        Comment("The subrun for the event to for the output histograms")
	};

    fhicl::Atom<int> Event {
        Name("Event"),
        Comment("The event number for the event to for the output histograms")
	};

};//struct
    
     using Parameters = art::EDAnalyzer::Table<Config>;

    explicit SSNetTest(Parameters const& config); 
  
    virtual void beginJob() override;
    virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;
    virtual void endJob() override;

  private:
  	
   //io stream to write to .txt file for EVD      
   std::ofstream out_stream;

   //io stream to read from .txt file containing single photon vertex, shower, and track params
   std::ifstream in_stream;
    
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
    std::string fInputVertexFile;//contains list of vertex for events
    int fRun_hist; //the run number to make histograms for
    int fSubrun_hist; //the subrun number to make histograms for
    int fEvent_hist; //the event number to make histograms for

// vector of shower-like hit indices
//   std::vector<size_t> _shrhits;

 // Pointers to the histograms we'll create. 
    TH1D* fPDGCodeHist;     ///< PDG code of all particles
    TH1D* fMomentumHist;    ///< momentum [GeV] of all selected particles
    TH1D* fTrackLengthHist; ///< true length [cm] of all selected particlesi
    TH1D* my_hist; ///sshits in ROI

    //Histograms for the event indicated in the fcl
    TH1D* fAngle_sshits; //the angular dist for all sshits in the ROI
    TH1D* fDist_sshits; //the radial dist for all sshits in the ROI
    TH1D* fAngle_sshits_noshower; //the angular dist for all sshits in the ROI minus the shower
    TH1D* fDist_sshits_noshower; //the radial dist for all sshits in the ROI minus the shower
   
    TH1D* fAngle_sshits_plane0;//same as above but only for hits on plane 0
    TH1D* fDist_sshits_plane0;
    TH1D* fAngle_sshits_noshower_plane0;
    TH1D* fDist_sshits_noshower_plane0;
    TH1D* fAngle_sshits_plane1;//same as above but only for hits on plane 1
    TH1D* fDist_sshits_plane1;
    TH1D* fAngle_sshits_noshower_plane1;
    TH1D* fDist_sshits_noshower_plane1;
    TH1D* fAngle_sshits_plane2;//same as above but only for hits on plane 2
    TH1D* fDist_sshits_plane2;
    TH1D* fAngle_sshits_noshower_plane2;
    TH1D* fDist_sshits_noshower_plane2;



    // The n-tuples we'll create.
    TTree* fmytree;     ///< tuple with information about all events
    TTree* fselectTree; ///< contains information about 1 track 1 shower events (each ROI)
    TTree* fROITree; ///<contains information each hit inside the ROI

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
   
    int fnum_sshits_ROI; //total number of sshits within each ROI, including both track and shower
    int fnum_sshits_ROI_no_shower; //total number of sshits within each ROI once only the shower is removed
    int fnum_sshits_ROI_no_track; //total number of sshits within each ROI once only the track is removed
    int fnum_sshits_ROI_no_track_no_shower; //total number of sshits within each ROI once both the track and shower are removed
       
    int fnum_sshits_shower; //number of sshits in a shower
    int fnum_hits_shower; //number of hits in a shower
    double fratio_hits_shower; //ratio of sshits to hits in a shower

    int fnum_sshits_track; //number of sshits in a track
    int fnum_hits_track; //number of hits in a track
    double fratio_hits_track; //ratio of sshits to hits in a track

    double fradial_dist_sshit_vtx; //the distance of each sshit in ROI to the vertex (all sshits)
    double fopening_angle_shower_sshit; //the angle on each plane between the shower direction and the vtx-sshit
    double fradial_dist_sshit_vtx_noshower; //the distance of each sshit in ROI to the vertex (shower removed)
    double fopening_angle_shower_sshit_noshower; //the angle on each plane between the shower direction and the vtx-sshit
    double fradial_dist_sshit_vtx_notrack; //the distance of each sshit in ROI to the vertex (track removed)
    double fopening_angle_shower_sshit_notrack; //the angle on each plane between the shower direction and the vtx-sshit
    double fradial_dist_sshit_vtx_noshowernotrack; //the distance of each sshit in ROI to the vertex (shower and track removed)
    double fopening_angle_shower_sshit_noshowernotrack; //the angle on each plane between the shower direction and the vtx-sshit
	 
    //int fnum_sshits_ROI_with_shower; 
    int n_inc;
    double angle_min;
    double angle_max;
    geo::GeometryCore const* fGeometry;
      //stuff to do the vertex to plane mapping
   // art::ServiceHandle<geo::Geometry>  geo;
  //  geo::Geometry const* fGeo;
    detinfo::DetectorProperties const* fDetprop;    
    //uiil::GeometryUtilities const* gser;


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
	, fInputVertexFile	       	  (config().InputVertexFile())
	, fRun_hist		       	  (config().Run())
	, fSubrun_hist		       	  (config().Subrun())
	, fEvent_hist		       	  (config().Event())
	
	// fInHitProducer   = p.get<std::string>("InHitProducer","gaushit");
        //fPxThresholdHigh = p.get<double>     ("PxThresholdHigh"        );
        //fPxThresholdLow  = p.get<double>     ("PxThresholdLow"         );
       // produces< std::vector< recob::Hit > >(Form("shrhit%zu",(size_t)(fPxThresholdHigh*100.)));
        //produces< std::vector< recob::Hit > >(Form("shrhit%zu",(size_t)(fPxThresholdLow*100.) ));


  {
    // get a pointer to the geometry service provider
     fGeometry = lar::providerFrom<geo::Geometry>();
     //fGeo =  lar::providerFrom<geo::Geometry>();
   // util::GeometryUtilities gser;

    //get a pointer to the detector properties service provider
     fDetprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
   //  gser = lar::providerFrom<util::GeometryUtilities>();

    //initialize io stream for printing to output
    out_stream.open("test.txt");
    if (!out_stream.is_open()){
	std::cout<<"ERROR output file not open"<<std::endl;
	exit(0);
	}

    //initialize io stream for reading from input, quit if file not opened properly
    in_stream.open(fInputVertexFile);
    if (!in_stream.is_open()){
      std::cout<<"ERROR input file not open: "<<fInputVertexFile<<std::endl;
      exit(0);
    }

 }


  
  //-----------------------------------------------------------------------
  void SSNetTest::beginJob()
  {
   // fb.open ("test.txt",std::ios::out);
  //  os(&fb);
    out_stream << "Test sentence\n";
    out_stream<<"radius "<<fRadius<<" wiretocm "<<fWireToCMConstant<<" timetocm "<<fTimeToCMConstant<<"\n";
    //fb.close();
   
    //util::GeometryUtilities gser;

    //double timetocm = (*gser).TimeToCm();
    //double wiretocm = (*gser).WireToCm();
   // double timetocm = (*gser).TimeToCm();
   // double wiretocm = (*gser).WireToCm();


  //  std::cout<<"time to cm = "<<timetocm<<", wire to cm = "<<wiretocm<<std::endl;

    // Get the detector length, to determine the maximum bin edge of one
    // of the histograms.
    const double detectorLength = DetectorDiagonal(*fGeometry);
  // auto const GeotimetoCM = gser.TimetoCM();
  // std::cout<<"the geometry time to cm is "<<GeotimetoCM<<std::endl
 
   //get the tpc and cryostat ID's  
    auto const TPC = (*fGeometry).begin_TPC();
    auto ID = TPC.ID();
    fCryostat = ID.Cryostat;
    fTPC = ID.TPC;
    std::cout<<TPC.ID()<<"= the beginning TPC ID" <<std::endl;
    std::cout<<"the cryostat id = "<<fCryostat<<std::endl;  
    std::cout<<"the tpc id = "<<fTPC<<std::endl;  
    //std::cout<<"this version has compiled"<<std::endl;
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

    n_inc = 100;
    angle_min = -180;
    angle_max = 180;

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;
  
    fPDGCodeHist     = 		tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
    fMomentumHist    = 		tfs->make<TH1D>("mom",     ";particle Momentum (GeV);",    100, 0.,    10.);
    fTrackLengthHist = 		tfs->make<TH1D>("length",  ";particle track length (cm);", 200, 0, detectorLength);
    my_hist = 			tfs->make<TH1D>("num sshits in roi",  ";;", 200, 0, 1000 );
    fAngle_sshits =  		tfs->make<TH1D>("angle_all",";Opening Angle all SSNet Hits Roi;",100, -3.3, 3.3);
    fDist_sshits =		tfs->make<TH1D>("dist_all",";Distance from Vertex all SSNet Hits Roi;", 100, 0, fRadius);
    fAngle_sshits_noshower = 	tfs->make<TH1D>("angle_noshower",";Opening Angle SSNet Hits Roi minus Shower;",100, -3.3, 3.3);
    fDist_sshits_noshower = 	tfs->make<TH1D>("dist_noshower", ";Opening Angle all SSNet Hits Roi;",100, 0, fRadius);
    fAngle_sshits_plane0 =  		tfs->make<TH1D>("angle_all_plane0",";Opening Angle all SSNet Hits Roi for Plane0;",n_inc, angle_min, angle_max);
    fDist_sshits_plane0 =		tfs->make<TH1D>("dist_all_plane0",";Distance from Vertex all SSNet Hits Roi for Plane0;",100, 0, fRadius);
    fAngle_sshits_noshower_plane0 = 	tfs->make<TH1D>("angle_noshower_plane0",";Opening Angle SSNet Hits Roi minus Shower for Plane0;",100, -3.3, 3.3);
    fDist_sshits_noshower_plane0  = 	tfs->make<TH1D>("dist_noshower_plane0", ";Opening Angle all SSNet Hits Roi for Plane0;",100, 0, fRadius);
    fAngle_sshits_plane1 =  		tfs->make<TH1D>("angle_all_plane1",";Opening Angle all SSNet Hits Roi for Plane1;",n_inc, angle_min, angle_max);
    fDist_sshits_plane1 =		tfs->make<TH1D>("dist_all_plane1",";Distance from Vertex all SSNet Hits Roi for Plane1;",100, 0, fRadius);
    fAngle_sshits_noshower_plane1 = 	tfs->make<TH1D>("angle_noshower_plane1",";Opening Angle SSNet Hits Roi minus Shower for Plane1;",100, -3.3, 3.3);
    fDist_sshits_noshower_plane1  = 	tfs->make<TH1D>("dist_noshower_plane1", ";Opening Angle all SSNet Hits Roi for Plane1;",100, 0, fRadius);
    fAngle_sshits_plane2 =  		tfs->make<TH1D>("angle_all_plane2",";Opening Angle all SSNet Hits Roi for Plane2;",n_inc, angle_min,angle_max);
    fDist_sshits_plane2 =		tfs->make<TH1D>("dist_all_plane2",";Distance from Vertex all SSNet Hits Roi for Plane2;",100, 0, fRadius);
    fAngle_sshits_noshower_plane2 = 	tfs->make<TH1D>("angle_noshower_plane2",";Opening Angle SSNet Hits Roi minus Shower for Plane2;",100, -3.3, 3.3);
    fDist_sshits_noshower_plane2  = 	tfs->make<TH1D>("dist_noshower_plane2", ";Opening Angle all SSNet Hits Roi for Plane2;",100, 0, fRadius);


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
    fselectTree->Branch("num_sshits_ROI_no_track_no_shower", &fnum_sshits_ROI_no_track_no_shower,   "num_sshits_ROI_no_track_no_shower/I");
    fselectTree->Branch("num_sshits_ROI_no_shower", &fnum_sshits_ROI_no_shower,   "num_sshits_ROI_no_shower/I");
    fselectTree->Branch("num_sshits_ROI_no_track", &fnum_sshits_ROI_no_track,   "num_sshits_ROI_no_track/I");
    fselectTree->Branch("num_sshits_ROI", &fnum_sshits_ROI,   "num_sshits_ROI/I");
    
    fselectTree->Branch("num_sshits_shower", &fnum_sshits_shower,   "num_sshits_shower/I");
    fselectTree->Branch("num_hits_shower", &fnum_hits_shower,   "num_hits_shower/I");
    fselectTree->Branch("ratio_hits_shower", &fratio_hits_shower,   "ratio_hits_shower/D");
    
    fselectTree->Branch("num_sshits_track", &fnum_sshits_track,   "num_sshits_track/I");
    fselectTree->Branch("num_hits_track", &fnum_hits_track,   "num_hits_track/I");
    fselectTree->Branch("ratio_hits_track", &fratio_hits_track,   "ratio_hits_track/D");


    fROITree->Branch("radial_dist_sshit_vtx", &fradial_dist_sshit_vtx, "radial_dist_sshit_vtx/D");
    fROITree->Branch("fopening_angle_shower_sshit", &fopening_angle_shower_sshit, "opening_angle_shower_sshit/D");
    fROITree->Branch("radial_dist_sshit_vtx_noshower", &fradial_dist_sshit_vtx_noshower, "radial_dist_sshit_vtx_noshower/D");
    fROITree->Branch("fopening_angle_shower_sshit_noshower", &fopening_angle_shower_sshit_noshower, "opening_angle_shower_sshit_noshower/D");
    fROITree->Branch("radial_dist_sshit_vtx_notrack", &fradial_dist_sshit_vtx_notrack, "radial_dist_sshit_vtx_notrack/D");
    fROITree->Branch("fopening_angle_shower_sshit_notrack", &fopening_angle_shower_sshit_notrack, "opening_angle_shower_sshit_notrack/D");
    fROITree->Branch("radial_dist_sshit_vtx_noshowernotrack", &fradial_dist_sshit_vtx_noshowernotrack, "radial_dist_sshit_vtx_noshowernotrack/D");
    fROITree->Branch("fopening_angle_shower_sshit_noshowernotrack", &fopening_angle_shower_sshit_noshowernotrack, "opening_angle_shower_sshit_noshowernotrack/D");

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
    //std::cout<<"------- starting event --------"<<std::endl;
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    
    std::cout<<"------- starting event: "<<fRun<<", "<<fSubRun<<", "<<fEvent<<" --------"<<std::endl;
     
 /*
 *
 * Find matching event in the input file and read in vertex, shower, and track
 *
 * */
   
    //initialize line
    int nlines = 0;
    std::string line;
    line.clear();
 

    //set to start of the input file
    in_stream.clear();
    in_stream.seekg(0, in_stream.beg);
    //std::cout<<"the position in the sequnce is "<<in_stream.tellg()<<std::endl; 
   
    //find corresponding line in .txt input file for this event
    bool matched= false;	
    std::vector<std::string> words_in_line;
    while ( getline (in_stream,line) )
    {
	words_in_line.clear();	
	int this_word= 0;
	std::istringstream iss(line);
//	std::vector<std::string> words_in_line;
	do {
        	std::string subs;
        	iss >> subs;
        	//if (word == 0 && subs!= std::to_string(fRun)){continue;}
		//std::cout << "Substring at word "<<this_word<<" is: " << subs <<std::endl;
		words_in_line.push_back(subs);	
		this_word++;
    	} while (iss);  

	if (words_in_line[0] == std::to_string(fRun) && words_in_line[1] == std::to_string(fSubRun) && words_in_line[2]==std::to_string(fEvent)){
		std::cout<<"Matched: "<<words_in_line[0]<<", "<<words_in_line[1]<<", "<<words_in_line[2]<<std::endl;
		matched = true;
		break;
		}
   	nlines++;
//	std::cout << line <<std::endl;
    }

    //if there isn't a match in the vertex file, skip this event
    if (matched == false){
	std::cout<<"ERROR: no match for this event was found: "<<std::endl;
	std::cout<<fRun<<", "<<fSubRun<<", "<<fEvent<<std::endl;
	std::cout<<"Skipping this event"<<std::endl;
//	break;
	//exit(0);
    }else{
	
    
    TVector3 single_vertex = TVector3(std::stof (words_in_line[3], 0), std::stof (words_in_line[4], 0), std::stof (words_in_line[5], 0));
    //double vertex_y = std::stof (words_in_line[4], 0);
    //double vertex_z = std::stof (words_in_line[5], 0);

    _xpos = single_vertex.X();
    _ypos = single_vertex.Y();
    _zpos = single_vertex.Z();

    TVector3 single_track_vertex_dir = TVector3(std::stof (words_in_line[6], 0), std::stof (words_in_line[7], 0), std::stof (words_in_line[8], 0));
    TVector3 single_shower_start = TVector3(std::stof (words_in_line[9], 0), std::stof (words_in_line[10], 0), std::stof (words_in_line[11], 0));
    TVector3 single_shower_dir = TVector3(std::stof (words_in_line[12], 0), std::stof (words_in_line[13], 0), std::stof (words_in_line[14], 0));

    std::cout<<"single photon vertex = "<<single_vertex.X()<<", "<<single_vertex.Y()<<", "<<single_vertex.Z()<<std::endl;
    std::cout<<"single photon track vertex = "<<single_track_vertex_dir.X()<<", "<<single_track_vertex_dir.Y()<<", "<<single_track_vertex_dir.Z()<<std::endl;
    std::cout<<"single photon shower start = "<<single_shower_start.X()<<", "<<single_shower_start.Y()<<", "<<single_shower_start.Z()<<std::endl;

    //close input .txt file
    //in_stream.close();
    std::cout<<"the number of events in the input file is: "<<nlines<<std::endl; 
//    std::cout<<"the number of characters in the input file is: "<<in_stream.gcount()<<std::endl;
  
    //write run subrun event to .txt
    //std::ostream os(&fb);
    out_stream <<"run subrun event "<<fRun<<" "<<fSubRun<<" "<<fEvent<<"\n";

/*
 *
 * Read in associations for Pandora objects
 *
 * */

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
std::list<std::pair<const recob::Hit*, const recob::Hit* >> _hitlist; //a list of all sshits<->hits in the event
std::list<std::pair<const recob::Hit*, const recob::Hit* >> _showerhitlist; //a list of sshits<->hits in the shower associated with vertex (not exclusively in ROI)
std::list<std::pair<const recob::Hit*, const recob::Hit* >> _trackhitlist; //a list of sshits<->hits in the track associated with vertex (not exclusively in ROI)
	
int n = 0; //the total number of sshits
fnum_total_sshits = sshitHandle->size();

// For every ssHit:
    for ( auto const& sshit : (*sshitHandle) )
      {
	int plane = sshit.View();
	double wire = sshit.WireID().Wire;
	double time = sshit.PeakTime();
	out_stream<<"sshit "<<n<<" plane "<<plane<<" wire time "<<wire <<" "<<time<<"\n";
	n++;
	//for every Hit
	for ( auto const& hit : (*hitHandle) ){
		//if (n==1){
		//	int plane = hit.View();
        	//	double wire = hit.WireID().Wire*fWireToCMConstant;
        	//	double time = hit.PeakTime()* fTimeToCMConstant;
		//}
		//if the wire, plane, and time match
		if (matches(&hit, &sshit)== true){	
			//make a pair of pointers to hits and add them to the list
			_hitlist.emplace(_hitlist.begin(), &hit, &sshit);
		}//if they match	
	}//for each Hit
 } // for each SSHit

//add hits to .txt file
int n2 = 0;
for ( auto const& hit : (*hitHandle) ){
	int plane = hit.View();
	double wire = hit.WireID().Wire;
        double time = hit.PeakTime();
	out_stream<<"hit "<<n2<<" plane "<<plane<<" wire time "<<wire <<" "<<time<<"\n";
	n2++;
}

fnum_total_matched_hits = _hitlist.size();
fratio_total_sshits_hits = (double)fnum_total_sshits/fnum_total_hits;
std::cout<<"number of hits = "<<fnum_total_hits<<", number of sshits = "<<n<<"/"<<fnum_total_sshits<<", number of matches = "<<fnum_total_matched_hits<<std::endl;
fmytree->Fill();
std::cout<<"The number of entries = "<<fmytree->GetEntries()<<std::endl;
 
/*
 * loop over tracks and showers in and find the 1 shower 1 track for each single photon vertex
 *
 */

  //std::vector<art::Ptr<recob::Track>> my_trks; //vector of pointers to tracks in case of 1 shower 1 track events
  //std::vector<art::Ptr<recob::Shower>> my_shrs; //vector of pointers to tracks in case of 1 shower 1 track events
  //std::vector<double*> my_vtxs; //vector of pointers to the first element in a double for the vertex XYZ position (3D)

 
 // std::cout<<"number of PFParticles: "<<PFPHandle->size()<<std::endl;
  //for each PFP
/*
 std::cout<<PFPHandle->size()<<std::endl;
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
		//fmytree->Fill();
        	std::cout<<"vertex pos xyz = "<<_xpos<<", "<<_ypos<<", "<< _zpos<<std::endl;
	//	std::cout<<"the address to the vertex is "<<xyz<<std::endl;
		//delete[] xyz;
	} //for each vertex
}//for each PFP
*/


//std::cout<<"the number of stored vertices, tracks, and showers = "<<my_vtxs.size()<<", "<<my_trks.size()<<", "<<my_shrs.size()<<std::endl;

/* 
 * Find the Pandora track and shower corresponding to the single photon vertex
 * and associate shrhits matched with hits in the shower and track
 *
 */

//std::cout<<"checking shower"<<std::endl;
//if(showerHandle->size() != 0) { //skip events with no shower
std::vector<TVector3> ShowerDir;
std::vector<double> ShowerLen;
std::vector<TVector3> ShowerStart;
std::vector<double> ShowerAng;

bool matched_showers =false;
int number_potential_shower_matches = 0;
int matched_shower_ind = -1; //the index of the track that best matches the single photon track
	
//for each shower in the event 
  for ( size_t shr_index = 0; shr_index != showerHandle->size(); ++shr_index ){
	auto const& shr = showerHandle->at(shr_index);
	//auto const& start  = shr.ShowerStart();
	//auto const& dir = shr.Direction();
	TVector3 start  = shr.ShowerStart();
	TVector3 dir = shr.Direction();
	
	double this_mag = (start - single_shower_start).Mag();
	double closest_mag = 0.5;
	if( this_mag <= closest_mag){
		number_potential_shower_matches++;
	//	track_vertex_dir = current_vtx_dir;
		std::cout<<"potential shower match at index = "<<shr_index<<", track vertex = "<<start.X()<<", "<<start.Y()<<", "<<start.Z()<<std::endl;
		closest_mag = this_mag;
		matched_shower_ind = shr_index;
		matched_showers= true;
	}
}//find single photon shower match
std::cout<<"current matched shower index = "<<matched_shower_ind<<std::endl;


if (matched_showers ==true){
	auto const& shr = showerHandle->at(matched_shower_ind);
	TVector3 start  = shr.ShowerStart();
	TVector3 dir = shr.Direction();
	
//	if (start!=single_shower_start || dir!= single_shower_dir|| matched_showers ==false){
//		continue;
//	}
	std::cout<<"matched showers!"<<std::endl;
	std::cout<<"pandora shower start = "<<start.X()<<", "<<start.Y()<<", "<<start.Z()<<std::endl;
	//matched_showers=true;

	//store info about the single photon shower
	ShowerStart.push_back(start);

	auto this_len = shr.Length();
	ShowerLen.push_back(this_len);

	ShowerDir.push_back(dir);

	auto this_ang = shr.OpenAngle();
	ShowerAng.push_back(this_ang);
	std::cout<<"the shower direction X component = "<<dir.X()<<std::endl;
	
	//store info about the hits in the shower
	//for (size_t s = 0; s != my_shrs.size(); ++s ){
	//	auto const& current_shr = *(my_shrs.at(s));
	//	auto const& current_start = current_shr.ShowerStart();
	//	if (start == current_start){
	//		std::cout<<"matched showers"<<std::endl;
	std::cout<<"the shower length is  = "<<shr.Length()<<", and the opening angle is = "<<shr.OpenAngle()<<std::endl;
	_showerhitlist.clear();	
	//get the associated hits for the shower
			 //const std::vector<art::Ptr<recob::Hit> > shr_hit_v = shr_hit_assn_v.at(shr_index);
	auto shr_hit_v = shr_hit_assn_v.at(matched_shower_ind);
	fnum_hits_shower = shr_hit_v.size();
	std::cout<<"the number of hits in the shower = "<<fnum_hits_shower<<std::endl;		 
	int num_sshit_in_shower = 0;
			
	
	//check if each associated hit is an ssnet hit
	for (size_t h=0; h < shr_hit_v.size(); h++){
		auto const hit = *(shr_hit_v.at(h)); //get the hit from the pointer in the vector
		const recob::Hit* this_hit = &hit; //get the address to the hit that was stored in the vector
			
		//in the hit list, add corresponding shr hits to the shower list
		for(auto const& item : _hitlist){
			auto const stored_hit = (std::get<0>(item));
			if(matches(this_hit, stored_hit)==true){
						//std::cout<<"removing shower hit"<<std::endl;	
						//_hitlist.remove(item);
				_showerhitlist.push_back(item);
				num_sshit_in_shower++;
				break;
			}	
		}
	}//for each hit
	fnum_sshits_shower = num_sshit_in_shower;

		//	if (num_sshit_in_shower > 300){
		//		std::cout<<"warning, very large number of sshits"<<std::endl;
		//	}
			//std::cout<<"number of remaining matched shr hits = "<<_hitlist.size()<<std::endl;
	std::cout<<"the number of sshits in the shower = "<<fnum_sshits_shower<<std::endl;
	std::cout<<"the number of sshits in the shower list = "<<_showerhitlist.size()<<std::endl;
	if ((unsigned int)fnum_sshits_shower != _showerhitlist.size()){
		std::cout<<"warning, number if sshits in the shower do not match"<<std::endl;
	}
			
	fratio_hits_shower = (double)fnum_sshits_shower/fnum_hits_shower;
			//std::cout<<"the ratio of sshits to hits is "<<fratio_hits_shower<<std::endl;
	//		fselectTree->Fill();
	//}//if the showers match
	
	//}//for each shower from a 1 shower 1 track topology 

	//auto const& length = shr.Length();
	//std::cout<<"shower length = "<<length<<std::endl;
	
//	auto const& start  = shr.ShowerStart();
//	std::cout<<"shower start xyz = "<<start.X()<<", "<<start.Y()<<", "<< start.Z()<<std::endl;
}//for each shower
if (matched_showers ==false){
	std::cout<<"there was no matched shower"<<std::endl;
}

std::vector<TVector3> TrackStart;
std::vector<TVector3> TrackEnd;

bool matched_track = false;
int number_matched_trk_hits = 0;
int number_potential_matches = 0;
int matched_track_ind = -1; //the index of the track that best matches the single photon track
	
//loop over each track to find best match
  for ( size_t trk_index = 0; trk_index != trackHandle->size(); ++trk_index ){
	auto const& trk = trackHandle->at(trk_index);
	TVector3 current_vtx_dir = trk.VertexDirection();
		//std::cout<<"this track first point = "<<current_trk.FirstPoint()<<std::endl;
		//std::cout<<"this track last point = "<<current_trk.LastPoint()<<std::endl;
		//std::cout<<"this track has N points = "<<current_trk.NPoints()<<std::endl;
		
	double this_mag = (current_vtx_dir - single_track_vertex_dir).Mag();
	double closest_mag = 0.5;
	if( this_mag <= closest_mag){
		number_potential_matches++;
	//	track_vertex_dir = current_vtx_dir;
		std::cout<<"potential track match at index = "<<trk_index<<", track vertex = "<<current_vtx_dir.X()<<", "<<current_vtx_dir.Y()<<", "<<current_vtx_dir.Z()<<std::endl;
		closest_mag = this_mag;
		matched_track_ind = trk_index;
		matched_track= true;
	}
}//get matched track

std::cout<<"the best matched track index is currently "<<matched_track_ind<<std::endl;

//for the matched track, read values
auto const& trk = trackHandle->at(matched_track_ind);
_trackhitlist.clear();
TVector3 current_vtx = trk.Vertex();
auto current_end = trk.End();

TrackStart.push_back(current_vtx);
TrackEnd.push_back(current_end);


//if(number_potential_matches > 1){
//	std::cout<<"ERROR: multiple track matches"<<std::endl;
//}

if(matched_track_ind != -1){
//	if (current_vtx_dir == single_track_vertex_dir){
		std::cout<<"matched tracks!"<<std::endl;
//		matched_track = true;	
		
		//get the associated hits for the track
		const std::vector<art::Ptr<recob::Hit> > trk_hit_v = trk_hit_assn_v.at(matched_track_ind);
		fnum_hits_track = trk_hit_v.size();
		std::cout<<"the number of hits in the track = "<<fnum_hits_track<<std::endl;		 
		
		//for each hit
		for (size_t h=0; h < trk_hit_v.size(); h++){
			auto const hit = *(trk_hit_v.at(h));
			const recob::Hit* this_hit = &hit;
			
		//if hit is in the hit map, put it in the track hit map
		for (auto const& item : _hitlist){
			auto const stored_hit = (std::get<0>(item));
			if(matches(this_hit, stored_hit)== true){
				number_matched_trk_hits++;
				_trackhitlist.push_back(item);
				break;
			}
		}	

	}//for each hit
	fnum_sshits_track = number_matched_trk_hits;
	fratio_hits_track = (double)fnum_sshits_track/fnum_hits_track;
	fselectTree->Fill();	
	
	std::cout<<"the number of matched shr hits in the track = "<<number_matched_trk_hits<<std::endl;
	std::cout<<"the number of matched shr hits in the track list = "<<_trackhitlist.size()<<std::endl;
	if ((unsigned int)number_matched_trk_hits !=_trackhitlist.size()){
		std::cout<<"warning, number of hits in track don't match"<<std::endl;
	}	

}//if the tracks match
		
//}//for each track from a 1 shower 1 track topology 
//}//for each track
if(matched_track == false){
	std::cout<<"there was no matched track"<<std::endl;
}

/*
 * Look at remaining shrhits in ROI
 *
 */

//calc the radius for the ROI
//double radius = getRadius(fRadius);
//
//
	
//for each vertex position in the vector
//
//
//
//
//
//
//
//
//
//for ( size_t vtx_index = 0; vtx_index != my_vtxs.size(); ++vtx_index ){
	//std::cout<<"looking at vertex # "<<vtx_index<<std::endl;
//	double* xyz_pos = my_vtxs.at(vtx_index);
	//std::cout<<"the vertex pointer is "<<xyz_pos<<std::endl;
	
	//double X, Y, Z; 
//	double X = *(xyz_pos);
//	double Y = *(++xyz_pos);
//	double Z = *(++xyz_pos);
int ind_in_hitlist = 0;
	
//if the shower or track is missing, skip this event
if(matched_track == true && matched_showers== true){
	double X = single_vertex.X();
	double Y = single_vertex.Y();
	double Z = single_vertex.Z();

	std::cout<<"the vertex XYZ = "<<X<<", "<<Y<<", "<<Z<<std::endl;
	
	// subset of hit <-> sshit pairs within ROI of the vertex
	std::list<std::pair<const recob::Hit*, const recob::Hit* >>  _ROIhitlist;
	std::list<std::pair<const recob::Hit*, const recob::Hit* >> _ROIhitlist_plane0;
	std::list<std::pair<const recob::Hit*, const recob::Hit* >> _ROIhitlist_plane1;
	std::list<std::pair<const recob::Hit*, const recob::Hit* >> _ROIhitlist_plane2;

	//get the associated shower direction
 	int vtx_index = 0;
	//std::cout<<"matched showers = "<<matched_showers<<std::endl;
	//std::cout<<"flag1"<<std::endl;	
	//std::cout<<"ShowerStart length: "<<ShowerStart.size()<<std::endl;
	TVector3 shower_start = ShowerStart.at(vtx_index);
	double shower_length = ShowerLen.at(vtx_index);
	TVector3 shower_dir = ShowerDir.at(vtx_index);
	double shower_angle = ShowerAng.at(vtx_index);
	//std::cout<<"flag2"<<std::endl;

	//create shower end point
	TVector3 shower_end = findEnd(&shower_start, &shower_length, &shower_dir);

	//create 3d cone of rim
	double shower_radius = shower_length * shower_angle;
	std::cout<<"the shower radius is "<<shower_radius<<std::endl;
        std::vector<std::array<double, 3>> coneRim = Circle3D(shower_end, shower_dir, shower_radius);

	std::cout<<"The shower start is = "<<shower_start.X()<<", "<<shower_start.Y()<<", "<<shower_start.Z()<<std::endl;
	std::cout<<"The shower length is = "<<shower_length<<std::endl;
	std::cout<<"The calculated shower end is = "<<shower_end.X()<<", "<<shower_end.Y()<<", "<<shower_end.Z()<<std::endl;

	std::vector<TVector2> shower_dir_plane;
	std::vector<TVector2> vertex_wire_time_plane;
	std::vector<TVector2> dist_ang_plane0;
	std::vector<TVector2> dist_ang_plane1;
	std::vector<TVector2> dist_ang_plane2;

 	std::vector<int> ROIhits_ind_plane0;
	std::vector<int> ROIhits_ind_plane1;
	std::vector<int> ROIhits_ind_plane2;
	
	//std::vector<TVector2> vertex_time_plane;

	//get track start and end point
	TVector3 track_start = TrackStart.at(vtx_index);
	TVector3 track_end = TrackEnd.at(vtx_index);
	
	int n_theta_small_pos = 0;
 	int n_theta_small_neg = 0;
	
	std::vector<TVector2> hits_all_angles_plane0;
	std::vector<TVector2> hits_all_angles_plane1;
	std::vector<TVector2> hits_all_angles_plane2;


	//for each plane
	for (int plane = 0; plane <fPlanes; ++plane){	
		//std::cout<<"starting plane "<<plane<<std::endl;
	
		//get the wire and time for each plane
		//grabbed this from lareventdisplay/EventDisplay/RecoBaseDrawer.cxx
		//int plane = 0;
		//double wire = geo->WireCoordinate(Y, Z, plane, fTPC, fCryostat);
		//double time = detprop->ConvertXToTicks(X, plane, fTPC,fCryostat);
		//std::cout<<"the time and wire on plane "<<plane<<" is ("<<time<<", "<<wire<<")"<<std::endl;
		//
	
		
		double wire = calcWire(Y, Z, plane, fTPC, fCryostat, *fGeometry);
		double time = calcTime(X, plane, fTPC,fCryostat, *fDetprop);
		vertex_wire_time_plane.push_back(TVector2(wire, time));
		//vertex_time_plane.push_back(time);	
	
		std::cout<<"the time and wire on plane "<<plane<<" is ("<<time<<", "<<wire<<")"<<std::endl;
		out_stream<<"vertex "<<vtx_index<<" plane "<<plane<<" wire time "<<wire<<" "<<time<<"\n"; 		
		
		//get track projections for this plane
		TVector2 track_start_plane = TVector2(calcTime(track_start.X(), plane, fTPC,fCryostat, *fDetprop), calcWire(track_start.Y(), track_start.Z(), plane, fTPC, fCryostat, *fGeometry));	
		TVector2 track_end_plane =  TVector2(calcTime(track_end.X(), plane, fTPC,fCryostat, *fDetprop), calcWire(track_end.Y(), track_end.Z(), plane, fTPC, fCryostat, *fGeometry));

		out_stream<<"track "<<vtx_index<<" plane "<<plane<<" wire time "<<track_start_plane.Y()<<" "<<track_start_plane.X()<<" "<<track_end_plane.Y()<<" "<<track_end_plane.X()<<"\n";


		//get shower projections for this plane
		TVector2 shower_start_plane = TVector2(calcTime(shower_start.X(), plane, fTPC,fCryostat, *fDetprop), calcWire(shower_start.Y(), shower_start.Z(), plane, fTPC, fCryostat, *fGeometry));	
		TVector2 shower_end_plane =  TVector2(calcTime(shower_end.X(), plane, fTPC,fCryostat, *fDetprop), calcWire(shower_end.Y(), shower_end.Z(), plane, fTPC, fCryostat, *fGeometry));

		//std::vector<TVector2>* min_max;
		//std::vector<TVector2> shower_plane = *getMinMaxShowerPlane(coneRim,plane, fTPC, fCryostat, *fGeometry, *fDetprop, shower_start_plane, min_max);
		std::vector<TVector2> min_max;
		std::vector<TVector2> shower_plane = getMinMaxShowerPlane(coneRim,plane, fTPC, fCryostat, *fGeometry, *fDetprop, shower_start_plane, min_max);
	//std::vector<TVector2> shower_plane = min_max;
//		TVector2 dist_min = shower_plane.at(1);
//		TVector2 dist_max = shower_plane.at(0);

//		out_stream<<"shower "<<vtx_index<<" plane "<<plane<<" startwt "<<shower_start_plane.Y()<<" "<<shower_start_plane.X()<<" "<<dist_min.Y()<<" "<<dist_min.X()<<" "<<dist_max.Y()<<" "<<dist_max.X()<<" \n";

		std::string my_out;
		for (auto const& item : shower_plane){
			//std::cout<<"this address = "<<&item<<std::endl;
			//std::cout<<"this item "<<item.Y()<<" "<<item.X()<<std::endl;
			my_out+= " " + std::to_string(item.Y()) + " " +  std::to_string(item.X()) + " ";
			//std::cout<<"my_out = "<<my_out<<std::endl;
		}				
		//std::cout<<"my_out = "<<my_out<<std::endl;

		out_stream<<"shower "<<vtx_index<<" plane "<<plane<<" startwt "<<shower_start_plane.Y()<<" "<<shower_start_plane.X()<<""<<my_out<<"\n";

		out_stream<<"ROI "<<vtx_index<<" plane "<<plane;
		std::vector<TVector2> points_ROI_this_plane = drawROI(fTimeToCMConstant, fWireToCMConstant,fRadius, time,wire);
		for (auto const& edge_point: points_ROI_this_plane){
			//std::cout<<"this 2d point on the ROI edge = "<<edge_point.X()<<" "<<edge_point.Y()<<std::endl;	
			out_stream<<" "<<edge_point.Y()<<" "<<edge_point.X();
		}
		out_stream<<"\n";
		
		TVector2 shower_direction_plane = calcNormVec(shower_start_plane, shower_end_plane); 
		shower_dir_plane.push_back(shower_direction_plane);
		
		//std::cout<<"the direction of the shower on this plane is "<< shower_direction_plane.X() <<" starting and ending from "<< shower_start_plane.X()<<shower_end_plane.X() <<std::endl;
		//for each remaining sshit
//		int ind_in_hitlist = 0;
		for(auto const& item : _hitlist){
			//if the sshit falls inside the ROI and on the plane
			auto const this_sshit = std::get<1>(item);
			int this_plane =  (*this_sshit).View();
			if(this_plane != plane){continue;}

			double dist = distToVtx(fTimeToCMConstant, fWireToCMConstant, fRadius, time, wire, this_sshit);
			//if (dist<100){
			//	std::cout<<"the dist of this hit to the vtx is "<<dist<<std::endl;
			//}
			if (inROI(fRadius, dist) == true){
			//	if(_ROIhitlist.size() < 3){
			//		std::cout<<"corresponding shower dir x = "<<shower_dir.X()<<std::endl;
			//	}
				_ROIhitlist.push_back(item);
                                fradial_dist_sshit_vtx = dist;
				fopening_angle_shower_sshit = angleFromShower(shower_direction_plane, this_sshit, time,wire);
				//std::cout<<"the angle from the shower is "<<angle<<std::endl;
				if (fradial_dist_sshit_vtx>fRadius){
					std::cout<<"ERROR: hit in ROIHitList outside of ROI, radial dist = "<< fradial_dist_sshit_vtx<<std::endl;
				}
				if ( fopening_angle_shower_sshit < 0.2){
					n_theta_small_pos++;
				}
				if ( fopening_angle_shower_sshit > 2.95){
					n_theta_small_neg++;
					//std::cout<<"the angle between the hit and the shower is "<<opening_angle<<std::endl;
				}

				if (plane==0){
					dist_ang_plane0.push_back(TVector2(fradial_dist_sshit_vtx,fopening_angle_shower_sshit));
					ROIhits_ind_plane0.push_back(ind_in_hitlist);
					_ROIhitlist_plane0.push_back(item);
				}	

				if (plane==1){
					dist_ang_plane1.push_back(TVector2(fradial_dist_sshit_vtx,fopening_angle_shower_sshit));
					ROIhits_ind_plane1.push_back(ind_in_hitlist);
					_ROIhitlist_plane1.push_back(item);
				}	

				if (plane==2){
					dist_ang_plane2.push_back(TVector2(fradial_dist_sshit_vtx,fopening_angle_shower_sshit));
					ROIhits_ind_plane2.push_back(ind_in_hitlist);
					 _ROIhitlist_plane2.push_back(item);
				}	
	
				//fROITree->Fill();	
			
			ind_in_hitlist++;
			} //if in ROI
	//	ind_in_hitlist++;
		}//each sshit

		//std::cout<<"ending plane "<<plane<<std::endl;
		if(plane==0){
		//	int n_inc = 100;
			hits_all_angles_plane0 =  getNumHitsSectors(n_inc,shower_end_plane, fRadius, time, wire, _ROIhitlist_plane0, fTimeToCMConstant,fWireToCMConstant);
		}
		if(plane==1){
			hits_all_angles_plane1 =  getNumHitsSectors(n_inc, shower_end_plane, fRadius, time, wire, _ROIhitlist_plane1, fTimeToCMConstant,fWireToCMConstant);
		}
		if(plane==2){
			hits_all_angles_plane2 =  getNumHitsSectors(n_inc, shower_end_plane, fRadius, time, wire, _ROIhitlist_plane2, fTimeToCMConstant,fWireToCMConstant);
		}



	}//for each plane
std::cout<<"the number of shower hits in the ROI on plane 0 = "<<ROIhits_ind_plane0.size()<<std::endl;
std::cout<<"the number of shower hits in the ROI on plane 1 = "<<ROIhits_ind_plane1.size()<<std::endl;
std::cout<<"the number of shower hits in the ROI on plane 2 = "<<ROIhits_ind_plane2.size()<<std::endl;

std::cout<<"the number of hits at small positive angles = "<<n_theta_small_pos<<" and the number for small negative = "<<n_theta_small_neg<<std::endl;

/*
for (auto item: ROIhits_ind_plane0){
	std::cout<<"plane 0 index = "<<item<<std::endl;
}
for (auto item: ROIhits_ind_plane1){
        std::cout<<"plane 1 index = "<<item<<std::endl;
}
for (auto item: ROIhits_ind_plane2){
        std::cout<<"plane 2 index = "<<item<<std::endl;
}
*/

/*
 *
 * Subtract out hits from the shower and track from the ROI
 *
 * */

//make list of sshits in ROI no track
//auto const _listnotrack = removeHitsROIList(_ROIhitlist, _trackhitlist);
auto const _listindnotrack = getIndHitsROIList(_ROIhitlist, _trackhitlist);
std::cout<<"the number of track hits in the ROI is "<<_listindnotrack.size()<<std::endl;
//auto dist_ang_notrack =  fillROITree( _listnotrack, vertex_wire_time_plane, shower_dir_plane, fPlanes, fTimeToCMConstant,  fWireToCMConstant, fRadius);
//for(auto const& item : dist_ang_notrack){
//	fradial_dist_sshit_vtx_notrack = item.X();
//	fopening_angle_shower_sshit_notrack = item.Y();
	//std::cout<< fradial_dist_sshit_vtx_notrack<<", "<<fopening_angle_shower_sshit_notrack<<std::endl;
	//fROITree->Fill();
//}


//make list of sshits in ROI no shower
//auto const _listnoshower = removeHitsROIList(_ROIhitlist, _showerhitlist);
auto const _listnoshower = getIndHitsROIList(_ROIhitlist, _showerhitlist);
std::cout<<"the number of shower hits in the ROI is "<<_listnoshower.size()<<std::endl;
//auto dist_ang_noshower =  fillROITree( _listnoshower, vertex_wire_time_plane, shower_dir_plane, fPlanes, fTimeToCMConstant,  fWireToCMConstant, fRadius);
//for (auto const& item : dist_ang_noshower){
//	if (item.X()>fRadius){
//		std::cout<<"ERROR before filling tree, ssnet hit in ROI minus shower outside radius, distance = "<<item.X()<<std::endl;
//	}
//}
	
//make list of sshits in ROI no track or shower
//auto const _listnoshowernotrack = getIndHitsROIList(_listnoshower, _trackhitlist);
std::vector<int> _listnoshowernotrack;
//add all hits in ROI minus the track
//for (auto const& item : _listindnotrack){
//	_listnoshowernotrack.push_back(item);
//}

//for each hit in the ROI minus the shower
//int ind = 0;
std::vector<int> indices_to_remove;
for (auto const& item : _listnoshower){
	bool contains = false;
	//for each hit in the ROI minus the track
	for (auto const& item2 : _listindnotrack){
		//if both the no shower and no track lists contains the same hit
		if (item == item2){
        		contains = true;
        		//remove the hit
			//_listnoshowernotrack.remove(item);
		}
	}
	//if the no track list doesn't contain the hit from the no shower list, remove it from the no shower no track list
	if(contains == true){ 
		_listnoshowernotrack.push_back(item);
       		//_listnoshowernotrack.erase(ind);
	}
	//ind++;
}
//std::cout<<"the number of shower hits to remove from the no track list is "<<indices_to_remove.size()<<std::endl;
//for (int ind : _listnoshower) {
//	_listnoshowernotrack.erase( _listnoshowernotrack.begin() + ind);
//}

//std::cout<<"the number of shower and track hits in the ROI is "<<_listnoshowernotrack.size()<<std::endl;
//auto dist_ang_noshower =  fillROITree( _listnoshower, vertex_wire_time_plane, shower_dir_plane, fPlanes, fTimeToCMConstant,  fWireToCMConstant, fRadius);
//auto dist_ang_noshowernotrack =  fillROITree( _listnoshowernotrack, vertex_wire_time_plane, shower_dir_plane, fPlanes, fTimeToCMConstant,  fWireToCMConstant, fRadius);

	
fnum_sshits_ROI = _ROIhitlist.size();
fnum_sshits_ROI_no_shower = _listnoshower.size();
fnum_sshits_ROI_no_track = _listindnotrack.size();
fnum_sshits_ROI_no_track_no_shower = _listnoshowernotrack.size();
//fnum_sshits_ROI_no_track_no_shower = (-1*fnum_sshits_ROI ) + (fnum_sshits_ROI_no_track + fnum_sshits_ROI_no_shower);

std::cout<<"num hits in ROI, -shower, -track, -shower and track"<<fnum_sshits_ROI<<", "<<fnum_sshits_ROI_no_shower<<", "<<fnum_sshits_ROI_no_track<<", "<<fnum_sshits_ROI_no_track_no_shower<<std::endl;

/*
 *
 *Write the ROI info to the tree and histograms 
 *
 * */


bool make_hists_this_event= false;
if (fRun == fRun_hist && fSubRun == fSubrun_hist && fEvent == fEvent_hist){
	std::cout<<"Making individual histograms for this event"<<std::endl;
	make_hists_this_event =true;	
}
std::cout<<make_hists_this_event<<std::endl;

int n = 0;
int n_plane0 = 0;
int n_plane1 = 0;
int n_plane2 = 0;
int n_shower_plane0 = 0;
//std::cout<<"the number of hits on plane 0 is "<<dist_ang_plane0.size()<<std::endl; 
//for(size_t n = 0; n != _ROIhitlist.size(); ++n){
for(auto const& item : _ROIhitlist){
	auto const this_sshit = std::get<1>(item);
	int plane =(*this_sshit).View();
	if (make_hists_this_event ==true){
//		std::cout<<"this sshit plane, time, wire = "<<plane<<", "<<(*this_sshit).PeakTime()<<", "<<(*this_sshit).WireID().Wire<<std::endl;
		//std::cout<<"the radial and angular dist for this sshit are = "<<dist_ang_plane0.at(n).X()<<", "<<dist_ang_plane0.at(n).Y()<<std::endl;
		}


	if(plane ==0){
		int ind_in_plane0 = contains_at_ind(n,ROIhits_ind_plane0);
		n_plane0++;	
		if(ind_in_plane0!=-1){
			fradial_dist_sshit_vtx = dist_ang_plane0.at(ind_in_plane0).X();
    			fopening_angle_shower_sshit = dist_ang_plane0.at(ind_in_plane0).Y();
//			std::cout<<"sshit on plane 0 at index (ind_in_plane0)"<<ind_in_plane0<<std::endl;
		//	std::cout<<"the dist and ang for this index is "<<fradial_dist_sshit_vtx<<", "<<fopening_angle_shower_sshit<<std::endl;
		}
	} 
	if(plane ==1){
		int ind_in_plane1 = contains_at_ind(n,ROIhits_ind_plane1);	
		//std::cout<<"On plane 1, the ROI ind is "<<n<<std::endl;
		n_plane1++;
		if(ind_in_plane1 != -1){
			fradial_dist_sshit_vtx = dist_ang_plane1.at(ind_in_plane1).X();
    			fopening_angle_shower_sshit = dist_ang_plane1.at(ind_in_plane1).Y();
//			std::cout<<"sshit on plane 1 at index (ind_in_plane1)"<<ind_in_plane1<<std::endl;
		}
	} 

	if(plane ==2){
		int ind_in_plane2 = contains_at_ind(n,ROIhits_ind_plane2);
		n_plane2++;
		if(contains_at_ind(n,ROIhits_ind_plane2) != -1){
			fradial_dist_sshit_vtx = dist_ang_plane2.at(ind_in_plane2).X();
    			fopening_angle_shower_sshit = dist_ang_plane2.at(ind_in_plane2).Y();
//			std::cout<<"sshit on plane 2 at index (ind_in_plane2)"<<ind_in_plane2<<std::endl;
		}
	} 
	
	if(make_hists_this_event ==true){
		//std::cout<<"the radial dist for this ssnet hit is "<<fradial_dist_sshit_vtx<<std::endl;
		fDist_sshits->Fill(fradial_dist_sshit_vtx);	
		fAngle_sshits->Fill(fopening_angle_shower_sshit);
		//std::cout<<"the dist and ang for this index is "<<fradial_dist_sshit_vtx<<", "<<fopening_angle_shower_sshit<<std::endl;
		if(plane == 0){
			fDist_sshits_plane0->Fill(fradial_dist_sshit_vtx);
                	//fAngle_sshits_plane0->Fill(fopening_angle_shower_sshit);
			for(unsigned int i= 0; i < hits_all_angles_plane0.size(); i++){
				fAngle_sshits_plane0->SetBinContent(i, hits_all_angles_plane0.at(i).Y()); 
			}
		//	for(auto item: hits_all_angles_plane0){
				//int bin = fAngle_sshits_plane0->GetBin(180*item.X()/M_PI);
		//		fAngle_sshits_plane0->SetBinContent(bin, item.Y());
		//	}

		}

		if(plane == 1){
			fDist_sshits_plane1->Fill(fradial_dist_sshit_vtx);
                	//fAngle_sshits_plane1->Fill(fopening_angle_shower_sshit);
			for(unsigned int i= 0; i < hits_all_angles_plane1.size(); i++){
                                fAngle_sshits_plane1->SetBinContent(i, hits_all_angles_plane1.at(i).Y());
                        }
		}
	
		if(plane == 2){
			fDist_sshits_plane2->Fill(fradial_dist_sshit_vtx);
                	//fAngle_sshits_plane2->Fill(fopening_angle_shower_sshit);
			for(unsigned int i= 0; i < hits_all_angles_plane2.size(); i++){
                                fAngle_sshits_plane2->SetBinContent(i, hits_all_angles_plane2.at(i).Y());
                        }
		}

		//std::cout<<"filling the histograms"<<std::endl;	
		//std::cout<<"for selected event, radial dist and angle are "<<fradial_dist_sshit_vtx<<", "<<fopening_angle_shower_sshit<<std::endl;
	}
	if(contains(n, _listnoshower)==true){
	 	//fradial_dist_sshit_vtx_noshower = dist_ang.at(n).X();
                //fopening_angle_shower_sshit_noshower = dist_ang.at(n).Y();
	 	fradial_dist_sshit_vtx_noshower = fradial_dist_sshit_vtx;
                fopening_angle_shower_sshit_noshower = fopening_angle_shower_sshit;
		if(make_hists_this_event ==true){	
			//std::cout<<"ind for no shower is "<<n<<std::endl;
   			//std::cout<<"filling no shower for specified event, plane = "<<plane<<std::endl;
	             	fDist_sshits_noshower->Fill(fradial_dist_sshit_vtx_noshower);
                	fAngle_sshits_noshower->Fill(fopening_angle_shower_sshit_noshower);
			if(plane ==0){
				//std::cout<<"filling no shower plane 0 at ind "<<n<<std::endl;
				fDist_sshits_noshower_plane0->Fill(fradial_dist_sshit_vtx_noshower);
	                        fAngle_sshits_noshower_plane0->Fill(fopening_angle_shower_sshit_noshower);		
			}

			if(plane ==1){
				fDist_sshits_noshower_plane1->Fill(fradial_dist_sshit_vtx_noshower);
	                        fAngle_sshits_noshower_plane1->Fill(fopening_angle_shower_sshit_noshower);		
			}

			if(plane ==2){
				fDist_sshits_noshower_plane2->Fill(fradial_dist_sshit_vtx_noshower);
	                        fAngle_sshits_noshower_plane2->Fill(fopening_angle_shower_sshit_noshower);		
			}

        	}
	}
	
	 if(contains(n, _listnoshower)==false){
		if(plane ==0){
			//std::cout<<"removing shower hit on plane 0 "<<n_shower_plane0<<"at ind"<<n<<std::endl;
			n_shower_plane0++;
		}
	}

 //   	if( n<unsigned(fnum_sshits_ROI_no_shower)){
//		fradial_dist_sshit_vtx_noshower = dist_ang_noshower.at(n).X(); 
//    		fopening_angle_shower_sshit_noshower = dist_ang_noshower.at(n).Y(); 
		//if (fradial_dist_sshit_vtx_noshower>fRadius){
		//	std::cout<<"ERROR: when filling tree ssnet hit in ROI minus shower outside radius, distance = "<<fradial_dist_sshit_vtx_noshower<<std::endl;
		//}
//	}
	if(contains(n, _listindnotrack)==true){    
		//fradial_dist_sshit_vtx_notrack = dist_ang_notrack.at(n).X(); 
     		//fopening_angle_shower_sshit_notrack = dist_ang_notrack.at(n).Y(); 
	        fradial_dist_sshit_vtx_notrack = fradial_dist_sshit_vtx; 
     		fopening_angle_shower_sshit_notrack = fopening_angle_shower_sshit; 
	} 
	
//	if(contains(n, _listnoshower)==true || contains(n, _listindnotrack)==true){
 //               fradial_dist_sshit_vtx_noshowernotrack = dist_ang.at(n).X();
  //              fopening_angle_shower_sshit_noshowernotrack = dist_ang.at(n).Y();
   //     }


	if(contains(n, _listnoshowernotrack)==true){   
		fradial_dist_sshit_vtx_noshowernotrack = fradial_dist_sshit_vtx; 
    		fopening_angle_shower_sshit_noshowernotrack =  fopening_angle_shower_sshit;	
	}
	n++;
	fROITree->Fill();
}

/**
 *
 *Run some fits 
 *
 **/

//fit all the ssnet hits on plane 0 with two gaussians
//F(mean1, sigma1, mean2,sigma2) = gauss(mean1,sigma1)+ gauss(mean2,sigma2)
//

//if(make_hists_this_event ==true){
//	TF1* f1= new TF1("f1","gaus([0], [1], [2]) + gaus([3], [4], [5])", angle_min, angle_max);
//	fAngle_sshits_plane0->Fit(f1);
//}

//double fitFunction(double *x, double *par){
//      return gaus(x,par) + gaus(x,&par[2]);
//}
if(make_hists_this_event ==true){

	//TF1 *func = new TF1("fit",fitf,angle_min,angle_max,6);
	//func->SetParNames("Norm1","Mean1","sigma1","Norm2","Mean2","Sigma2");
	//func->SetParLimits(1,-40,40); //the first peak should be close to 0
	//func->SetParLimits(2,-50,50); //the width shouldn't be too big
	//func->SetParLimits(4,-50,50); //the width shouldn't be too big

	TF1* func= new TF1("fitf",fitf,angle_min,angle_max,6);
	//func->SetParNames("Norm1","Mean1","sigma1");
	func->SetParNames("Norm1","Mean1","sigma1","Norm2","Mean2","Sigma2");
	func->SetParLimits(0, 0, 100);
	func->SetParLimits(1,-30,30); //the first peak should be close to 0
        func->SetParLimits(2,-30,30); //the width should be <50
	func->SetParameter(0, 50.0);
	func->SetParameter(1, 0.01);
	func->SetParameter(2, 10.0);
	func->SetParLimits(3, 0, 100);
	func->SetParLimits(4,angle_min,angle_max); //the first peak should be close to 0
        func->SetParLimits(5,-30,30); //the width should be <50
	func->SetParameter(3, 50.0);
	func->SetParameter(4, 30);
	func->SetParameter(5, 10.0);
	

	TFitResultPtr r0 = fAngle_sshits_plane0->Fit("fitf","MSVR");

 	TFitResultPtr r1 = fAngle_sshits_plane1->Fit("fitf","MSVR");	
	TFitResultPtr r2 = fAngle_sshits_plane2->Fit("fitf","MSVR");	
}
//double chi2   = r->Chi2();                  // to retrieve the fit chi2
//double par0   = r->Parameter(0);            // retrieve the value for the parameter 0
//double err0   = r->ParError(0);             // retrieve the error for the parameter 0

//std::cout<<"the chi2 for the fit is "<<chi2<<std::endl;

//store the chi-square for each event 


std::cout<<"the number of hits on plane0-2 when filling the tree was "<<n_plane0<<", "<<n_plane1<<", "<<n_plane2<<std::endl;

//if(make_hists_this_event ==true){
  // fDist_sshits_noshower->Draw();
  // fAngle_sshits_noshower->Draw();
// }


//my_hist->Fill(fnum_sshits_ROI_no_track_no_shower);
std::cout<<"the number of shower hits within the ROI before removing track or shower =  "<<fnum_sshits_ROI<<std::endl;

//if ((unsigned int)fnum_sshits_ROI != _listnotrack.size()){
//	std::cout<<"warning, some track sshits removed from ROI"<<std::endl;
//}
fselectTree->Fill();

//}//loop over vertices
}//if there is a matched track and shower
else{
	std::cout<<"ERROR: skipping this event, missing track and/or shower"<<std::endl;
}

//total_num_vertices+= my_vtxs.size();
total_num_vertices++;
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
 }//if there's a matching single photon event
 } // SSNetTest::analyze()
 
void SSNetTest::endJob(){
	//fb.close();
	out_stream.close();
	in_stream.close();
	std::cout<<"Ending job"<<std::endl;		
}
  
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
	double dist_from_vtx = sqrt((diff_t*diff_t) + (diff_w*diff_w));
	return dist_from_vtx;
}//distToVtx

//double distToVtx2d(TVector2 my_point, TVector2 vertex){
//	TVector2 dist = vertex-my_point;
//	return dist.Mod();
//}

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

//for a given vertex, ROI radius length, and plane, returns a set of 2D points outlining the ROI
std::vector<TVector2> drawROI(double fTimeToCMConstant, double fWireToCMConstant, double radius, double vertex_time, double vertex_wire){
	//constexpr unsigned short nRimPts = 16;
	//std::vector<TVector2, (nRimPts + 1)> rimPts;
	double rad_time = fTimeToCMConstant*radius;
	double rad_wire = fWireToCMConstant* radius;
	std::vector<TVector2> points = drawEllipse(vertex_wire, vertex_time, rad_wire, rad_time);
	return points;
}//draw ROI


std::vector<TVector2> drawEllipse(double start_x, double start_y, double ax_x, double ax_y){
	int nRimPts = 16;
	std::vector<TVector2> my_pts; 
	double x_inc = 2*ax_x/nRimPts;
	//std::cout<<"a= "<<ax_x<<", b= "<<ax_y<<std::endl;
	//std::cout<<"the x inrement is "<<x_inc<<std::endl;
	double coeff = ax_y/ax_x;
	//std::cout<<"the b/a coeff where a>b = "<<coeff<<std::endl;
	for (int i = -nRimPts/2; i <= nRimPts/2; ++i){
		double this_x = start_x + i*x_inc;
		double x_minus_start = this_x - start_x;
		double sqrt_term1 = pow(ax_y,2);
		double sqrt_term2 = pow(coeff, 2)*pow(x_minus_start,2);
		//std::cout<<"the first term under the radical = "<<sqrt_term1<<" and the second term under the radical = "<<sqrt_term2<<std::endl;
		double term2 = sqrt(sqrt_term1 - sqrt_term2);
		double y_plus = start_y + term2;
		double y_minus= start_y - term2;
		//std::cout<<"term 1 = "<<this_x<<" and term 2 = "<<term2<<std::endl;
		my_pts.push_back(TVector2(this_x, y_plus));
		my_pts.push_back(TVector2(this_x, y_minus));
	}
	return my_pts;
}//draw ellipse

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
 	TVector2 norm_hit_vtx_direction_plane  = hit_vtx_direction_plane.Unit();
        TVector2 norm_shower_direction_plane = shower_direction_plane.Unit();
	double dot =  (norm_hit_vtx_direction_plane.X()*norm_shower_direction_plane.X()) + (norm_hit_vtx_direction_plane.Y()*norm_shower_direction_plane.Y());
        //std::cout<<"the dot product between the hit vector and the shower vector is "<< dot<<std::endl;
	        								
        double opening_angle = acos(dot);
	if (norm_hit_vtx_direction_plane.X() < norm_shower_direction_plane.X()){
        	opening_angle = -opening_angle;
	}
	return opening_angle;
}

//takes a vector of (time, wire) and returns a vector of (cm, cm)
TVector2 calcCM(TVector2 vec_time_wire, double fTimeToCMConstant, double fWireToCMConstant){
	return TVector2(vec_time_wire.X()* fTimeToCMConstant, vec_time_wire.Y()*fWireToCMConstant);
}

//returns true if v2 is clockwise of v1 in a circle
bool areClockwise(TVector2 v1, TVector2 v2){
	if ( -v1.X()*v2.Y() + v1.Y()*v2.X() > 0){
		return true;
	} else{
		return false;
	}
}

//for a given plane, takes vectors marking the boundaries of a sector in the ROI (for fixed angular incements)
//if the hit is in the sector, returns true
//vectors are in cm, vertex info is is time+wire
bool inSector(TVector2 sector_start_cm, TVector2 sector_end_cm,  const recob::Hit*  hit, double vertex_time, double vertex_wire, double fTimetoCMConstant, double fWiretoCMConstant){
 	double hit_time = (*hit).PeakTime();
        double hit_wire = (*hit).WireID().Wire;
	
	double rel_time = hit_time - vertex_time;
	double rel_wire = hit_wire - vertex_wire;
	
	TVector2 rel_point = TVector2(rel_time, rel_wire);

	//convert the time and wire to CM
	TVector2 rel_point_cm = calcCM(rel_point, fTimetoCMConstant, fWiretoCMConstant);
	//TVector2 sector_start_cm  = calcCM(sector_start);
	//TVector2 sector_end_cm  = calcCM(sector_end);

	//if this point is clockwise from the end and counterclockwise from the start, it's in the sector
	if (!areClockwise(sector_start_cm, rel_point_cm) && areClockwise(sector_end_cm, rel_point_cm) ){
		return true;
	} else {
		return false;
	}
}


//given a start vector direction (given the shower direction), a vertex position 
//iterates 360 over a circle for a given number of increments
//fills a vector with the number of hits in each sector
std::vector<TVector2> getNumHitsSectors(int n_inc,TVector2 shower_end, double radius, double vertex_time, double vertex_wire, std::list<std::pair<const recob::Hit*, const recob::Hit* >> _ROIlist, double fTimetoCMConstant, double fWiretoCMConstant){

	std::vector<TVector2> hits_per_bin;
	double theta_seg = 2*M_PI/n_inc;

	std::vector<int> hit_has_been_filled(_ROIlist.size(), 0);
	
	TVector2 shower_end_cm =  calcCM(shower_end, fTimetoCMConstant, fWiretoCMConstant);
	TVector2 vertex = calcCM(TVector2(vertex_time, vertex_wire), fTimetoCMConstant, fWiretoCMConstant);	
	std::cout<<"the vertex on this plane is "<<vertex.X()<<", "<<vertex.Y()<<std::endl;
	TVector2 sec_start = shower_end_cm - vertex;
	sec_start = sec_start * -1;
	std::cout<<"the sector start for sector 0 is "<<sec_start.X()<<", "<<sec_start.Y()<<std::endl;
	//TVector2 sec_start = vertex + TVector2(0, radius);//define the 0 point wrt the vertex in CM 
	TVector2 sec_end;
	//Tvector2 sec_end = getSectorEnd(sec_start, theta_seg);
	//calc number in the first sector

	double this_theta = 0;

	//push back number
	//hits_per_bin[0] = n_sec;
	int n_hits_ROI_plane = 0;
	//for each sector
	for(int i = 0; i<n_inc; i++){
		int n_sec = 0; //the number of hits in currect sector
		//set the start vector(using previous if n!= 0)
		if (i > 0){
			sec_start = sec_end;
		}
		//sec_start = sec_end; //the new start is the end to the previous sector
		//calc the end vector
		sec_end = getSectorEnd(sec_start, theta_seg);
		
		//std::cout<<"the start vector is "<<sec_start.X()<<", "<<sec_start.Y()<<"for sector "<<i<<std::endl;
		//std::cout<<"the end vector is "<<sec_end.X()<<", "<<sec_end.Y()<<std::endl;		
		//for each hit in the list
		int ind_in_hitlist = 0;
		for(auto const& item : _ROIlist){
			auto const this_hit = (std::get<0>(item));
			
			//if hit is in sector
			if(inSector(sec_start, sec_end, this_hit, vertex_time, vertex_wire,  fTimetoCMConstant, fWiretoCMConstant) ==true){
				n_sec++; //increment number of hits in this sector
				hit_has_been_filled[ind_in_hitlist]++;				
			}
			ind_in_hitlist++;
		}//for each hit
		//push back the number of hits? add to a tree? some form of storage
		hits_per_bin.push_back(TVector2(this_theta, n_sec));
		n_hits_ROI_plane += n_sec;
		this_theta += theta_seg;
		//std::cout<<"The angle for this sector is "<<this_theta<<std::endl;
		//std::cout<<"The number of hits in sector "<<i<<" is "<<n_sec<<std::endl;
	}//for each sector
	std::cout<<"the total hits in the ROI in this plane = "<<n_hits_ROI_plane<<std::endl;

	for (auto item: hit_has_been_filled ){
		//std::cout<<"this his was read "<<item<<" times"<<std::endl;
		if(item != 1){
			std::cout<<"ERROR, this hit was counted "<<item<<" times"<<std::endl;
		}
	}
	return hits_per_bin;
}

//for a given vector (in cm!!!) marking the start of a sector and a given theta value, calculates the end of the sector assuming the end
//is counterclockwise of the start
//centered at the vertex on the plane
TVector2 getSectorEnd(TVector2 sector_start_cm, double theta){
	//TVector2 sector_end_cm;
	double end_time;
	double end_wire;

	//TVector2 vertex = TVector2(vertex_time, vertex_wire);
	//TVector2 vertex_cm = calcCM(vertex, fTimeToCMConstant, fWireToCMConstant);

	//do a 2d rotation	
	end_time = sector_start_cm.X() * cos(theta) - sector_start_cm.Y() * sin(theta);
	end_wire = sector_start_cm.X() * sin(theta) + sector_start_cm.Y() * cos(theta);

	return TVector2(end_time, end_wire);
}
/*
//takes two lists, one of SShits in the ROI and one of all shower or track hits in the event
//returns a new list of the ROI hits with the ones included in the track/shower removed
std::list<std::pair<const recob::Hit*, const recob::Hit* >> removeHitsROIList(std::list<std::pair<const recob::Hit*, const recob::Hit* >> _ROIlist, std::list<std::pair<const recob::Hit*, const recob::Hit* >> _showerlist){
	int num_hits_obj_in_ROI = 0;
	//make an empty list
	std::list<std::pair<const recob::Hit*, const recob::Hit* >> _listcopy;
	
	//for each hit in the ROI
	for(auto const& item : _ROIlist){
		//copy it to the new list
		_listcopy.push_back(item);
		auto const this_hit = (std::get<0>(item));
	
		//for each hit in the shower/track
		for(auto const& this_item : _showerlist){
			 auto const stored_hit = (std::get<0>(this_item));
                         //if the shower/track hit matches a hit in the ROI
			 if(matches(this_hit, stored_hit)== true){
				//remove the hit from the list
				 _listcopy.remove(this_item);
				num_hits_obj_in_ROI++;
			}
		}
	}
	std::cout<<"the number of hits from the track/shower object in the ROI is "<<num_hits_obj_in_ROI<<std::endl;
	return _listcopy;
}
*/

//takes two lists, one of SShits in the ROI and one of all shower or track hits in the event
//returns a list of indices in the ROI list corresponding to the hits in the ROI minus those in the track/shower 
std::vector<int> getIndHitsROIList(std::list<std::pair<const recob::Hit*, const recob::Hit* >> _ROIlist, std::list<std::pair<const recob::Hit*, const recob::Hit* >> _showerlist){
	int num_hits_obj_in_ROI = 0;
	//make an empty list
	std::vector<int> indices;	

	int this_ind = 0;
	//for each hit in the ROI
	for(auto const& item : _ROIlist){
		auto const this_hit = (std::get<0>(item));
		bool match = false;	
		//for each hit in the shower/track
		for(auto const& this_item : _showerlist){
			 auto const stored_hit = (std::get<0>(this_item));
                         //if the shower/track hit matches a hit in the ROI
			 if(matches(this_hit, stored_hit)== true){
				num_hits_obj_in_ROI++;
				match = true;
			//}else{//store the index of the ROI hits that don't match with the track/shower
			}
		}
		//if there is no match for this hit in the shower/track, save the index
		if (match==false){
			indices.push_back(this_ind);
		}
		this_ind++;
	}
	std::cout<<"the number of hits from the track/shower object in the ROI is "<<num_hits_obj_in_ROI<<std::endl;
	return indices;
}


//draw a circle in 3D given the length, direction, and  radius
std::vector<std::array<double, 3>> Circle3D(const TVector3& centerPos, const TVector3& axisDir, const double& radius)  {
   // B. Baller Create a polyline circle in 3D 
      
      // Make the rotation matrix to transform into the circle coordinate system
     TRotation r;
     r.RotateX(axisDir.X());
     r.RotateY(axisDir.Y());
     r.RotateZ(axisDir.Z());
     constexpr unsigned short nRimPts = 16;
     std::vector<std::array<double, 3>> rimPts(nRimPts + 1);
     for(unsigned short iang = 0; iang < nRimPts; ++iang) {
       double rimAngle = iang * 2 * M_PI / (float)nRimPts;
        TVector3 rim = {0, 0, 1};
        rim.SetX(radius * cos(rimAngle));
        rim.SetY(radius * sin(rimAngle));
        rim.SetZ(0);
        rim.Transform(r);
        rim += centerPos;
        for(unsigned short ixyz = 0; ixyz < 3; ++ixyz) rimPts[iang][ixyz] = rim[ixyz];
      } // iang
     // close the circle
     rimPts[nRimPts] = rimPts[0];
     return rimPts;
}
/*
std::vector<TVector2> fillROITree(std::list<std::pair<const recob::Hit*, const recob::Hit* >> _listnotrack, std::vector<TVector2>vertex_wire_time_plane, std::vector<TVector2> shower_dir_plane, int fPlanes, double fTimeToCMConstant, double fWireToCMConstant, double fRadius){
	std::vector<TVector2> dist_ang;
	for (auto const& item : _listnotrack){
		auto const this_sshit = std::get<1>(item);
		for(int plane = 0; plane <fPlanes; ++plane){
			if ((*this_sshit).View()!= plane){
				continue;
			}	
			auto wire =(vertex_wire_time_plane.at(plane).X());
			auto time =(vertex_wire_time_plane.at(plane).Y());
			//std::cout<<"the plane is "<<plane<<" and the vertex position is (wire, time) = "<<wire<<", "<<time<<std::endl;		
	
			auto this_shower_dir_plane = shower_dir_plane.at(plane);
			double fradial_dist_sshit_vtx_notrack = distToVtx(fTimeToCMConstant, fWireToCMConstant, fRadius, time, wire, this_sshit);
			if (fradial_dist_sshit_vtx_notrack> fRadius){
				std::cout<<"ERROR: in fillROITree ssnet hit in ROI minus shower outside radius, distance = "<<fradial_dist_sshit_vtx_notrack<<std::endl;	
			}
       			double fopening_angle_shower_sshit_notrack = angleFromShower(this_shower_dir_plane, this_sshit, time,wire);
			dist_ang.push_back(TVector2(fradial_dist_sshit_vtx_notrack,  fopening_angle_shower_sshit_notrack));
			//(&fROITree)->Fill();
		}//for each plane
	}//for each hit
	return dist_ang;
}//fill ROI Tree
*/

//returns the index of the value in the vector if contained, -1 otherwise
int contains_at_ind(int n, std::vector<int> list){
	int my_ind = -1;
	int this_ind = 0;
	for (auto const& item : list){
		if (item ==n){
			//std::cout<<"match, the item in the list is "<<item<<" and the index in the _ROIHitLits is "<<n<<std::endl;
			my_ind = this_ind;
			break;
		}
		this_ind++;
	}
	return my_ind;
}

//returns true if the int is an item in the vector
bool contains(int n, std::vector<int> _listnotrack){
	bool contains =false;
	for(auto const& item : _listnotrack){
		if(item == n){
			contains = true;
			//break;
		}
	}	
	return contains;
}

//double fitFunction(double *x, double *par){
//      return TF1::gaus(x,par) + TF1::gaus(x,&par[2]);
//}

double mygaus(double x, double A, double mean, double sigma){
	return A*exp(-(x-mean)*(x-mean)/(2.0*sigma*sigma));
}

double fitf(Double_t *x,Double_t *par) {
        return  mygaus(x[0], par[0],par[1],par[2])+ mygaus(x[0], par[3],par[4],par[5]);
//      return  mygaus(x[0], par[0],par[1],par[2]);
  
}




//gets the two other vertices of a shower for a given plane by computing a circle in 3D
//and then projecting each point to a given plane
//The two points furtherest away from the shower start in either direction are taken to be the vertices
std::vector<TVector2> getMinMaxShowerPlane(std::vector<std::array<double, 3>> coneRim,int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo, detinfo::DetectorProperties const& detprop, TVector2 shower_start_plane, std::vector<TVector2> min_max){
//	TVector2 first =  TVector2(1, 1);
//	std::vector<TVector2> plane_points{first};
	std::vector<TVector2> plane_points;	
	//std::cout<<"flag.seg.0"<<std::endl; 
	for(std::array<double, 3> point : coneRim ){
	//	std::cout<<"flag.seg.1"<<std::endl;
		//calc 2d point of cone rim
		double x = point.at(0);
		double y = point.at(1);
		double z = point.at(2);
	//	std::cout<<"flag.seg.2"<<std::endl;
		TVector2 this_point_2d = TVector2(calcTime(x, plane, fTPC,fCryostat, detprop), calcWire(y, z, plane, fTPC, fCryostat, geo));
	//	std::cout<<"this 2d point on the circle = "<<this_point_2d.X()<<" "<<this_point_2d.Y()<<std::endl;
	//	std::cout<<"flag.seg.3"<<std::endl;
		plane_points.push_back(this_point_2d);
	//	int ind = plane_points.size()-1;
	//	if (ind>=0){
	//		std::cout<<"at index "<<ind<<" the address is "<<&plane_points[ind]<<std::endl;
	//	}		
	}
	//for (size_t vtx_index = 0; vtx_index != coneRim.size(); ++vtx_index){
	//	std::vector<TVector2>* plane_points_copy = plane_points;
	//	std::cout<<"this address = "<<std::to_string(plane_points_copy)<<std::endl;
	//	plane_points_copy++;
	//}
	//std::cout<<"this address = "<<&item<<std::endl;
	return plane_points;
	//min_max = &plane_points;
	//return *min_max;
	//return plane_points;
	//find min and max of rim
	//TVector2 min = plane_points.at(0);
	//TVector2 max = plane_points.at(0);
//	TVector2 max_x_vec = plane_points.at(0);
//	TVector2 max_y_vec = plane_points.at(0);

	//double dist_min = distToVtx2d(min, shower_start_plane);
	//double dist_max = distToVtx2d(max, shower_start_plane);
//	double start_x_plane = shower_start_plane.X();
//	double start_y_plane = shower_start_plane.Y();


//	double dist_x_xmax =  abs(start_x_plane - max_x_vec.X());
	//double dist_y_xmax = ;
	//double dist_x_ymax = ;
//	double dist_y_ymax =  abs(start_y_plane - max_y_vec.Y());

//	for(TVector2 vec:plane_points){
//		double this_x = vec.X();
//		double this_x_dist = abs(start_x_plane - this_x);

//		double this_y = vec.Y();
//		double this_y_dist = abs(start_y_plane - this_y);

		//double this_dist = distToVtx2d(vec, shower_start_plane);
		//std::cout<<"the dist of this point on the circle to the shower start is "<<this_dist<<std::endl;
//		if (this_y_dist>dist_y_ymax){
//			dist_y_ymax=this_y_dist;
//			std::cout<<"updating y max dist to "<<dist_y_ymax<<std::endl;
//			max_y_vec = vec;
//		}	
//		if (this_x_dist>dist_x_xmax){
//			dist_x_xmax=this_x_dist;
//			std::cout<<"updating x max dist to "<<dist_x_xmax<<std::endl;
//			max_x_vec= vec;
//		}	
//	}
	
//	std::vector<TVector2> this_min_max;
//	this_min_max.push_back(max_x_vec);
//	this_min_max.push_back(max_y_vec);
//	min_max = &this_min_max;
//	return min_max;
}

//converts radius for ROI in units of time and wire from sm
//double getRadius(double rad_in_cm){
	//convert cm to time
	//convert cm to wire
	//r = sqrt(t^2 + w^2)
//	return 500;
//}//getRadius


} // local namespace
