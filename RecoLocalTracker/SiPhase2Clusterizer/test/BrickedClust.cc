#include <map>
#include <vector>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include <TH2F.h>
#include <TH1F.h>
#include <THStack.h>
#include <TProfile.h>
#include <TProfile2D.h>
using Phase2TrackerGeomDetUnit = PixelGeomDetUnit;
struct ClusterHistos {
  // use TH1D instead of TH1F to avoid stauration at 2^31
  // above this increments with +1 don't work for float, need double

  TH1D* numberClusters[3];
  TH1D* clusterSize[3];
  TH1D* clusterCharge[3];
  TH1D* clusterSizeX[3];
  TH1D* clusterSizeY[3];
  TProfile2D* clusterSize2D[3];
  TProfile* clusterSizeXvsX[3];
  TProfile* clusterSizeYvsY[3];
  TProfile* clusterDXvsX[3];
  TProfile* clusterDYvsY[3];
  TH1I* clusterCol[3];
  TH1I* clusterRow[3];

  TH2F* globalPosXY[3];
  TH2F* localPosXY[3];

  TH1F* deltaX[3];
  TH1F* deltaY[3];
  TH1F* deltaX_P[3];
  TH1F* deltaY_P[3];

  TH1D* primarySimHits[3];
  TH1D* otherSimHits[3];
};

class BrickedClust : public edm::EDAnalyzer {
public:
  typedef std::map<unsigned int, std::vector<PSimHit> > SimHitsMap;
  typedef std::map<unsigned int, SimTrack> SimTracksMap;

  explicit BrickedClust(const edm::ParameterSet&);
  ~BrickedClust();
  void beginJob() override;
  void endJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  std::vector<edm::EDGetTokenT<edm::PSimHitContainer>> simHitTokensB_;
  std::vector<edm::EDGetTokenT<edm::PSimHitContainer>> simHitTokensE_;
  std::map<unsigned int, ClusterHistos>::iterator createLayerHistograms(unsigned int);
  std::vector<unsigned int> getSimTrackId(const edm::Handle<edm::DetSetVector<PixelDigiSimLink> >&,
                                          const DetId&,
                                          unsigned int);

  edm::ParameterSet config_;
  //edm::EDGetTokenT<Phase2TrackerCluster1DCollectionNew> tokenClusters_;
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> tokenClusters_;
  edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink> > tokenLinks_;
  edm::EDGetTokenT<edm::PSimHitContainer> tokenSimHitsB_;
  edm::EDGetTokenT<edm::PSimHitContainer> tokenSimHitsE_;
  edm::EDGetTokenT<edm::SimTrackContainer> tokenSimTracks_;

  bool catECasRings_;
  double simtrackminpt_;
  
  TH2F* trackerLayout_;
  TH2F* trackerLayoutXY_;
  TH2F* trackerLayoutXYBar_;
  TH2F* trackerLayoutXYEC_;

  std::map<unsigned int, ClusterHistos> histograms_;
  const std::pair<double, double> pixel_cell_transformation_(const std::pair<double, double> &pos,
                                                             unsigned int icell,
                                                             const std::pair<double, double> &pitch);
  const std::pair<double, double> pixel_cell_transformation_(const MeasurementPoint &pos,
                                                             unsigned int icell,
                                                             const std::pair<double, double> &pitch);
};

BrickedClust::BrickedClust(const edm::ParameterSet& conf)
    : config_(conf),
      //tokenClusters_(consumes<Phase2TrackerCluster1DCollectionNew>(conf.getParameter<edm::InputTag>("src"))),
      tokenLinks_(consumes<edm::DetSetVector<PixelDigiSimLink> >(conf.getParameter<edm::InputTag>("links"))),
      //tokenSimHitsB_(consumes<edm::PSimHitContainer>(conf.getParameter<edm::InputTag>("simhitsbarrel"))),
      //tokenSimHitsE_(consumes<edm::PSimHitContainer>(conf.getParameter<edm::InputTag>("simhitsendcap"))),
      tokenSimTracks_(consumes<edm::SimTrackContainer>(conf.getParameter<edm::InputTag>("simtracks"))),
      catECasRings_(conf.getParameter<bool>("ECasRings")),
      simtrackminpt_(conf.getParameter<double>("SimTrackMinPt")) {

      tokenClusters_ = consumes<edmNew::DetSetVector<SiPixelCluster>>(conf.getParameter<edm::InputTag>("src"));

  const std::vector<edm::InputTag> psimhit_B(config_.getParameter<std::vector<edm::InputTag>>("PSimHitSourceB"));

  for (const auto& itag : psimhit_B) {
    simHitTokensB_.push_back(consumes<edm::PSimHitContainer>(itag));
  }


  const std::vector<edm::InputTag> psimhit_E(config_.getParameter<std::vector<edm::InputTag>>("PSimHitSourceE"));

  for (const auto& itag : psimhit_E) {
    simHitTokensE_.push_back(consumes<edm::PSimHitContainer>(itag));
  }
							}


BrickedClust::~BrickedClust() {}

void BrickedClust::beginJob() {
  edm::Service<TFileService> fs;
  fs->file().cd("/");
  TFileDirectory td = fs->mkdir("Common");
  // Create common histograms
  trackerLayout_ = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 1200, 0.0, 120.0);
  trackerLayoutXY_ = td.make<TH2F>("XVsY", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
  trackerLayoutXYBar_ = td.make<TH2F>("XVsYBar", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
  trackerLayoutXYEC_ = td.make<TH2F>("XVsYEC", "x vs. y position", 2400, -120.0, 120.0, 2400, -120.0, 120.0);
}

void BrickedClust::endJob() {}

void BrickedClust::analyze(const edm::Event& event, const edm::EventSetup& eventSetup) {
  /*
     * Get the needed objects
     */

  // Get the clusters
  //edm::Handle<Phase2TrackerCluster1DCollectionNew> clusters;
  edm::Handle<edmNew::DetSetVector<SiPixelCluster>> clusters;
  //edm::Handle<SiPixelCluster> clusters;
  
  event.getByToken(tokenClusters_, clusters);

  
  //const edmNew::DetSetVector<SiPixelCluster> &input = *clusters;


  // Get the PixelDigiSimLinks
  edm::Handle<edm::DetSetVector<PixelDigiSimLink> > pixelSimLinks;
  event.getByToken(tokenLinks_, pixelSimLinks);

  // Get the SimHits
  
  /*edm::Handle<edm::PSimHitContainer> simHitsRaw[2];
  event.getByToken(tokenSimHitsB_, simHitsRaw[0]);
  event.getByToken(tokenSimHitsE_, simHitsRaw[1]);*/

  
  


  

  // Get the SimTracks
  edm::Handle<edm::SimTrackContainer> simTracksRaw;
  event.getByToken(tokenSimTracks_, simTracksRaw);

  // Get the geometry
  edm::ESHandle<TrackerGeometry> geomHandle;
  eventSetup.get<TrackerDigiGeometryRecord>().get("idealForDigi", geomHandle);
  const TrackerGeometry* tkGeom = &(*geomHandle);

  edm::ESHandle<TrackerTopology> tTopoHandle;
  eventSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  const TrackerTopology* tTopo = tTopoHandle.product();

  /*
     * Rearrange the simTracks
     */

  // Rearrange the simTracks for ease of use <simTrackID, simTrack>
  SimTracksMap simTracks;
  for (edm::SimTrackContainer::const_iterator simTrackIt(simTracksRaw->begin()); simTrackIt != simTracksRaw->end();
       ++simTrackIt) {
    if (simTrackIt->momentum().pt() > simtrackminpt_) {
      simTracks.emplace(simTrackIt->trackId(), *simTrackIt);
    }
  }

  /*
     * Validation   
     */

  // Number of clusters
  std::map<unsigned int, unsigned int> nClusters[3];
  std::map<unsigned int, unsigned int> nPrimarySimHits[3];
  std::map<unsigned int, unsigned int> nOtherSimHits[3];

  // Loop over modules
  //for (Phase2TrackerCluster1DCollectionNew::const_iterator DSViter = clusters->begin(); DSViter != clusters->end();
    //   ++DSViter) 
  for (edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter = clusters->begin(); DSViter != clusters->end();
       ++DSViter) {
    // Get the detector unit's id
    unsigned int rawid(DSViter->detId());
    DetId detId(rawid);
    unsigned int layer = (tTopo->side(detId) != 0) * 1000;  // don't split up endcap sides
    
    TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);

    
    layer += tTopo->layer(detId);
    
    //layer+= mType.isBarrel()*100;   
   // layer+= dtype().isBarrel()*100;   

 
   
    unsigned int det = 0;
    if (mType == TrackerGeometry::ModuleType::Ph2PXB  /*PSP*/) {
      det = 1;
    } else if (mType == TrackerGeometry::ModuleType::Ph2PXF  /*PSS || mType == TrackerGeometry::ModuleType::Ph2SS*/) {
      det = 2;
    } else {
      std::cout << "UNKNOWN DETECTOR TYPE!" << std::endl;
    }

     std::cout << layer << std::endl;
    // Get the geomdet
    const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
    const Phase2TrackerGeomDetUnit* tkDetUnit = dynamic_cast<const Phase2TrackerGeomDetUnit*>(geomDetUnit); 
    if (!geomDetUnit)
      continue;
   
    layer += 100*(geomDetUnit->type().isBarrel());


    // initialize the nhit counters if they don't exist for this layer
    auto nhitit(nClusters[det].find(layer));
    


    if (nhitit == nClusters[det].end()) {
      nClusters[det].emplace(layer, 0);
      nPrimarySimHits[det].emplace(layer, 0);
      nOtherSimHits[det].emplace(layer, 0);
    }



    // Create histograms for the layer if they do not yet exist
    std::map<unsigned int, ClusterHistos>::iterator histogramLayer(histograms_.find(layer));
    if (histogramLayer == histograms_.end())
      histogramLayer = createLayerHistograms(layer);


     std::cout << layer << " " << 3 << std::endl;

    // Loop over the clusters in the detector unit
    //for (edmNew::DetSet<Phase2TrackerCluster1D>::const_iterator clustIt = DSViter->begin(); clustIt != DSViter->end();
      //   ++clustIt) {
    for (edmNew::DetSet<SiPixelCluster>::const_iterator clustIt = DSViter->begin(); clustIt != DSViter->end();
         ++clustIt) {
      // determine the position
      //MeasurementPoint mpClu(clustIt->center(), clustIt->column() + 0.5);
      const std::vector<SiPixelCluster::Pixel> &pixelsVec = clustIt->pixels();
      histogramLayer->second.clusterSizeX[det]->Fill(clustIt->sizeX());
      histogramLayer->second.clusterSizeY[det]->Fill(clustIt->sizeY());
    


      const auto pitch = tkDetUnit->specificTopology().pitch();
 
      double cluster_tot = 0;
      std::pair<double, double> cluster_position({0.0, 0.0});

      for (unsigned int i = 0; i < pixelsVec.size(); ++i) {  // loop over pixels
        
	int pixx = pixelsVec[i].x;  // index as float=iteger, row index, 0-159
        int pixy = pixelsVec[i].y;  // same, col index, 0-415
        double adc = pixelsVec[i].adc;
	std::cout<<pixx<<" "<<pixy<<std::endl;
      	histogramLayer->second.clusterCharge[det]->Fill(adc);
      	histogramLayer->second.clusterCol[det]->Fill(pixx);
      	histogramLayer->second.clusterRow[det]->Fill(pixy);
         
	cluster_tot += adc;
          // Use the center of the pixel
        cluster_position.first += adc * (pixx + 0.5);
          
	  //std::cout<<(current_digi.row()%2)<<std::endl; 

	if (tkDetUnit->specificTopology().isBricked() && (int(pixx)%2))
	  cluster_position.second += adc * (pixy + 0.5+ 0.5);
		
	 //if (tkDetUnit->specificTopology().isBricked() && (current_digi.row()%2)) 
          
	else cluster_position.second += adc * (pixy + 0.5);
		
				}  // end if subdet (pixel loop)


      cluster_position.first /= double(cluster_tot);
      cluster_position.second /= double(cluster_tot);

      //histogramLayer->second.clusterSize2D[det]->Fill(cluster_position.first,cluster_position.second,clustIt->size());
 
      //MeasurementPoint mpClu(clustIt->x(), clustIt->y() + 0.5);
      MeasurementPoint mpClu(cluster_position.first, cluster_position.second);
      //Local3DPoint localPosClu = geomDetUnit->topology().localPosition(mpClu);
      //LocalPoint localPosClu = tkDetUnit->specificTopology().localPosition(mpClu);
      Local3DPoint localPosClu = tkDetUnit->specificTopology().localPosition(mpClu);
     
      unsigned int s = 2; 
      const std::pair<double, double> icell_psh = pixel_cell_transformation_(mpClu, s, pitch);
      histogramLayer->second.clusterSize2D[det]->Fill(icell_psh.first, icell_psh.second,clustIt->size());
      std::cout<<icell_psh.first<<" "<<icell_psh.second<<std::endl;

      histogramLayer->second.clusterSizeXvsX[det]->Fill(icell_psh.first, clustIt->sizeX());
      histogramLayer->second.clusterSizeYvsY[det]->Fill(icell_psh.second, clustIt->sizeY());
      

      //Global3DPoint globalPosClu = geomDetUnit->surface().toGlobal(localPosClu);
      Global3DPoint globalPosClu = geomDetUnit->surface().toGlobal(localPosClu);

      // Get all the simTracks that form the cluster
      std::vector<unsigned int> clusterSimTrackIds;
      for (unsigned int i = 0; i < pixelsVec.size(); ++i) {  // loop over pixels
      //for (unsigned int i(0); i < (unsigned int)clustIt->size(); ++i) 
        //unsigned int channel(Phase2TrackerDigi::pixelToChannel(clustIt->firstRow() + i, clustIt->column()));
        //unsigned int channel(PixelDigi::pixelToChannel(clustIt->firstRow() + i, clustIt->column()));
	        
	int chax = int(pixelsVec[i].x);  // index as float=iteger, row index, 0-159
        int chay = int(pixelsVec[i].y);  // same, col index, 0-415
                  

	unsigned int channel(PixelDigi::pixelToChannel(chax, chay));
      

        std::vector<unsigned int> simTrackIds(getSimTrackId(pixelSimLinks, detId, channel));
        for (auto it : simTrackIds) {
          bool add = true;
          for (unsigned int j = 0; j < clusterSimTrackIds.size(); ++j) {
            // only save simtrackids that are not present yet
            if (it == clusterSimTrackIds.at(j))
              add = false;
          }
          if (add)
            clusterSimTrackIds.push_back(it);
        }
      }
      std::sort(clusterSimTrackIds.begin(), clusterSimTrackIds.end());

      // find the closest simhit
      // this is needed because otherwise you get cases with simhits and clusters being swapped
      // when there are more than 1 cluster with common simtrackids
      const PSimHit* simhit = 0;  // bad naming to avoid changing code below. This is the closest simhit in x

     std::cout << simHitTokensB_.size() << std::endl;
  for (int i= 0, j = 0; i < (int)simHitTokensB_.size() && j < (int)simHitTokensE_.size(); i++, j++){

  //for (const auto& sh_token : simHitTokensB_) 
  
    edm::Handle<edm::PSimHitContainer> simHitsRaw[2];
    event.getByToken(simHitTokensB_[i], simHitsRaw[0]);
    event.getByToken(simHitTokensE_[j], simHitsRaw[1]);
     std::cout << simHitTokensB_.size() << std::endl;

      float minx = 10000;
      for (unsigned int simhitidx = 0; simhitidx < 2; ++simhitidx) {  // loop over both barrel and endcap hits
        for (auto simhitIt : *simHitsRaw[simhitidx]) {
          if (rawid == simhitIt.detUnitId()) {
            //std::cout << "=== " << rawid << " " << &simhitIt << " " << simhitIt.trackId() << " " << simhitIt.localPosition().x() << " " << simhitIt.localPosition().y() << std::endl;
            auto it = std::lower_bound(clusterSimTrackIds.begin(), clusterSimTrackIds.end(), simhitIt.trackId());
            if (it != clusterSimTrackIds.end() && *it == simhitIt.trackId()) {
              if (!simhit || fabs(simhitIt.localPosition().x() - localPosClu.x()) < minx) {
                minx = fabs(simhitIt.localPosition().x() - localPosClu.x());
                simhit = &simhitIt;
              }
            }
          }
        }
      }
    }
      if (!simhit){
        continue;}



      // only look at simhits from highpT tracks
      auto simTrackIt(simTracks.find(simhit->trackId()));
      if (simTrackIt == simTracks.end()){
        continue;}

      /*
             * Cluster related variables
             */

      // cluster size
     
      ++(nClusters[det].at(layer));
      ++(nOtherSimHits[det].at(layer));

     
       // cluster size
      
      histogramLayer->second.clusterSize[det]->Fill(clustIt->size());

      // Fill the position histograms
      trackerLayout_->Fill(globalPosClu.z(), globalPosClu.perp());
      trackerLayoutXY_->Fill(globalPosClu.x(), globalPosClu.y());
      if (fabs(layer) < 1000) {
        trackerLayoutXYBar_->Fill(globalPosClu.x(), globalPosClu.y());
      } else {
        trackerLayoutXYEC_->Fill(globalPosClu.x(), globalPosClu.y());
      }

      histogramLayer->second.localPosXY[det]->Fill(localPosClu.x(), localPosClu.y());
      if (fabs(layer) < 1000) {
        histogramLayer->second.globalPosXY[det]->Fill(globalPosClu.z(), globalPosClu.perp());
      } else {
        histogramLayer->second.globalPosXY[det]->Fill(globalPosClu.x(), globalPosClu.y());
      }

      // now get the position of the closest hit
      Local3DPoint localPosHit(simhit->localPosition());

      histogramLayer->second.deltaX[det]->Fill((localPosClu.x() - localPosHit.x()));
      histogramLayer->second.deltaY[det]->Fill((localPosClu.y() - localPosHit.y()));
      histogramLayer->second.clusterDXvsX[det]->Fill(icell_psh.first, (localPosClu.x() - localPosHit.x()));
      histogramLayer->second.clusterDYvsY[det]->Fill(icell_psh.second, (localPosClu.y() - localPosHit.y()));

      // Primary particles only
      unsigned int procT(simhit->processType());
      if (simTrackIt->second.vertIndex() == 0 and
          (procT == 2 || procT == 7 || procT == 9 || procT == 11 || procT == 13 || procT == 15)) {
        ++(nPrimarySimHits[det].at(layer));
        --(nOtherSimHits[det].at(layer));  // avoid double counting
        histogramLayer->second.deltaX_P[det]->Fill(localPosClu.x() - localPosHit.x());
        histogramLayer->second.deltaY_P[det]->Fill(localPosClu.y() - localPosHit.y());
      }
    }
  }


  // fill the counter histos per layer
  for (unsigned int det = 1; det < 3; ++det) {
    for (auto it : nClusters[det]) {
      auto histogramLayer(histograms_.find(it.first));
      if (histogramLayer == histograms_.end())
        std::cout << "*** SL *** No histogram for an existing counter! This should not happen!" << std::endl;
      histogramLayer->second.numberClusters[det]->Fill(it.second);
    }
    for (auto it : nPrimarySimHits[det]) {
      auto histogramLayer(histograms_.find(it.first));
      if (histogramLayer == histograms_.end())
        std::cout << "*** SL *** No histogram for an existing counter! This should not happen!" << std::endl;
      histogramLayer->second.primarySimHits[det]->Fill(it.second);
    }
    for (auto it : nOtherSimHits[det]) {
      auto histogramLayer(histograms_.find(it.first));
      if (histogramLayer == histograms_.end())
        std::cout << "*** SL *** No histogram for an existing counter! This should not happen!" << std::endl;
      histogramLayer->second.otherSimHits[det]->Fill(it.second);
    }
  }

     std::cout <<  " " << 9 << std::endl;

}

// Create the histograms
std::map<unsigned int, ClusterHistos>::iterator BrickedClust::createLayerHistograms(
    unsigned int ival) {
  std::ostringstream fname1, fname2;

  edm::Service<TFileService> fs;
  fs->file().cd("/");

  std::string tag;
  unsigned int id;
  if (ival < 1000) {
    id = ival%100;
    fname1 << "Barrel";
    fname2 << "Layer_" << id;
    tag = "";
  } else {
    int side = ival / 1000;
    id = ival - side * 1000;
    id = id%100;
    fname1 << "Forward";
    fname2 << "Layer_" << id;
      tag = "";
  }


  TFileDirectory td1 = fs->mkdir(fname1.str().c_str());
  TFileDirectory td = td1.mkdir(fname2.str().c_str());

  ClusterHistos local_histos;

  std::ostringstream histoName;

  /*
     * Number of clusters
     */

  //Homemade

  histoName.str("");
  histoName << "Cluster_Charge_Barrel" << tag.c_str() << id;
  local_histos.clusterCharge[1] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 100, 0 , 1000);
  
  histoName.str("");
  histoName << "Cluster_Charge_Forward" << tag.c_str() << id;
  local_histos.clusterCharge[2] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 100, 0 , 1000);


  histoName.str("");
  histoName << "Cluster_SizeX_Barrel" << tag.c_str() << id;
  local_histos.clusterSizeX[1] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 20, 0, 20);

  histoName.str("");
  histoName << "Cluster_SizeX_Forward" << tag.c_str() << id;
  local_histos.clusterSizeX[2] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 20, 0, 20);

  
  histoName.str("");
  histoName << "Cluster_SizeY_Barrel" << tag.c_str() << id;
  local_histos.clusterSizeY[1] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 20,0, 20);

  histoName.str("");
  histoName << "Cluster_SizeY_Forward" << tag.c_str() << id;
  local_histos.clusterSizeY[2] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 20, 0, 20);
  
  histoName.str("");
  histoName << "Cluster_SizeXvsX_Barrel" << tag.c_str() << id;
  local_histos.clusterSizeXvsX[1] = td.make<TProfile>(histoName.str().c_str(), histoName.str().c_str(), 50, 0 , 50);

  histoName.str("");
  histoName << "Cluster_SizeXvsX_Forward" << tag.c_str() << id;
  local_histos.clusterSizeXvsX[2] = td.make<TProfile>(histoName.str().c_str(), histoName.str().c_str(), 50, 0, 50);
 
  histoName.str("");
  histoName << "Cluster_SizeYvsY_Barrel" << tag.c_str() << id;
  local_histos.clusterSizeYvsY[1] = td.make<TProfile>(histoName.str().c_str(), histoName.str().c_str(), 200, 0 , 200);
  
  histoName.str("");
  histoName << "Cluster_SizeYvsY_Forward" << tag.c_str() << id;
  local_histos.clusterSizeYvsY[2] = td.make<TProfile>(histoName.str().c_str(), histoName.str().c_str(), 200, 0 , 200);



  histoName.str("");
  histoName << "Cluster_DXvsX_Barrel" << tag.c_str() << id;
  local_histos.clusterDXvsX[1] = td.make<TProfile>(histoName.str().c_str(), histoName.str().c_str(), 50, 0 , 50);
  histoName.str("");
  histoName << "Cluster_DXvsX_Forward" << tag.c_str() << id;
  local_histos.clusterDXvsX[2] = td.make<TProfile>(histoName.str().c_str(), histoName.str().c_str(), 50, 0 , 50);


  histoName.str("");
  histoName << "Cluster_DYvsY_Barrel" << tag.c_str() << id;
  local_histos.clusterDYvsY[1] = td.make<TProfile>(histoName.str().c_str(), histoName.str().c_str(), 200, 0 , 200);
  histoName.str("");
  histoName << "Cluster_DYvsY_Forward" << tag.c_str() << id;
  local_histos.clusterDYvsY[2] = td.make<TProfile>(histoName.str().c_str(), histoName.str().c_str(), 200, 0 , 200);
  

  histoName.str("");
  histoName << "Cluster_Size2D_Barrel" << tag.c_str() << id;
  local_histos.clusterSize2D[1] = td.make<TProfile2D>(histoName.str().c_str(), histoName.str().c_str(), 50, 0 , 50, 200, 0 , 200);
  histoName.str("");
  histoName << "Cluster_Size2D_Forward" << tag.c_str() << id;
  local_histos.clusterSize2D[2] = td.make<TProfile2D>(histoName.str().c_str(), histoName.str().c_str(), 50, 0 , 50, 200, 0 , 200);

  histoName.str("");
  histoName << "Cluster_Col_Barrel" << tag.c_str() << id;
  local_histos.clusterCol[1] = td.make<TH1I>(histoName.str().c_str(), histoName.str().c_str(), 1500,0, 1500);
  histoName.str("");
  histoName << "Cluster_Col_Forward" << tag.c_str() << id;
  local_histos.clusterCol[2] = td.make<TH1I>(histoName.str().c_str(), histoName.str().c_str(), 1500,0, 1500);
  histoName.str("");
  histoName << "Cluster_Row_Barrel" << tag.c_str() << id;
  local_histos.clusterRow[1] = td.make<TH1I>(histoName.str().c_str(), histoName.str().c_str(), 1500,0, 1500);
  histoName.str("");
  histoName << "Cluster_Row_Forward" << tag.c_str() << id;
  local_histos.clusterRow[2] = td.make<TH1I>(histoName.str().c_str(), histoName.str().c_str(), 1500,0, 1500);
//End homemade







  histoName.str("");
  histoName << "Number_Clusters_Barrel" << tag.c_str() << id;
  local_histos.numberClusters[1] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 50, 0., 0.);
  histoName.str("");
  histoName << "Number_Clusters_Forward" << tag.c_str() << id;
  local_histos.numberClusters[2] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 50, 0., 0.);

  /*
     * Cluster size
     */

  histoName.str("");
  histoName << "Cluster_Size_Barrel" << tag.c_str() << id;
  local_histos.clusterSize[1] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 21, -0.5, 20.5);
  histoName.str("");
  histoName << "Cluster_Size_Forward" << tag.c_str() << id;
  local_histos.clusterSize[2] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 21, -0.5, 20.5);

  /*
     * Local and Global positions
     */

  histoName.str("");
  histoName << "Local_Position_XY_Barrel" << tag.c_str() << id;
  local_histos.localPosXY[1] =
      td.make<TH2F>(histoName.str().c_str(), histoName.str().c_str(), 500, 0., 0., 500, 0., 0.);

  histoName.str("");
  histoName << "Local_Position_XY_Forward" << tag.c_str() << id;
  local_histos.localPosXY[2] =
      td.make<TH2F>(histoName.str().c_str(), histoName.str().c_str(), 500, 0., 0., 500, 0., 0.);

  histoName.str("");
  histoName << "Global_Position_XY_Barrel" << tag.c_str() << id;
  local_histos.globalPosXY[1] =
      td.make<TH2F>(histoName.str().c_str(), histoName.str().c_str(), 500, -120.0, 120.0, 2400, -120.0, 120.0);

  histoName.str("");
  histoName << "Global_Position_XY_Forward" << tag.c_str() << id;
  local_histos.globalPosXY[2] =
      td.make<TH2F>(histoName.str().c_str(), histoName.str().c_str(), 500, -120.0, 120.0, 2400, -120.0, 120.0);

  /*
     * Delta positions with SimHits
     */

  histoName.str("");
  histoName << "Delta_X_Barrel" << tag.c_str() << id;
  local_histos.deltaX[1] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 100, -0.02, 0.02);

  histoName.str("");
  histoName << "Delta_X_Forward" << tag.c_str() << id;
  local_histos.deltaX[2] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 100, -0.02, 0.02);

  

  histoName.str("");
  histoName << "Delta_Y_Barrel" << tag.c_str() << id;
  local_histos.deltaY[1] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 100, -0.2, 0.2);

  histoName.str("");
  histoName << "Delta_Y_Forward" << tag.c_str() << id;
  local_histos.deltaY[2] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 100, -3., 3.);

  /*
     * Delta position with simHits for primary tracks only
     */


  histoName.str("");
  histoName << "Delta_X_P" << tag.c_str() << id;
  local_histos.deltaX_P[1] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 100, -0.02, 0.02);

  histoName.str("");
  histoName << "Delta_X_P" << tag.c_str() << id;
  local_histos.deltaX_P[2] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 100, -0.02, 0.02);

  
  histoName.str("");
  histoName << "Delta_Y_P" << tag.c_str() << id;
  local_histos.deltaY_P[1] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 100, -0.2, 0.2);

  histoName.str("");
  histoName << "Delta_Y_P" << tag.c_str() << id;
  local_histos.deltaY_P[2] = td.make<TH1F>(histoName.str().c_str(), histoName.str().c_str(), 100, -3., 3.);

  /*
     * Information on the Digis per cluster
     */

  histoName.str("");
  histoName << "Primary_Digis_Barrel" << tag.c_str() << id;
  local_histos.primarySimHits[1] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 100, 0., 0.);
  histoName.str("");
  histoName << "Primary_Digis_Forward" << tag.c_str() << id;
  local_histos.primarySimHits[2] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 100, 0., 0.);

  histoName.str("");
  histoName << "Other_Digis_Barrel" << tag.c_str() << id;
  local_histos.otherSimHits[1] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 100, 0., 0.);
  histoName.str("");
  histoName << "Other_Digis_Forward" << tag.c_str() << id;
  local_histos.otherSimHits[2] = td.make<TH1D>(histoName.str().c_str(), histoName.str().c_str(), 100, 0., 0.);





  /*
     * End
     */

  std::pair<std::map<unsigned int, ClusterHistos>::iterator, bool> insertedIt(
      histograms_.insert(std::make_pair(ival, local_histos)));
  fs->file().cd("/");

  return insertedIt.first;
}

std::vector<unsigned int> BrickedClust::getSimTrackId(
    const edm::Handle<edm::DetSetVector<PixelDigiSimLink> >& pixelSimLinks, const DetId& detId, unsigned int channel) {
  std::vector<unsigned int> retvec;
  edm::DetSetVector<PixelDigiSimLink>::const_iterator DSViter(pixelSimLinks->find(detId));
  if (DSViter == pixelSimLinks->end())
    return retvec;
  for (edm::DetSet<PixelDigiSimLink>::const_iterator it = DSViter->data.begin(); it != DSViter->data.end(); ++it) {
    if (channel == it->channel()) {
      retvec.push_back(it->SimTrackId());
    }
  }
  return retvec;
}

const std::pair<double, double> BrickedClust::pixel_cell_transformation_(
    const MeasurementPoint& pos, unsigned int icell, const std::pair<double, double>& pitch) {
  return pixel_cell_transformation_(std::pair<double, double>({pos.x(), pos.y()}), icell, pitch);
}

const std::pair<double, double> BrickedClust::pixel_cell_transformation_(
    const std::pair<double, double>& pos, unsigned int icell, const std::pair<double, double>& pitch) {
  // Get the position modulo icell
  // Case the whole detector (icell==0)
  double xmod = pos.first;
  double ymod = pos.second;


  // Actual pixel cells
  if (icell != 0) {
    xmod = std::fmod(pos.first, icell);
    ymod = std::fmod(pos.second, icell);
  }

  const double xcell = xmod * pitch.first*1e4;
  const double ycell = ymod * pitch.second*1e4;

  return std::pair<double, double>({xcell, ycell});
}

DEFINE_FWK_MODULE(BrickedClust);
