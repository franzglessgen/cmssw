#include <iostream>
#include <cmath>

#include "SimTracker/SiPhase2Digitizer/plugins/PixelBrickedDigitizerAlgorithm.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineSimService.h"

#include "CondFormats/SiPixelObjects/interface/GlobalPixel.h"
#include "CondFormats/DataRecord/interface/SiPixelQualityRcd.h"
#include "CondFormats/DataRecord/interface/SiPixelFedCablingMapRcd.h"
#include "CondFormats/DataRecord/interface/SiPixelLorentzAngleSimRcd.h"

// Geometry
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

using namespace edm;
using namespace sipixelobjects;

void PixelBrickedDigitizerAlgorithm::init(const edm::EventSetup& es) {
  if (use_ineff_from_db_)  // load gain calibration service fromdb...
    theSiPixelGainCalibrationService_->setESObjects(es);

  if (use_deadmodule_DB_)
    es.get<SiPixelQualityRcd>().get(SiPixelBadModule_);

  if (use_LorentzAngle_DB_)  // Get Lorentz angle from DB record
    es.get<SiPixelLorentzAngleSimRcd>().get(SiPixelLorentzAngle_);

  // gets the map and geometry from the DB (to kill ROCs)
  es.get<SiPixelFedCablingMapRcd>().get(fedCablingMap_);
  es.get<TrackerDigiGeometryRecord>().get(geom_);
}

PixelBrickedDigitizerAlgorithm::PixelBrickedDigitizerAlgorithm(const edm::ParameterSet& conf)
    : Phase2TrackerDigitizerAlgorithm(conf.getParameter<ParameterSet>("AlgorithmCommon"),
                                      conf.getParameter<ParameterSet>("PixelBrickedDigitizerAlgorithm")),
      odd_row_interchannelCoupling_next_row_(conf.getParameter<ParameterSet>("PixelBrickedDigitizerAlgorithm")
                                                 .getParameter<double>("Odd_row_interchannelCoupling_next_row")),
      even_row_interchannelCoupling_next_row_(conf.getParameter<ParameterSet>("PixelBrickedDigitizerAlgorithm")
                                                  .getParameter<double>("Even_row_interchannelCoupling_next_row")),
      odd_column_interchannelCoupling_next_column_(
          conf.getParameter<ParameterSet>("PixelBrickedDigitizerAlgorithm")
              .getParameter<double>("Odd_column_interchannelCoupling_next_column")),
      even_column_interchannelCoupling_next_column_(
          conf.getParameter<ParameterSet>("PixelBrickedDigitizerAlgorithm")
              .getParameter<double>("Even_column_interchannelCoupling_next_column")),
      apply_timewalk_(conf.getParameter<ParameterSet>("PixelBrickedDigitizerAlgorithm").getParameter<bool>("ApplyTimewalk")),
      timewalk_model_(
          conf.getParameter<ParameterSet>("PixelBrickedDigitizerAlgorithm").getParameter<edm::ParameterSet>("TimewalkModel")) {
  pixelFlag_ = true;
  LogDebug("PixelBrickedDigitizerAlgorithm") << "Algorithm constructed "
                                      << "Configuration parameters:"
                                      << "Threshold/Gain = "
                                      << "threshold in electron Endcap = " << theThresholdInE_Endcap_
                                      << "threshold in electron Barrel = " << theThresholdInE_Barrel_ << " "
                                      << theElectronPerADC_ << " " << theAdcFullScale_
                                      << " The delta cut-off is set to " << tMax_ << " pix-inefficiency "
                                      << addPixelInefficiency_;
}
PixelBrickedDigitizerAlgorithm::~PixelBrickedDigitizerAlgorithm() { LogDebug("PixelBrickedDigitizerAlgorithm") << "Algorithm deleted"; }
void PixelBrickedDigitizerAlgorithm::accumulateSimHits(std::vector<PSimHit>::const_iterator inputBegin,
                                                std::vector<PSimHit>::const_iterator inputEnd,
                                                const size_t inputBeginGlobalIndex,
                                                const uint32_t tofBin,
                                                const Phase2TrackerGeomDetUnit* pixdet,
                                                const GlobalVector& bfield) {
  // produce SignalPoint's for all SimHit's in detector
  // Loop over hits
  uint32_t detId = pixdet->geographicalId().rawId();
  size_t simHitGlobalIndex = inputBeginGlobalIndex;  // This needs to be stored to create the digi-sim link later

  // find the relevant hits
  std::vector<PSimHit> matchedSimHits;
  std::copy_if(inputBegin, inputEnd, std::back_inserter(matchedSimHits), [detId](auto const& hit) -> bool {
    return hit.detUnitId() == detId;
  });
  // loop over a much reduced set of SimHits
  for (auto const& hit : matchedSimHits) {
    LogDebug("PixelBrickedDigitizerAlgorithm") << hit.particleType() << " " << hit.pabs() << " " << hit.energyLoss() << " "
                                        << hit.tof() << " " << hit.trackId() << " " << hit.processType() << " "
                                        << hit.detUnitId() << hit.entryPoint() << " " << hit.exitPoint();

    std::vector<DigitizerUtility::EnergyDepositUnit> ionization_points;
    std::vector<DigitizerUtility::SignalPoint> collection_points;

    double signalScale = 1.0;
    // fill collection_points for this SimHit, indpendent of topology
    if (select_hit(hit, (pixdet->surface().toGlobal(hit.localPosition()).mag() * c_inv), signalScale)) {
      primary_ionization(hit, ionization_points);  // fills ionization_points

      // transforms ionization_points -> collection_points
      drift(hit, pixdet, bfield, ionization_points, collection_points);

      // compute induced signal on readout elements and add to _signal
      // hit needed only for SimHit<-->Digi link
      induce_signal(hit, simHitGlobalIndex, tofBin, pixdet, collection_points);
    }
    ++simHitGlobalIndex;
  }
}
//
// -- Select the Hit for Digitization
//
bool PixelBrickedDigitizerAlgorithm::select_hit(const PSimHit& hit, double tCorr, double& sigScale) {
  double time = hit.tof() - tCorr;
  return (time >= theTofLowerCut_ && time < theTofUpperCut_);
}

// ======================================================================
//
//  Add  Cross-talk contribution
//
// ======================================================================
void PixelBrickedDigitizerAlgorithm::add_cross_talk(const Phase2TrackerGeomDetUnit* pixdet) {
  if (!pixelFlag_)
    return;

  const Phase2TrackerTopology* topol = &pixdet->specificTopology();

  // cross-talk calculation valid for the case of 25x100 pixels
  const float pitch_first = 0.0025;
  const float pitch_second = 0.0100;

  // 0.5 um tolerance when comparing the pitch to accommodate the small changes in different TK geometrie (temporary fix)
  const double pitch_tolerance(0.0005);

  if (std::abs(topol->pitch().first - pitch_first) > pitch_tolerance ||
      std::abs(topol->pitch().second - pitch_second) > pitch_tolerance)
    return;

  uint32_t detID = pixdet->geographicalId().rawId();
  signal_map_type& theSignal = _signal[detID];
  signal_map_type signalNew;

  int numRows = topol->nrows();
  int numColumns = topol->ncolumns();

  for (auto& s : theSignal) {
    float signalInElectrons = s.second.ampl();  // signal in electrons

    auto hitChan = PixelDigi::channelToPixel(s.first);

    float signalInElectrons_odd_row_Xtalk_next_row = signalInElectrons * odd_row_interchannelCoupling_next_row_;
    float signalInElectrons_even_row_Xtalk_next_row = signalInElectrons * even_row_interchannelCoupling_next_row_;
    float signalInElectrons_odd_column_Xtalk_next_column =
        signalInElectrons * odd_column_interchannelCoupling_next_column_;
    float signalInElectrons_even_column_Xtalk_next_column =
        signalInElectrons * even_column_interchannelCoupling_next_column_;

    // subtract the charge which will be shared
    s.second.set(signalInElectrons - signalInElectrons_odd_row_Xtalk_next_row -
                 signalInElectrons_even_row_Xtalk_next_row - signalInElectrons_odd_column_Xtalk_next_column -
                 signalInElectrons_even_column_Xtalk_next_column);

    if (hitChan.first != 0) {
      auto XtalkPrev = std::make_pair(hitChan.first - 1, hitChan.second);
      int chanXtalkPrev = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second);
      if (hitChan.first % 2 == 1)
        signalNew.emplace(chanXtalkPrev,
                          DigitizerUtility::Amplitude(signalInElectrons_even_row_Xtalk_next_row, nullptr, -1.0));
      else
        signalNew.emplace(chanXtalkPrev,
                          DigitizerUtility::Amplitude(signalInElectrons_odd_row_Xtalk_next_row, nullptr, -1.0));
    }
    if (hitChan.first < numRows - 1) {
      auto XtalkNext = std::make_pair(hitChan.first + 1, hitChan.second);
      int chanXtalkNext = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkNext.first, XtalkNext.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkNext.first, XtalkNext.second);
      if (hitChan.first % 2 == 1)
        signalNew.emplace(chanXtalkNext,
                          DigitizerUtility::Amplitude(signalInElectrons_odd_row_Xtalk_next_row, nullptr, -1.0));
      else
        signalNew.emplace(chanXtalkNext,
                          DigitizerUtility::Amplitude(signalInElectrons_even_row_Xtalk_next_row, nullptr, -1.0));
    }

    if (hitChan.second != 0) {
      auto XtalkPrev = std::make_pair(hitChan.first, hitChan.second - 1);
      int chanXtalkPrev = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second);
      if (hitChan.second % 2 == 1)
        signalNew.emplace(chanXtalkPrev,
                          DigitizerUtility::Amplitude(signalInElectrons_even_column_Xtalk_next_column, nullptr, -1.0));
      else
        signalNew.emplace(chanXtalkPrev,
                          DigitizerUtility::Amplitude(signalInElectrons_odd_column_Xtalk_next_column, nullptr, -1.0));
    }
    if (hitChan.second < numColumns - 1) {
      auto XtalkNext = std::make_pair(hitChan.first, hitChan.second + 1);
      int chanXtalkNext = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkNext.first, XtalkNext.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkNext.first, XtalkNext.second);
      if (hitChan.second % 2 == 1)
        signalNew.emplace(chanXtalkNext,
                          DigitizerUtility::Amplitude(signalInElectrons_odd_column_Xtalk_next_column, nullptr, -1.0));
      else
        signalNew.emplace(chanXtalkNext,
                          DigitizerUtility::Amplitude(signalInElectrons_even_column_Xtalk_next_column, nullptr, -1.0));
    }
  }
  for (auto const& l : signalNew) {
    int chan = l.first;
    auto iter = theSignal.find(chan);
    if (iter != theSignal.end()) {
      iter->second += l.second.ampl();
    } else {
      theSignal.emplace(chan, DigitizerUtility::Amplitude(l.second.ampl(), nullptr, -1.0));
    }
  }
}

PixelBrickedDigitizerAlgorithm::TimewalkCurve::TimewalkCurve(const edm::ParameterSet& pset)
    : x_(pset.getParameter<std::vector<double>>("charge")), y_(pset.getParameter<std::vector<double>>("delay")) {
  if (x_.size() != y_.size())
    throw cms::Exception("Configuration")
        << "Timewalk model error: the number of charge values does not match the number of delay values!";
}

double PixelBrickedDigitizerAlgorithm::TimewalkCurve::operator()(double x) const {
  auto it = std::lower_bound(x_.begin(), x_.end(), x);
  if (it == x_.begin())
    return y_.front();
  if (it == x_.end())
    return y_.back();
  int index = std::distance(x_.begin(), it);
  double x_high = *it;
  double x_low = *(--it);
  double p = (x - x_low) / (x_high - x_low);
  return p * y_[index] + (1 - p) * y_[index - 1];
}

PixelBrickedDigitizerAlgorithm::TimewalkModel::TimewalkModel(const edm::ParameterSet& pset) {
  threshold_values = pset.getParameter<std::vector<double>>("ThresholdValues");
  const auto& curve_psetvec = pset.getParameter<std::vector<edm::ParameterSet>>("Curves");
  if (threshold_values.size() != curve_psetvec.size())
    throw cms::Exception("Configuration")
        << "Timewalk model error: the number of threshold values does not match the number of curves.";
  for (const auto& curve_pset : curve_psetvec)
    curves.emplace_back(curve_pset);
}

double PixelBrickedDigitizerAlgorithm::TimewalkModel::operator()(double q_in, double q_threshold) const {
  auto index = find_closest_index(threshold_values, q_threshold);
  return curves[index](q_in);
}

std::size_t PixelBrickedDigitizerAlgorithm::TimewalkModel::find_closest_index(const std::vector<double>& vec,
                                                                       double value) const {
  auto it = std::lower_bound(vec.begin(), vec.end(), value);
  if (it == vec.begin())
    return 0;
  else if (it == vec.end())
    return vec.size() - 1;
  else {
    auto it_upper = it;
    auto it_lower = --it;
    auto closest = (value - *it_lower > *it_upper - value) ? it_upper : it_lower;
    return std::distance(vec.begin(), closest);
  }
}
//
// -- Compare Signal with Threshold
//
bool PixelBrickedDigitizerAlgorithm::isAboveThreshold(const DigitizerUtility::SimHitInfo* hitInfo, float charge, float thr) {
  if (charge < thr)
    return false;
  if (apply_timewalk_ && hitInfo) {
    float corrected_time = hitInfo->time();
    double time = corrected_time + timewalk_model_(charge, thr);
    return (time >= theTofLowerCut_ && time < theTofUpperCut_);
  } else
    return true;
}



void PixelBrickedDigitizerAlgorithm::induce_signal(
    const PSimHit& hit,
    const size_t hitIndex,
    const uint32_t tofBin,
    const Phase2TrackerGeomDetUnit* pixdet,
    const std::vector<DigitizerUtility::SignalPoint>& collection_points) {
  // X  - Rows, Left-Right, 160, (1.6cm)   for barrel
  // Y  - Columns, Down-Up, 416, (6.4cm)
  const Phase2TrackerTopology* topol = &pixdet->specificTopology();
  uint32_t detID = pixdet->geographicalId().rawId();
  signal_map_type& theSignal = _signal[detID];
  const GlobalPoint& pos_glob = pixdet->position();
  

  // local map to store pixels hit by 1 Hit.
  using hit_map_type = std::map<int, float, std::less<int> >;
  hit_map_type hit_signal;

  // Assign signals to readout channels and store sorted by channel number
  // Iterate over collection points on the collection plane
  for (auto const& v : collection_points) {
    float CloudCenterX = v.position().x();  // Charge position in x
    float CloudCenterY = v.position().y();  //                 in y
    float SigmaX = v.sigma_x();             // Charge spread in x
    float SigmaY = v.sigma_y();             //               in y
    float Charge = v.amplitude();           // Charge amplitude


    // Find the maximum cloud spread in 2D plane , assume 3*sigma
    float CloudRight = CloudCenterX + clusterWidth_ * SigmaX;
    float CloudLeft = CloudCenterX - clusterWidth_ * SigmaX;
    float CloudUp = CloudCenterY + clusterWidth_ * SigmaY;
    float CloudDown = CloudCenterY - clusterWidth_ * SigmaY;

    // Define 2D cloud limit points
    LocalPoint PointRightUp = LocalPoint(CloudRight, CloudUp);
    LocalPoint PointLeftDown = LocalPoint(CloudLeft, CloudDown);

    // This points can be located outside the sensor area.
    // The conversion to measurement point does not check for that
    // so the returned pixel index might be wrong (outside range).
    // We rely on the limits check below to fix this.
    // But remember whatever we do here THE CHARGE OUTSIDE THE ACTIVE
    // PIXEL ARE IS LOST, it should not be collected.

    // Convert the 2D points to pixel indices
    MeasurementPoint mp = topol->measurementPosition(PointRightUp);
    //MeasurementPoint mp_bricked = topol->measurementPosition(PointRightUp);
    int IPixRightUpX = static_cast<int>(std::floor(mp.x()));  // cast reqd.
    //int IPixRightUpY = static_cast<int>(std::floor(mp.y()));

    int numColumns = topol->ncolumns();  // det module number of cols&rows
    int numRows = topol->nrows();
    IPixRightUpX = numRows > IPixRightUpX ? IPixRightUpX : numRows - 1;

    //Specific to bricked geometry 
    int IPixRightUpY = static_cast<int>(mp.y() - 0.5*(IPixRightUpX%2) );

    mp = topol->measurementPosition(PointLeftDown);
    
    int IPixLeftDownX = static_cast<int>(std::floor(mp.x()));

    IPixLeftDownX = 0 < IPixLeftDownX ? IPixLeftDownX : 0;

    //Specific to bricked geometry 
    int IPixLeftDownY = static_cast<int>(mp.y() - 0.5*(IPixLeftDownX%2));//changed in case negative value	

    IPixRightUpY =  numColumns > IPixRightUpY    ? IPixRightUpY : numColumns - 1;
    IPixLeftDownY = 0 < IPixLeftDownY ? IPixLeftDownY : 0;


    // First integrate charge strips in x
    hit_map_type x;
    for (int ix = IPixLeftDownX; ix <= IPixRightUpX; ++ix) {  // loop over x index
      float xLB, LowerBound;
      // Why is set to 0 if ix=0, does it meen that we accept charge
      // outside the sensor?
      if (ix == 0 || SigmaX == 0.) {  // skip for surface segemnts
        LowerBound = 0.;
      } else {
        mp = MeasurementPoint(ix, 0.0);
        xLB = topol->localPosition(mp).x();
        LowerBound = 1 - calcQ((xLB - CloudCenterX) / SigmaX);
      }

      float xUB, UpperBound;
      if (ix == numRows - 1 || SigmaX == 0.) {
        UpperBound = 1.;
      } else {
        mp = MeasurementPoint(ix + 1, 0.0);
        xUB = topol->localPosition(mp).x();
        UpperBound = 1. - calcQ((xUB - CloudCenterX) / SigmaX);
      }
      float TotalIntegrationRange = UpperBound - LowerBound;  // get strip
      x.emplace(ix, TotalIntegrationRange);                   // save strip integral
    }

    // Now integrate strips in y. Two maps will be filled: y and y_bricked which will both be used for the induced signal.
    

	int IPixLeftDownY_bricked = IPixLeftDownY;
	int IPixRightUpY_bricked = IPixRightUpY;

     //Specific to bricked geometry
     IPixRightUpY = std::min( IPixRightUpY + int((IPixRightUpX%2)), numColumns-1);

     //This map will be twice as large as the non-bricked hit map in y to harbor both the integrated charge from the bricked and non-bricked columns.
     hit_map_type y;
    for (int iy = IPixLeftDownY; iy <= IPixRightUpY; ++iy) {  // loop over y index
      float yLB, LowerBound;
      if (iy == 0 || SigmaY == 0.) {
        LowerBound = 0.;
      } else {
        mp = MeasurementPoint(0.0, iy);
        yLB = topol->localPosition(mp).y();
        LowerBound = 1. - calcQ((yLB - CloudCenterY) / SigmaY);
      

	}

      float yUB, UpperBound;
      if (iy >= numColumns - 1 || SigmaY == 0.) {
        UpperBound = 1.;
      } else {
        

	mp = MeasurementPoint(0.0, iy + 1);
        yUB = topol->localPosition(mp).y();
        UpperBound = 1. - calcQ((yUB - CloudCenterY) / SigmaY);
      


	}

      float TotalIntegrationRange = UpperBound - LowerBound;
	
      //Even indices correspond to the non-bricked columns
      y.emplace(2*iy, TotalIntegrationRange);  // save strip integral

    }




   IPixLeftDownY_bricked = std::max( IPixLeftDownY_bricked - int((!(IPixLeftDownX%2))), 0);

    for (int iy = IPixLeftDownY_bricked; iy <= IPixRightUpY_bricked; ++iy) {  // loop over y index
      float yLB, LowerBound;
      if (iy == 0 || SigmaY == 0.) {
        LowerBound = 0.;
      } else {

	mp = MeasurementPoint(0.0, iy + 0.5);
        yLB = topol->localPosition(mp).y();
        LowerBound = 1. - calcQ((yLB - CloudCenterY) / SigmaY);


	}

      float yUB, UpperBound;
      if (iy >= numColumns  || SigmaY == 0.) { // This was changed for bricked pixels
        UpperBound = 1.;
      } else {

	mp = MeasurementPoint(0.0, iy + 1.5 );
        yUB = topol->localPosition(mp).y();
        UpperBound = 1. - calcQ((yUB - CloudCenterY) / SigmaY);

		}

      

      float TotalIntegrationRange = UpperBound - LowerBound;
      //Odd indices correspond to bricked columns
      y.emplace(2*iy + 1, TotalIntegrationRange);  // save strip integral
    

	}//loop over y index

    
    // Get the 2D charge integrals by folding x and y strips
    for (int ix = IPixLeftDownX; ix <= IPixRightUpX; ++ix) {    // loop over x index

        for (int iy = std::max(0,IPixLeftDownY - int((ix%2)))   ; iy <= IPixRightUpY ; ++iy) {  // loop over y index

	int iy_considered = iy*2 + ix%2;
	float ChargeFraction = Charge * x[ix] * y[iy_considered];

	int chanFired = -1;
        if (ChargeFraction > 0.) {
          chanFired =
              pixelFlag_ ? PixelDigi::pixelToChannel(ix, iy) : Phase2TrackerDigi::pixelToChannel(ix, iy);  // Get index
          // Load the amplitude
          hit_signal[chanFired] += ChargeFraction;
        }
			}  
			      }//x loop
			}//collection loop



   



  // Fill the global map with all hit pixels from this event
  float corr_time = hit.tof() - pixdet->surface().toGlobal(hit.localPosition()).mag() * c_inv;
  for (auto const& hit_s : hit_signal) {
    int chan = hit_s.first;
    theSignal[chan] +=
        (makeDigiSimLinks_ ? DigitizerUtility::Amplitude(hit_s.second, &hit, hit_s.second, corr_time, hitIndex, tofBin)
                           : DigitizerUtility::Amplitude(hit_s.second, nullptr, hit_s.second));
  }
}




















