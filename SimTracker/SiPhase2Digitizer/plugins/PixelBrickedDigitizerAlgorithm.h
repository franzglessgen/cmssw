#ifndef _SimTracker_SiPhase2Digitizer_PixelBrickedDigitizerAlgorithm_h
#define _SimTracker_SiPhase2Digitizer_PixelBrickedDigitizerAlgorithm_h

#include "SimTracker/SiPhase2Digitizer/plugins/Phase2TrackerDigitizerAlgorithm.h"

class PixelBrickedDigitizerAlgorithm : public Phase2TrackerDigitizerAlgorithm {
private:
  // A list of 2d points
  class TimewalkCurve {
  public:
    // pset must contain "charge" and "delay" of type vdouble
    TimewalkCurve(const edm::ParameterSet& pset);

    // linear interpolation
    double operator()(double x) const;

  private:
    std::vector<double> x_;
    std::vector<double> y_;
  };

  // Holds the timewalk model data
  class TimewalkModel {
  public:
    TimewalkModel(const edm::ParameterSet& pset);

    // returns the delay for given input charge and threshold
    double operator()(double q_in, double q_threshold) const;

  private:
    std::size_t find_closest_index(const std::vector<double>& vec, double value) const;

    std::vector<double> threshold_values;
    std::vector<TimewalkCurve> curves;
  };

public:
  PixelBrickedDigitizerAlgorithm(const edm::ParameterSet& conf);
  ~PixelBrickedDigitizerAlgorithm() override;

  // initialization that cannot be done in the constructor
  void init(const edm::EventSetup& es) override;

  // void initializeEvent();
  // run the algorithm to digitize a single det
  void accumulateSimHits(const std::vector<PSimHit>::const_iterator inputBegin,
                         const std::vector<PSimHit>::const_iterator inputEnd,
                         const size_t inputBeginGlobalIndex,
                         const uint32_t tofBin,
                         const Phase2TrackerGeomDetUnit* pixdet,
                         const GlobalVector& bfield) override;
  bool select_hit(const PSimHit& hit, double tCorr, double& sigScale) override;
  bool isAboveThreshold(const DigitizerUtility::SimHitInfo* hitInfo, float charge, float thr) override;
  void add_cross_talk(const Phase2TrackerGeomDetUnit* pixdet) override;

  // Addition four xtalk-related parameters to PixelBrickedDigitizerAlgorithm specific parameters initialized in Phase2TrackerDigitizerAlgorithm
  const double odd_row_interchannelCoupling_next_row_;
  const double even_row_interchannelCoupling_next_row_;
  const double odd_column_interchannelCoupling_next_column_;
  const double even_column_interchannelCoupling_next_column_;
  
  // Specific for bricked pixel
  void induce_signal(const PSimHit& hit,
                     const size_t hitIndex,
                     const unsigned int tofBin,
                     const Phase2TrackerGeomDetUnit* pixdet,
                     const std::vector<DigitizerUtility::SignalPoint>& collection_points);

  // Timewalk parameters
  bool apply_timewalk_;
  const TimewalkModel timewalk_model_;
};
#endif
