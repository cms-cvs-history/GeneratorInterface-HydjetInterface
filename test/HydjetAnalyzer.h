#ifndef HydjetAnalyzer_H
#define HydjetAnalyzer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"

// forward declarations
class TFile;
class TH1D;


class HydjetAnalyzer : public edm::EDAnalyzer
{
 public:
  explicit HydjetAnalyzer(const edm::ParameterSet& );
  virtual ~HydjetAnalyzer() {} 

  virtual void analyze(const edm::Event&, const edm::EventSetup& );
  virtual void beginJob(const edm::EventSetup& );
  virtual void endJob();

 private:
 
  TH1D*        phdNdEta;           // histogram for dN/deta
  TH1D*        phdNdY;             // histogram for dN/dy
  TH1D*        phdNdPt;            // histogram for dN/dpt
  TH1D*        phdNdPhi;           // histogram for dN/dphi

  std::string  modLabel_;
  int          nevents_;

};

#endif
//module to analyze hydjet events
