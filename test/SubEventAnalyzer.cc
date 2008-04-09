// -*- C++ -*-
//
// Package:    SubEventAnalyzer
// Class:      SubEventAnalyzer
// 
/**\class SubEventAnalyzer SubEventAnalyzer.cc yetkin/SubEventAnalyzer/src/SubEventAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Yetkin Yilmaz
//         Created:  Tue Dec 18 09:44:41 EST 2007
// $Id: SubEventAnalyzer.cc,v 1.1 2008/01/28 16:08:46 yilmaz Exp $
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"


#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

#include "TH1D.h"
#include "TNtuple.h"

using namespace std;


//
// class decleration
//

class SubEventAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SubEventAnalyzer(const edm::ParameterSet&);
      ~SubEventAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

   TNtuple* subs;
   TNtuple* events;
   TNtuple* vertices;

   edm::Service<TFileService> fs;
   std::string modLabel_;
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
SubEventAnalyzer::SubEventAnalyzer(const edm::ParameterSet& iConfig)
  : modLabel_(iConfig.getUntrackedParameter<string>("moduleLabel","source"))
{
   //now do what ever initialization is needed

}


SubEventAnalyzer::~SubEventAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SubEventAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  /* Events, Sub Events and Vertices are looped over seperately for different 
     consistency checks. Particles are looped over for three times; this is to
     check the consistency of different ways to reach the same info.
  */

   using namespace edm;
   using namespace HepMC;
  
   TH1D dNdEtaEv("h2","dNdEta of Whole Event",30,-6,6);

   Handle<HepMCProduct> mc;
// iEvent.getByLabel("source",mc);
   iEvent.getByLabel(modLabel_,mc);

   double N = 0;
   double Nall = 0;
   double Etot = 0;
   double Con = 0;
   double Pthigh = 0;

   const GenEvent* evt = mc->GetEvent();

   GenEvent::particle_const_iterator beginpar = evt->particles_begin();
   GenEvent::particle_const_iterator endpar = evt->particles_end();
   for(GenEvent::particle_const_iterator p = beginpar; p != endpar; ++p){

      HepMC::GenParticle* par = *p;

      int stat = par->status();
      if(stat != 1) continue;
      double eta = par->momentum().eta();
      if(eta > 2) continue;
      double pt = par->momentum().perp();
      if(pt > Pthigh) Pthigh = pt;
   }

   GenEvent::vertex_const_iterator begin = evt->vertices_begin();
   GenEvent::vertex_const_iterator end = evt->vertices_end();
   for(GenEvent::vertex_const_iterator v = begin; v != end; ++v){
      HepMC::GenVertex* vert = *v;
   }

   Handle<GenHIEvent> genhi;
// iEvent.getByLabel("source",genhi);
   iEvent.getByLabel(modLabel_,genhi);
   int nsub = genhi->getNsubs();

   for(int isub = 0; isub < nsub; isub++){

      TH1D dNdEta("h1","dNdEta of the Sub Event",30,-6,6);
      double con = 0;
      int index = 0;

      SubEvent sub = genhi->getSubEvent(isub);
      GenVertex* v = sub.getVertex(*evt);
      GenVertex::vertex_iterator vbeg = v->vertices_begin(HepMC::relatives);
      GenVertex::vertex_iterator vend = v->vertices_end(HepMC::relatives);
      for(GenVertex::vertex_iterator vit = vbeg; vit != vend; ++vit){
	
      bool fill = 1;
      GenVertex* vert = *vit;
      HepMC::GenVertex::particle_iterator pf;
      HepMC::GenVertex::particle_iterator startf = vert->particles_begin( HepMC::family ); //                                                       
      HepMC::GenVertex::particle_iterator endf = vert->particles_end( HepMC::family );
      for ( pf = startf; pf != endf; ++pf ) {
	 GenParticle* par1 = *pf;
	 int idf = par1->pdg_id();
	 if(idf == 92) fill = 0; // not to include pythia "strings" in the momentum consistency check
      }

      double cn = vert->check_momentum_conservation();
      int pin = vert->particles_in_size();
      int pout = vert->particles_out_size();
      if(fill)vertices->Fill(cn,pin,pout);
      con = con + cn;

   }

   vector<GenParticle*> subparts = sub.getParticles(*evt);

   double pthigh = 0;
   double etot = 0;
   double npsall = subparts.size();
   double nps = 0;
   for(int ipar = 0; ipar < npsall; ++ipar){
      int stat = subparts[ipar]->status();
      if(stat > 10) continue;
      nps = nps + 1;
      double e = subparts[ipar]->momentum().e();
      etot = etot + e;
      double eta =  subparts[ipar]->momentum().eta();
      dNdEta.Fill(eta);
      dNdEtaEv.Fill(eta);
      if(eta > 2) continue;
      double pt = subparts[ipar]->momentum().perp();
      if(pt > pthigh) pthigh = pt;
   }

   N = N + nps;
   Nall = Nall + npsall;
   Etot = Etot + etot;

   double width = dNdEta.GetRMS();
   subs->Fill(pthigh,con,nps,npsall,etot,width);

   }

   double wid = dNdEtaEv.GetRMS();
   events->Fill(Pthigh,Con,N,Nall,Etot,wid,nsub);

}


// ------------ method called once each job just before starting event loop  ------------
void 
SubEventAnalyzer::beginJob(const edm::EventSetup&)
{

   vertices = fs->make<TNtuple>("vertices","","mom:npi:npo");
   subs = fs->make<TNtuple>("subs","","pthigh:mom:npr:np:etot:width");
   events = fs->make<TNtuple>("events","","pthigh:mom:npr:np:etot:width:nsub");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
SubEventAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SubEventAnalyzer);
