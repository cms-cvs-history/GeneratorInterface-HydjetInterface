/*
 *
 * Generates HYDJET HepMC events
 *
 * Original Author: Camelia Mironov
*/

#include <iostream>
#include <cstdio>
#include "time.h"

#include "GeneratorInterface/HydjetInterface/interface/HydjetSource.h"
#include "GeneratorInterface/HydjetInterface/interface/HydjetWrapper.h"
#include "GeneratorInterface/CommonInterface/interface/PythiaCMS.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "CLHEP/HepMC/include/PythiaWrapper6_2.h"
#include "CLHEP/HepMC/ConvertHEPEVT.h"
#include "CLHEP/HepMC/CBhepevt.h"

using namespace edm;
using namespace std;

HepMC::ConvertHEPEVT conv;

HydjetSource :: HydjetSource(const ParameterSet & pset, 
			     InputSourceDescription const& desc):
GeneratedInputSource(pset, desc), evt(0), 
abeamtarget_(pset.getUntrackedParameter<double>("aBeamTarget",207.)),
bfixed_(pset.getUntrackedParameter<double>("bFixed",0.)),
bmax_(pset.getUntrackedParameter<double>("bMax",0.)),
bmin_(pset.getUntrackedParameter<double>("bMin",0.)),
cflag_(pset.getUntrackedParameter<int>("cFlag",0)),
comenergy(pset.getUntrackedParameter<double>("comEnergy",5500.)),
hyMode(HydjetSource::kHydroQJets),
maxEventsToPrint_(pset.getUntrackedParameter<int>("maxEventsToPrint",1)),
nhard_(0),
nmultiplicity_(pset.getUntrackedParameter<int>("nMultiplicity",30000)),
nsoft_(0),
pythiaPylistVerbosity_(pset.getUntrackedParameter<int>("pythiaPylistVerbosity",0))
{
  // Default constructor

  // PYLIST Verbosity Level
  // Valid PYLIST arguments are: 1, 2, 3, 5, 7, 11, 12, 13
  pythiaPylistVerbosity_ = pset.getUntrackedParameter<int>("pythiaPylistVerbosity",0);
  cout << "Pythia PYLIST verbosity level = " << pythiaPylistVerbosity_ << endl;

  //Max number of events printed on verbosity level 
  maxEventsToPrint_ = pset.getUntrackedParameter<int>("maxEventsToPrint",0);
  cout << "Number of events to be printed = " << maxEventsToPrint_ << endl;


  //initialize pythia
  hyjpythia_init();

  //initialize hydro part
  hyjhydro_init();

  // Read PYTHIA parameters
   ParameterSet pythia_params = 
    pset.getParameter<ParameterSet>("PythiaParameters") ;
   vector<string> pars = 
      pythia_params.getParameter<vector<string> >("pythia");
    
    // Loop over all parameters and stop in case of mistake
    for( vector<string>::const_iterator  
	   itPar = pars.begin(); itPar != pars.end(); ++itPar ) 
      {
	static string sRandomValueSetting("MRPY(1)");
	if( 0 == itPar->compare(0,sRandomValueSetting.size(),sRandomValueSetting) ) 
	  {
	    throw edm::Exception(edm::errors::Configuration,"PythiaError")
	      <<" attempted to set random number seed.";
	  }
	if( ! call_pygive(*itPar) ) 
	  {
	    throw edm::Exception(edm::errors::Configuration,"PythiaError") 
	      <<" pythia did not accept the following \""<<*itPar<<"\"";
	  }
      }


 // Read HYDJET parameters if they exist
    ParameterSet hydjet_params = 
     pset.getParameter<ParameterSet>("HydjetParameters") ;
    
  // Read the HYDJET parameters from the set
    vector<string> pars_hyj = 
      hydjet_params.getParameter<vector<string> >("hydjet");
    
  // Loop over all parameters and stop in case of mistake
    for( vector<string>::const_iterator  
	   itPar = pars_hyj.begin(); itPar != pars_hyj.end(); ++itPar ) 
      {     
	if( ! call_hyjgive(*itPar) ) 
	  {
	    throw edm::Exception(edm::errors::Configuration,"HYDJET Error") 
	      <<" HYDJET did not accept the following \""<<*itPar<<"\"";
	  }
      }
    
    
    if( hyMode != HydjetSource::kHydroOnly ) 
      { 
	call_pyinit("CMS", "p", "p", comenergy);
      }
    
    cout<<endl;
    produces<HepMCProduct>();

}


//_____________________________________________________________________
HydjetSource::~HydjetSource()
{
  // distructor
  call_pystat(1);

  clear();
}


//________________________________________________________________
HepMC::GenParticle* HydjetSource::build_particle(int index) 				
{
  // Builds particle object corresponding to index in lujets


  HepMC::GenParticle* p = 
    new HepMC::GenParticle( HepLorentzVector(lujets.p[0][index],// px
					     lujets.p[1][index],// py
					     lujets.p[2][index],// pz
					     lujets.p[3][index] // E
					     ),
			    lujets.k[1][index],// id
			    lujets.k[0][index] // status
			    );
  p->suggest_barcode(index);

  return p;
}


//____________________________________________________________________
bool HydjetSource::build_vertices(int i, 
					   vector<HepMC::GenParticle*>& luj_entries,
					   HepMC::GenEvent* evt)
{
  // Build a production vertex for particle with index i in lujets
  // add the vertex to the event

  // fix: need to fix to look for the second mothers in case of flavor 
  // K(I,2)=91-94; cluster, string, indep, CMshower

  HepMC::GenParticle* pi       = luj_entries[i]; 
  HepMC::GenVertex* prod_vtx_i = pi->production_vertex();
  int mother_i                 = lujets.k[2][i];    
    
  if ( !prod_vtx_i && mother_i > 0 )
    {
      prod_vtx_i = luj_entries[mother_i]->end_vertex(); //decay vertex of the mother
      if (prod_vtx_i) 
	{
	  // if the decay vertex of its mother exists
	  // assign it to the particle, as the production vertex
	  prod_vtx_i->add_particle_out( pi );
	} 
    }
	 
  HepLorentzVector prod_pos( lujets.v[0][mother_i],
			     lujets.v[1][mother_i],
			     lujets.v[2][mother_i],
			     lujets.v[4][mother_i]
			     ); 
  if ( !prod_vtx_i && (mother_i > 0 || prod_pos != HepLorentzVector(0,0,0,0)) )
    {
      prod_vtx_i = new HepMC::GenVertex();
      prod_vtx_i->add_particle_out( pi );	
      evt->add_vertex( prod_vtx_i );
    }
	   
  if ( prod_vtx_i && prod_vtx_i->position()==HepLorentzVector(0,0,0,0) )
    {
      prod_vtx_i->set_position( prod_pos );
    }
	    
  //  check the consistency of the end_vertices
  if ( prod_vtx_i && mother_i > 0 )
    {
      if ( !luj_entries[mother_i]->end_vertex() )
	{
	  // if end vtx of the  mother isn't specified, do it now
	  prod_vtx_i->add_particle_in( luj_entries[mother_i] );
	} else if( luj_entries[mother_i]->end_vertex() != prod_vtx_i )
	  {
	    //error. the decay vtx of the mother is different from the daughter production vtx
	    cerr << "HydjetSource::build_production_vertex: "<<
	      "inconsistent mother/daughter produced vtx in event!" << endl;
	    luj_entries[mother_i]->end_vertex()->print(cerr);
	    prod_vtx_i->print(cerr);

	  }
    }
  return true;
}


//_____________________________________________________________________
bool HydjetSource::call_hyjgive(const std::string& iParm ) 
{
  // Set Hydjet parameters

  bool accepted = true;

  if( !strncmp(iParm.c_str(),"nhsel",5))
      hyjpar.nhsel = atoi(&iParm[strcspn(iParm.c_str(),"=")+1]);
  else if( !strncmp(iParm.c_str(),"ptmin",5))
      hyjpar.ptmin = atof(&iParm[strcspn(iParm.c_str(),"=")+1]);
  else if( !strncmp(iParm.c_str(),"fpart",5))
      hyflow.fpart = atof(&iParm[strcspn(iParm.c_str(),"=")+1]);
  else if( !strncmp(iParm.c_str(),"ylfl",4))
      hyflow.ylfl  = atof(&iParm[strcspn(iParm.c_str(),"=")+1]);
  else if( !strncmp(iParm.c_str(),"ytfl",4))
      hyflow.ytfl  = atof(&iParm[strcspn(iParm.c_str(),"=")+1]);
  else accepted = false;

  return accepted;
}


//______________________________________________________________________
bool HydjetSource::call_pygive(const std::string& iParm ) 
{
  // Set Pythia parameters

  int numWarn = pydat1.mstu[26]; //# warnings
  int numErr = pydat1.mstu[22];// # errors
  // call the fortran routine pygive with a fortran string
  PYGIVE( iParm.c_str(), iParm.length() );  
  // if an error or warning happens it is problem
  return pydat1.mstu[26] == numWarn && pydat1.mstu[22] == numErr;   
}


//____________________________________________________________________
void HydjetSource::clear()
{

}


//_____________________________________________________________________
bool HydjetSource::get_hydjet_particles(HepMC::GenEvent* evt)
{
  // Hard particles. The first nhard_ lines form lujets array.
  // It corresponds to hard multijet part of the event: hard 
  // pythia/pyquen sub-events (sub-collisions) for a given event
  // PYTHIA/PYQUEN-induced, initial protons and partons, final partons, 
  // strings, unstable and stable hadrons - full multijet story a la pythia-pyjets

  // Soft particles. The last nsoft_ lines of lujets
  // It corresponds to HYDRO-induced, hadrons only

  // return T/F if succes/failure

  int lujetsEntries= nhard_+nsoft_;

  // create a particle instance for each lujets entry and fill a map
  // create a vector which maps from the lujets particle index to the 
  // GenParticle address

  vector<HepMC::GenParticle*> luj_entries(lujetsEntries);
  for ( int i1 = 0; i1<lujetsEntries; i1++ )
    {     
      luj_entries[i1] = build_particle(i1);
    }

  // loop over particles again to create vertices
  for ( int i2 = 0; i2<lujetsEntries; i2++ )
    {
      build_vertices( i2,luj_entries,evt );
    } 

  // handle the case with particles comming from nowhere 
  // no mothers, no daughters
  for ( int i3 = 0; i3<lujetsEntries; i3++ )
    {
      if ( !luj_entries[i3]->end_vertex() &&
	   !luj_entries[i3]->production_vertex() )
	{
	  HepMC::GenVertex* prod_vtx_i3 = new  HepMC::GenVertex();
	  prod_vtx_i3->add_particle_out( luj_entries[i3] ) ;
	  evt->add_vertex( prod_vtx_i3 );
	} 
    }
  return true;
}


//______________________________________________________________
bool HydjetSource::hyjhydro_init()
{
  //initialize hydjet HYDRO part

  //hydjet mode	
  call_hyjgive("nhsel=2");

  //minimum pT hard
  call_hyjgive("ptmin=-1.");

  // fraction of soft (hydro-induced) hadronic multiplicity, 
  // proportional to the number of nucleons-participants
  call_hyjgive("fpart=0.");

  //max longitudinal flow rapidity
  call_hyjgive("ylfl=5.");

  //max transverse flow rapidity
  call_hyjgive("ytfl=1.");

  return true;
}


//____________________________________________________________________
bool HydjetSource::hyjpythia_init()
{
  //initialize PYTHIA

  //random number seed
  edm::Service<RandomNumberGenerator> rng;
  uint32_t seed = rng->mySeed();
  ostringstream sRandomSet;
  sRandomSet <<"MRPY(1)="<<seed;
  call_pygive(sRandomSet.str());

  // QCD dijet production
  call_pygive("MSEL=1");

  // to avoid stopping run
  call_pygive("MSTU(21)=1");

  // tolerance parameter to adjust fragmentation
  call_pygive("PARU(14)=1.");

  // pp multiple scattering off
  call_pygive("MSTP(81)=0");

  return true;
}


//_____________________________________________________________________
bool HydjetSource::produce(Event & e)
{
  // generate single event

  cout<<".";
  
  nsoft_    = 0;
  nhard_    = 0;

  HYDRO(abeamtarget_,cflag_,bmin_,bmax_,bfixed_,nmultiplicity_);
  nsoft_    = hyfpar.nhyd;
  nhard_    = hyfpar.npyt;

  // event information
  HepMC::GenEvent* evt      = new HepMC::GenEvent();
  get_hydjet_particles(evt); 

  evt->set_signal_process_id(pypars.msti[0]);      // type of the process
  evt->set_event_scale(pypars.pari[16]);           // Q^2
  evt->set_event_number(numberEventsInRun() - remainingEvents() - 1);

  if(evt) 
    {
      auto_ptr<HepMCProduct> bare_product(new HepMCProduct());
      bare_product->addHepMCData(evt );
      e.put(bare_product);

      // print PYLIST info
      if(event() <= maxEventsToPrint_ && pythiaPylistVerbosity_) 	
	call_pylist(pythiaPylistVerbosity_);      
    }
  
  return true;
}


//________________________________________________________________

