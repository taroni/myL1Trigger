#ifndef STDPRESELECTIONSELECTOR
#define STDPRESELECTIONSELECTOR

#include <memory>
#include <algorithm>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

struct StdPreselectionSelector{
  public:
  explicit StdPreselectionSelector(const edm::ParameterSet& iConfig) {}

  virtual ~StdPreselectionSelector() {}

  // Collections to be selected
  typedef reco::GsfElectronCollection collection;
  typedef std::vector<reco::GsfElectronRef> container ; 
  //typedef std::vector<reco::GsfElectron> container ; 
  typedef container::const_iterator const_iterator;

  //define iterators with above typedef
  const_iterator begin () const { return selected_.begin () ; }
  const_iterator end () const { return  selected_.end () ; }

  void select(edm::Handle<reco::GsfElectronCollection> electrons, 
	      const edm::Event& iEvent , 
	      const edm::EventSetup& iEs)
  {
    selected_.clear();
    // Loop over electrons
    unsigned int i = 0 ;
    for ( reco::GsfElectronCollection::const_iterator eleIt = electrons->begin () ;
	  eleIt != electrons->end () ;
	  ++eleIt )
      {
	edm::Ref<reco::GsfElectronCollection> electronRef(electrons,i);
	if (selection(electronRef))
	  selected_.push_back (electronRef) ;
	//selected_.push_back (*eleIt) ;
	++i;
      }
  }
  
  private:
  container selected_ ;
   
  bool selection(edm::Ref<reco::GsfElectronCollection> eleRef);
  
};

#endif

