#ifndef WireHIT_HH
#define WireHIT_HH

// Base Class Headers ----------------
#include "GFRecoHitIfc.h"
#include "GFWireHitPolicy.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op

// Collaborating Class Declarations --

typedef GFRecoHitIfc<GFWireHitPolicy> WireRecoHit;

class WireHit : public WireRecoHit {
public:

  // Constructors/Destructors ---------
  WireHit();

  WireHit(TVector3 point1,TVector3 point2, double rdrift, double resR, bool smearing=true);
  virtual ~WireHit();

  virtual GFAbsRecoHit* clone();
  
  // Operations ----------------------
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);
public:
  double get_smeared_rdrift();


private:

  // Private Data Members ------------
  static const int NparHitRep = 7;
  

  // Private Methods -----------------

public:
  ClassDef(WireHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
