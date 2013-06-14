#include <stdlib.h> // exit
// This Class' Header ------------------
#include "WireHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "GeaneTrackRep2.h"
#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"

// Class Member definitions -----------

ClassImp(WireHit)


WireHit::~WireHit()
{}

WireHit::WireHit()
  : WireRecoHit(NparHitRep)
{}

WireHit::WireHit(TVector3 point1,TVector3 point2, double rdrift, double resR, bool smearing)
  : WireRecoHit(NparHitRep){

     double x1 = point1.X();
     double y1 = point1.Y();
     double z1 = point1.Z();
     double x2 = point2.X();
     double y2 = point2.Y();
     double z2 = point2.Z();

     fHitCov[0][0] = 0.0;
     fHitCov[1][1] = 0.0;
     fHitCov[2][2] = 0.0;
     fHitCov[3][3] = 0.0;
     fHitCov[4][4] = 0.0;
     fHitCov[5][5] = 0.0;
     fHitCov[6][6] = resR*resR;


     fHitCoord[0][0] = x1;
     fHitCoord[1][0] = y1;
     fHitCoord[2][0] = z1;
     fHitCoord[3][0] = x2;
     fHitCoord[4][0] = y2;
     fHitCoord[5][0] = z2;
     // used for the check with vacuum
     if (smearing) {
        fHitCoord[6][0] = gRandom->Gaus(rdrift,resR); // it might be negative, but leave this.
     } else {
        fHitCoord[6][0] = rdrift;
     }

     //printf("fHitCoord: %lf %lf %lf\n",fHitCoord[0][0],fHitCoord[1][0],fHitCoord[2][0]);
  }

GFAbsRecoHit* 
WireHit::clone(){
   return new WireHit(*this);
}


   TMatrixT<double>
WireHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
   if ((dynamic_cast<const RKTrackRep*>(stateVector) != NULL)) {
      fprintf(stderr,"TKTrackRep is not considered\n");
      exit(1);
   }
   else if ((dynamic_cast<const GeaneTrackRep2*>(stateVector) != NULL)) {
      //I know, since this is the same everytime, it could be done in the
      //the constructor, but I do it here anyway, to make clear that in the
      //case of several track-reps per hit, it would have to be done here
      //    fHMatrix.ResizeTo(NparHitRep,5);
      TMatrixT<double> HMatrix(1,6);
      //TMatrixT<double> HMatrix(NparHitRep,6); // 2013/3/12

      HMatrix[0][0] = 0.;  // q/p
      HMatrix[0][1] = 0.; 
      HMatrix[0][2] = 0.; 
      HMatrix[0][3] = 1.;
      HMatrix[0][4] = 0.;
      HMatrix[0][5] = 0.;

      return HMatrix;
   }
   else {
      std::cerr << "WireHit can only handle state"
         << " vectors of type GeaneTrackRep2 -> abort" << std::endl;
      throw;
   }

}
double WireHit::get_smeared_rdrift()
{
   return fHitCoord[6][0];
}


