// This Class' Header ------------------
#include "PointHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "GeaneTrackRep2.h"
#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"

// Class Member definitions -----------

ClassImp(PointHit)


PointHit::~PointHit()
{}

PointHit::PointHit()
  : SpacepointRecoHit(NparHitRep)
{}

PointHit::PointHit(TVector3 point,double resXY, double resZ)
  : SpacepointRecoHit(NparHitRep){

     /*
  fHitCov[0][0] = res*res;
  fHitCov[1][1] = 1.2*1.2*res*res;
  fHitCov[2][2] = 4*res*res;
  GFDetPlane d;


  fHitCoord[0][0] = gRandom->Gaus(point.X(),res);
  fHitCoord[1][0] = gRandom->Gaus(point.Y(),1.2*res);
  fHitCoord[2][0] = gRandom->Gaus(point.Z(),2*res);
  */
  fHitCov[0][0] = resXY*resXY;
  fHitCov[1][1] = resXY*resXY;
  fHitCov[2][2] = resZ*resZ;
  //fHitCov[2][2] = 0.0;
  //GFDetPlane d;


  // used for the check with vacuum
  fHitCoord[0][0] = gRandom->Gaus(point.X(),resXY);
  fHitCoord[1][0] = gRandom->Gaus(point.Y(),resXY);
  fHitCoord[2][0] = gRandom->Gaus(point.Z(),resZ);
  //printf("fHitCoord: %lf %lf %lf\n",fHitCoord[0][0],fHitCoord[1][0],fHitCoord[2][0]);
  //fHitCoord[2][0] = point.Z();
  
  //fHitCoord[0][0] = point.X();
  //fHitCoord[1][0] = point.Y();
  //fHitCoord[2][0] = point.Z();



}
GFAbsRecoHit* 
PointHit::clone(){
  return new PointHit(*this);
}


TMatrixT<double>
PointHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ((dynamic_cast<const RKTrackRep*>(stateVector) != NULL)) {
    //I know, since this is the same everytime, it could be done in the
    //the constructor, but I do it here anyway, to make clear that in the
    //case of several track-reps per hit, it would have to be done here
   //    fHMatrix.ResizeTo(NparHitRep,5);
   TMatrixT<double> HMatrix(2,5);

    HMatrix[0][0] = 0.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 1.;
    HMatrix[0][4] = 0.;

    HMatrix[1][0] = 0.;
    HMatrix[1][1] = 0.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    HMatrix[1][4] = 1.;

    return HMatrix;
  }
  else if ((dynamic_cast<const GeaneTrackRep2*>(stateVector) != NULL)) {
    //I know, since this is the same everytime, it could be done in the
    //the constructor, but I do it here anyway, to make clear that in the
    //case of several track-reps per hit, it would have to be done here
   //    fHMatrix.ResizeTo(NparHitRep,5);
   TMatrixT<double> HMatrix(2,6);

    HMatrix[0][0] = 0.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 1.;
    HMatrix[0][4] = 0.;
    HMatrix[0][5] = 0.;

    HMatrix[1][0] = 0.;
    HMatrix[1][1] = 0.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    HMatrix[1][4] = 1.;
    HMatrix[1][5] = 0.;

    return HMatrix;
  }
 else {
   std::cerr << "PointHit can only handle state"
	     << " vectors of type GeaneTrackRep2 -> abort" << std::endl;
   throw;
 }
 
}


