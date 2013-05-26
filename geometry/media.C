{
  TGeoElementTable *table = gGeoManager->GetElementTable();

  TGeoElement* elemH  = table->GetElement(1.0);
  TGeoElement* elemHe = table->GetElement(2.0);
  TGeoElement* elemC  = table->GetElement(6.0);
  TGeoElement* elemN  = table->GetElement(7.0);
  TGeoElement* elemO  = table->GetElement(8.0);
  TGeoElement* elemAl = table->GetElement(13.0);
  TGeoElement* elemAr = table->GetElement(18.0);
  TGeoElement* elemFe = table->GetElement(26.0);

  // mass = g/mol
  double massH = elemH->A();
  double massHe = elemHe->A();
  double massC = elemC->A();
  double massN = elemN->A();
  double massO = elemO->A();
  double massAl = elemAl->A();
  double massAr = elemAr->A();

  // Vaccum
  TGeoMixture *vacuumMat = new TGeoMixture("vacuumMat",3,1.2e-15);
  vacuumMat->AddElement(elemN, 78);
  vacuumMat->AddElement(elemO, 21);
  vacuumMat->AddElement(elemAr, 1);

  // Air
  TGeoMixture *airMat = new TGeoMixture("airMat",3,0.0012); // g/cm3
  airMat->AddElement(elemN, 78);
  airMat->AddElement(elemO, 21);
  airMat->AddElement(elemAr, 1);

  // Iron
  TGeoMaterial *ironMat = new TGeoMaterial("ironMat",elemFe,7.874);

  // Scintillator
  TGeoMixture *scintiMat = new TGeoMixture("scintiMat",2,1.032);
  scintiMat->AddElement(elemC, 8);
  scintiMat->AddElement(elemH, 8);

  // Aluminum
  TGeoMaterial *AlMat = new TGeoMaterial("AlMat",elemAl,2.70);

  // Mylar
  TGeoMixture *mylarMat = new TGeoMixture("mylarMat",3,1.39); // PDG p104
  mylarMat->AddElement(elemH, 8);
  mylarMat->AddElement(elemC, 10);
  mylarMat->AddElement(elemO, 4);

  // Gas ....
  double rho_ethane = 1.356e-3;
  double rho_Ar = 1.782e-3;
  double rho_He = 0.1786e-3;
  double rho_Ar_ethane_50_50 = (rho_Ar+rho_ethane)/2.0;
  double rho_He_ethane_50_50 = (rho_He+rho_ethane)/2.0;

  // Ethane
  TGeoMixture *ethaneMat = new TGeoMixture("ethaneMat",2,rho_ethane);
  ethaneMat->AddElement(elemH, 6);
  ethaneMat->AddElement(elemC, 2);
  
  // Ar
  TGeoMaterial *ArMat = new TGeoMaterial("ArMat",elemAr,rho_Ar);

  // He
  TGeoMaterial *HeMat = new TGeoMaterial("HeMat",elemHe,rho_He);
  
  // CFRP (CDC cylinder material)
  TGeoMaterial *CFRPMat = new TGeoMaterial("CFRP",elemC,1.6);

  // Ar + Ethane
  TGeoMixture *ArEthaneMat = new TGeoMixture("ARGON-ETHANE(50,50)",2, rho_Ar_ethane_50_50);
  ArEthaneMat->AddElement(ArMat,    0.5);
  ArEthaneMat->AddElement(ethaneMat,0.5);

  // He + Ethane
  TGeoMixture *HeEthaneMat = new TGeoMixture("ARGON-ETHANE(50,50)",2, rho_He_ethane_50_50);
  HeEthaneMat->AddElement(HeMat,    0.5);
  HeEthaneMat->AddElement(ethaneMat,0.5);


  /* ****************
   *
   *  Medium 
   *
   *****************/

  unsigned int medInd(0);
  Double_t mPar[10];
  //TAG sollte wieder 0 werden sens flag
  mPar[0]=0.;//sensitive volume flag

  mPar[1]=1.;//magnetic field flag
  //mPar[1]=2.;
  //mPar[1]=0.; // this works

  mPar[2]=30.;//max fiel in kGauss

  mPar[3]=18.0;//maximal angular dev. due to field <= 20 degree is the maxium defined by Geant3
  //mPar[3]=360.0;//maximal angular dev. due to field
  //mPar[3]=1.1;//maximal angular dev. due to field
  //mPar[3]=0.1;//maximal angular dev. due to field
  //mPar[3]=0.001;//maximal angular dev. due to field

  mPar[4]=101.0;//max step allowed (in cm)
  //mPar[4]=1.e-3 ;//max step allowed (in cm)
  // 1mm=1000um = 0.1cm
  // 10um = 0.001
  
  // used by 2012-02-28
  //mPar[4]=0.01;//max step allowed (in cm)

  //mPar[4]=0.001;//max step allowed (in cm) // not changed so much

  mPar[5]=1.e-5;//max fractional energy loss
  mPar[6]=1.e-3;//boundary crossing accuracy
  mPar[7]=1.e-5;//minimum step
  mPar[8]=0.;//not defined
  mPar[9]=0.;//not defined


  TGeoMedium *vacuum   = new TGeoMedium("vacuum",medInd++,vacuumMat,mPar);  
  TGeoMedium *air      = new TGeoMedium("air",   medInd++,airMat,   mPar);  
  TGeoMedium *iron     = new TGeoMedium("iron",  medInd++,ironMat,mPar);
  TGeoMedium *scinti   = new TGeoMedium("scinti",medInd++,scintiMat,mPar);  
  TGeoMedium *alumi    = new TGeoMedium("alumi", medInd++,AlMat,    mPar);  
  TGeoMedium *mylar    = new TGeoMedium("mylar", medInd++,mylarMat, mPar);  
  TGeoMedium *CFRP     = new TGeoMedium("CFRP",  medInd++,CFRPMat,mPar);  
  TGeoMedium *ArEthane = new TGeoMedium("ARGON-ETHANE(50,50)",medInd++,ArEthaneMat,mPar);
  TGeoMedium *HeGas    = new TGeoMedium("HeGas",medInd++,HeMat,mPar);  
  TGeoMedium *HeEthane = new TGeoMedium("HELIUM-ETHANE(50,50)",medInd++,HeEthaneMat,mPar);

}
