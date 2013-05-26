#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGeoManager.h"
#include "TPluginManager.h"

#include "GFException.h"
#include "GFAbsTrackRep.h"
#include "GeaneTrackRep2.h" // <= using TGeant3
#include "RKTrackRep.h"

#include "GFConstField.h"
#include "GFFieldManager.h"

#include "PointHit.h"
#include "StripHit.h"
#include "GFTrack.h"
#include "GFKalman.h"
#include "GFDaf.h"

#include "TGeoManager.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h" // TGeoTranslation

#include "unistd.h" //getopt
#include "stdlib.h"
/* arguments */
int arg_seed;
char arg_genfitGeom_name[128];
char arg_output_root[128];
int arg_total;
char arg_event_type[128];
char arg_hit_type[128];
char arg_track_type[128];
char arg_extrap_pos[128];
char arg_input_root[128];
FILE* arg_writeout_track_fp;
double arg_posini_x;
double arg_posini_y;
double arg_posini_z;
double arg_momini_x;
double arg_momini_y;
double arg_momini_z;
/* end of arguments */

/* config */

char g_inwall_material[32];
double g_inwall_first_pos;
double g_inwall_thickness;

char g_chamber_type[32];
char g_chamber_gas[32];
double g_chamber_first_pos;
double g_chamber_last_pos;
double g_chamber_spacing;
int g_chamber_num_layer;

char g_scinti_type[32];
char g_scinti_material[32];
double g_scinti_thickness;
double g_scinti_first_pos;

char g_target_type[32];
char g_target_material[32];
double g_target_thickness;
double g_target_offset_z;
double g_target_radius;
double g_target_spacing;
int g_target_number_of_disk;
int g_target_center_disk_number;

char g_solenoid_material[32];
char g_solenoid_bfld_type[32]; // uniform, gradient
double g_solenoid_bfld_tesla;

double g_dio_startpoint;// MeV/c
/* end configure */


TMatrixT<double> *stMCT;
TMatrixT<double> *covMCT;
TMatrixT<double> *stREC;
TMatrixT<double> *covREC;

//TMatrixT<double> *stREC;
//TMatrixT<double> *covREC;

//int PDGcode = 211; // pi+
//int PDGcode = -211; // pi-
int PDGcode = 11; // e-
//int PDGcode = -11; // e+
//int PDGcode = 2212; // proton

inline double sqrt2(double a, double b) { return sqrt(a*a+b*b); }
inline double sqrt3(double a, double b, double c) { return sqrt(a*a+b*b+c*c); }

int g_nhits_fit;



#define MAX_HIT 1500
struct hits
{
   int nhits;
   int layer[MAX_HIT];
   double posx[MAX_HIT];
   double posy[MAX_HIT];
   double posz[MAX_HIT];
   double momx[MAX_HIT];
   double momy[MAX_HIT];
   double momz[MAX_HIT];
   double tofg[MAX_HIT];
};
void print_hits(int idx, struct hits* hits)
{
   printf("idx %d: layer %d (pos) %lf %lf %lf (mom) %lf %lf %lf (tofg) %lf\n",
         idx,
         hits->layer[idx],
         hits->posx[idx],
         hits->posy[idx],
         hits->posz[idx],
         hits->momx[idx],
         hits->momy[idx],
         hits->momz[idx],
         hits->tofg[idx]);
}

struct hits g_hits_biw;
struct hits g_hits_det;
struct hits g_hits_abs;
struct hits g_hits_tgt;
int g_nhits_tgt_before_chamber;


//pulls
double invmom;
double sigmasqustate;
double sigma_p;
double momSi;
double momRe;
double momTr;
double momPu;

double qopSi;
double qopRe;
double qopTr;
double qopPu;

double upSi;
double upRe;
double upTr;
double upPu;
double vpSi;
double vpRe;
double vpTr;
double vpPu;

double uSi;
double uRe;
double uTr;
double uPu;
double vSi;
double vRe;
double vTr;
double vPu;

struct trackpar
{
   double x;
   double y;
   double z;
   double px;
   double py;
   double pz; // longitudinal momentum
   //double pt; // transverse momentum
   //double pa; // all(total) momentum
};
struct tree_value
{
   int iev;
   int error;
   int nfail;
   int ndf;
   double prob;
   double chi2;
   int ini_disk_number;
   double elapsed_time; // second
   struct trackpar ini;  // hit at start point
   struct trackpar hit;  // first hit
   struct trackpar fit;  // extrapolated at first hit
};
void print_tree_value(const char* prefix, FILE*fp, struct tree_value* tv)
{
   fprintf(fp,"=============\n");
   fprintf(fp,"%s\n",prefix);
   fprintf(fp,"=============\n");
   fprintf(fp,"iev %d\n",tv->iev);
   fprintf(fp,"error %d\n",tv->error);
   fprintf(fp,"nfail %d\n",tv->nfail);
   fprintf(fp,"ndf %d\n",tv->ndf);

   //fprintf(fp,"nhits %d\n",g_nhits);
   fprintf(fp,"nhits %d\n",g_hits_det.nhits);
   fprintf(fp,"nhits_tgt %d\n",g_hits_tgt.nhits);
   fprintf(fp,"numhits_tgt_before_chamber %d\n",g_nhits_tgt_before_chamber);

   fprintf(fp,"nhits_fit %d\n",g_nhits_fit);
   fprintf(fp,"chi2 %lf\n",tv->chi2);
   fprintf(fp,"prob %lf\n",tv->prob);
   fprintf(fp,"elapsted_time %lf\n",tv->elapsed_time);
   fprintf(fp,"ini_disk_num %d\n",tv->ini_disk_number);
   double ipt = sqrt2(tv->ini.px,tv->ini.py);
   double ipa = sqrt3(tv->ini.px,tv->ini.py,tv->ini.pz);
   fprintf(fp,"ini_(x,y,z)    %10.5f %10.5f %10.5f\n",tv->ini.x,tv->ini.y,tv->ini.z);
   fprintf(fp,"ini_(px,py,pz) %10.5f %10.5f %10.5f : pt %10.5lf pa %10.5lf\n",tv->ini.px,tv->ini.py,tv->ini.pz, ipt, ipa);
   fprintf(fp,"x(cm)    :  %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.x,  tv->fit.x  , tv->fit.x - tv->hit.x);
   fprintf(fp,"y(cm)    :  %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.y,  tv->fit.y  , tv->fit.y - tv->hit.y);
   fprintf(fp,"z(cm)    :  %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.z,  tv->fit.z  , tv->fit.z - tv->hit.z);
   fprintf(fp,"px(GeV/c):  %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.px, tv->fit.px,  tv->fit.px- tv->hit.px);
   fprintf(fp,"py(GeV/c):  %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.py, tv->fit.py,  tv->fit.py- tv->hit.py);
   fprintf(fp,"pz(GeV/c):  %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.pz, tv->fit.pz,  tv->fit.pz- tv->hit.pz);
   double hpt = sqrt2(tv->hit.px,tv->hit.py);
   double hpa = sqrt3(tv->hit.px,tv->hit.py,tv->hit.pz);
   double fpt = sqrt2(tv->fit.px,tv->fit.py);
   double fpa = sqrt3(tv->fit.px,tv->fit.py,tv->fit.pz);
   fprintf(fp,"pt(GeV/c):  %10.5f => %10.5f : diff(%10.5f)\n",hpt, fpt,  fpt- hpt);
   fprintf(fp,"pa(GeV/c):  %10.5f => %10.5f : diff(%10.5f)\n",hpa, fpa,  fpa- hpa);
   fprintf(fp,"=============\n");
}
void fill_tree_config(const char* tname)
{
   TTree* t = new TTree(tname,tname);
   /* arguments */
   t->Branch("seed",&arg_seed,"seed/I");
   t->Branch("event_type",arg_event_type,"event_type/C");
   t->Branch("hit_type",arg_hit_type,"hit_type/C");
   t->Branch("track_type",arg_track_type,"track_type/C");
   t->Branch("extrap_pos",arg_extrap_pos,"extrap_pos/C");

   /* contents of config.txt */
   t->Branch("inwall_material",g_inwall_material,"inwall_material/C");
   t->Branch("inwall_first_pos",&g_inwall_first_pos,"inwall_first_pos/D");
   t->Branch("inwall_thickness",&g_inwall_thickness,"inwall_thickness/D");

   t->Branch("chamber_type"      ,g_chamber_type,"chamber_type/C");
   t->Branch("chamber_gas"       ,g_chamber_gas, "chamber_gas/C");
   t->Branch("chamber_first_pos" ,&g_chamber_first_pos,"chamber_first_pos/D");
   t->Branch("chamber_last_pos"  ,&g_chamber_last_pos,"chamber_last_pos/D");
   t->Branch("chamber_spacing"   ,&g_chamber_spacing, "chamber_spacing/D");
   t->Branch("chamber_num_layer" ,&g_chamber_num_layer,"chamber_num_layer/I");

   t->Branch("scinti_type",g_scinti_type,"g_scinti_type/C");
   t->Branch("scinti_material",g_scinti_material,"g_scinti_material/C");
   t->Branch("scinti_thickness",&g_scinti_thickness,"g_scinti_thickness/D");
   t->Branch("scinti_first_pos",&g_scinti_first_pos,"g_scinti_first_pos/D");

   t->Branch("target_type",g_target_type,"g_target_type/C");
   t->Branch("target_material",g_target_material,"g_target_material/C");
   t->Branch("target_thickness",&g_target_thickness,"g_target_thickness/D");
   t->Branch("target_offset_z",&g_target_offset_z,"g_target_offset_z/D");
   t->Branch("target_radius",&g_target_radius,"g_target_radius/D");
   t->Branch("target_spacing",&g_target_spacing,"g_target_spacing/D");

   t->Branch("solenoid_material",g_solenoid_material,"g_solenoid_material/C");
   t->Branch("solenoid_bfld_type",g_solenoid_bfld_type,"g_solenoid_bfld_type/C");
   t->Branch("solenoid_bfld_tesla",&g_solenoid_bfld_tesla,"g_solenoid_bfld_tesla/D");

   t->Fill();
   t->Write();
}
TTree* make_branch(const char* tname, struct tree_value* tv)
{
   TTree* t = new TTree(tname,tname);
   //t->Branch("stMCT","TMatrixT<double>",&stMCT);
   //t->Branch("covMCT","TMatrixT<double>",&covMCT);
   //t->Branch("stREC","TMatrixT<double>",&stREC);
   //t->Branch("covREC","TMatrixT<double>",&covREC);
   t->Branch("iev",&tv->iev,"iev/I");
   t->Branch("error",&tv->error,"error/I");
   t->Branch("nfail",&tv->nfail,"nfail/I");
   t->Branch("ndf",&tv->ndf,"ndf/I");
   t->Branch("chi2",&tv->chi2,"chi2/D");
   t->Branch("prob",&tv->prob,"prob/D"); // chi2-probability
   t->Branch("nhits_fit",&g_nhits_fit,"nhits_fit/I"); // number of (first series of) hits used for fitting_

   // time for making hits
   t->Branch("elapsed_time",&tv->elapsed_time,"elapsed_time/D");

   // start disk number
   t->Branch("ini_disk_number",&tv->ini_disk_number,"ini_disk_number/I");

   // number of hits at target (Al disk etc.) before reaching at chamber
   t->Branch("nhits_tgt_before_chamber",&g_nhits_tgt_before_chamber,"nhits_tgt_before_chamber/I");

   // hit at start point
   t->Branch("ini_x",&tv->ini.x,"ini_x/D");
   t->Branch("ini_y",&tv->ini.y,"ini_y/D");
   t->Branch("ini_z",&tv->ini.z,"ini_z/D");
   t->Branch("ini_px",&tv->ini.px,"ini_px/D");
   t->Branch("ini_py",&tv->ini.py,"ini_py/D");
   t->Branch("ini_pz",&tv->ini.pz,"ini_pz/D");
   t->SetAlias("ini_pt","sqrt(ini_px*ini_px+ini_py*ini_py)");
   t->SetAlias("ini_pa","sqrt(ini_pt*ini_pt+ini_pz*ini_pz)");
   t->SetAlias("ini_r","sqrt(ini_x*ini_x+ini_y*ini_y)");

   // first hit
   t->Branch("hit_x",&tv->hit.x,"hit_x/D");
   t->Branch("hit_y",&tv->hit.y,"hit_y/D");
   t->Branch("hit_z",&tv->hit.z,"hit_z/D");
   t->Branch("hit_px",&tv->hit.px,"hit_px/D");
   t->Branch("hit_py",&tv->hit.py,"hit_py/D");
   t->Branch("hit_pz",&tv->hit.pz,"hit_pz/D");
   t->SetAlias("hit_pt","sqrt(hit_px*hit_px+hit_py*hit_py)");
   t->SetAlias("hit_pa","sqrt(hit_pt*hit_pt+hit_pz*hit_pz)");
   t->SetAlias("hit_r","sqrt(hit_x*hit_x+hit_y*hit_y)");

   // extraplated at first hit
   t->Branch("fit_x",&tv->fit.x,"fit_x/D");
   t->Branch("fit_y",&tv->fit.y,"fit_y/D");
   t->Branch("fit_z",&tv->fit.z,"fit_z/D");
   t->Branch("fit_px",&tv->fit.px,"fit_px/D");
   t->Branch("fit_py",&tv->fit.py,"fit_py/D");
   t->Branch("fit_pz",&tv->fit.pz,"fit_pz/D");
   t->SetAlias("fit_pt","sqrt(fit_px*fit_px+fit_py*fit_py)");
   t->SetAlias("fit_pa","sqrt(fit_pt*fit_pt+fit_pz*fit_pz)");
   t->SetAlias("fit_r","sqrt(fit_x*fit_x+fit_y*fit_y)");

   // Residuals
   t->SetAlias("res_x","(hit_x-fit_x)");
   t->SetAlias("res_y","(hit_y-fit_y)");
   t->SetAlias("res_r","(hit_r-fit_r)");
   t->SetAlias("res_px","(hit_px-fit_px)");
   t->SetAlias("res_py","(hit_py-fit_py)");
   t->SetAlias("res_pz","(hit_pz-fit_pz)");
   t->SetAlias("res_pt","(hit_pt-fit_pt)");
   t->SetAlias("res_pa","(hit_pa-fit_pa)");

   // hits at target
   t->Branch("nhits_tgt",&g_hits_tgt.nhits,"nhits_tgt/I");
   t->Branch("hits_tgt_disk_number",g_hits_tgt.layer,"hits_tgt_disk_number[nhits_tgt]/I");
   t->Branch("hits_tgt_x",g_hits_tgt.posx,"hits_tgt_x[nhits_tgt]/D");
   t->Branch("hits_tgt_y",g_hits_tgt.posy,"hits_tgt_y[nhits_tgt]/D");
   t->Branch("hits_tgt_z",g_hits_tgt.posz,"hits_tgt_z[nhits_tgt]/D");
   t->Branch("hits_tgt_px",g_hits_tgt.momx,"hits_tgt_px[nhits_tgt]/D");
   t->Branch("hits_tgt_py",g_hits_tgt.momy,"hits_tgt_py[nhits_tgt]/D");
   t->Branch("hits_tgt_pz",g_hits_tgt.momz,"hits_tgt_pz[nhits_tgt]/D");
   t->Branch("hits_tgt_tofg",g_hits_tgt.tofg,"hits_tgt_tofg[nhits_tgt]/D");
   t->SetAlias("hits_tgt_pt","sqrt(hits_tgt_px*hits_tgt_px+hits_tgt_py*hits_tgt_py)");
   t->SetAlias("hits_tgt_pa","sqrt(hits_tgt_pt*hits_tgt_pt+hits_tgt_pz*hits_tgt_pz)");
   t->SetAlias("hits_tgt_r","sqrt(hits_tgt_x*hits_tgt_x+hits_tgt_y*hits_tgt_y)");

   // hits at absorber (scinti_layer)
   t->Branch("nhits_abs",&g_hits_abs.nhits,"nhits_abs/I");
   t->Branch("hits_abs_x",g_hits_abs.posx,"hits_abs_x[nhits_abs]/D");
   t->Branch("hits_abs_y",g_hits_abs.posy,"hits_abs_y[nhits_abs]/D");
   t->Branch("hits_abs_z",g_hits_abs.posz,"hits_abs_z[nhits_abs]/D");
   t->Branch("hits_abs_px",g_hits_abs.momx,"hits_abs_px[nhits_abs]/D");
   t->Branch("hits_abs_py",g_hits_abs.momy,"hits_abs_py[nhits_abs]/D");
   t->Branch("hits_abs_pz",g_hits_abs.momz,"hits_abs_pz[nhits_abs]/D");
   t->Branch("hits_abs_tofg",g_hits_abs.tofg,"hits_abs_tofg[nhits_abs]/D");
   t->SetAlias("hits_abs_pt","sqrt(hits_abs_px*hits_abs_px+hits_abs_py*hits_abs_py)");
   t->SetAlias("hits_abs_pa","sqrt(hits_abs_pt*hits_abs_pt+hits_abs_pz*hits_abs_pz)");
   t->SetAlias("hits_abs_r","sqrt(hits_abs_x*hits_abs_x+hits_abs_y*hits_abs_y)");

   // hits at detector (straw tubes or drift chamber)
   t->Branch("nhits",&g_hits_det.nhits,"nhits/I");
   t->Branch("hits_layer",g_hits_det.layer,"hits_layer[nhits]/I");
   t->Branch("hits_x",g_hits_det.posx,"hits_x[nhits]/D");
   t->Branch("hits_y",g_hits_det.posy,"hits_y[nhits]/D");
   t->Branch("hits_z",g_hits_det.posz,"hits_z[nhits]/D");
   t->Branch("hits_px",g_hits_det.momx,"hits_px[nhits]/D");
   t->Branch("hits_py",g_hits_det.momy,"hits_py[nhits]/D");
   t->Branch("hits_pz",g_hits_det.momz,"hits_pz[nhits]/D");
   t->Branch("hits_tofg",g_hits_det.tofg,"hits_tofg[nhits]/D");
   t->SetAlias("hits_pt","sqrt(hits_px*hits_px+hits_py*hits_py)");
   t->SetAlias("hits_pa","sqrt(hits_pt*hits_pt+hits_pz*hits_pz)");
   t->SetAlias("hits_r","sqrt(hits_x*hits_x+hits_y*hits_y)");

   // hits at virtual plane before Inwall
   t->Branch("nhits_biw",&g_hits_biw.nhits,"nhits_biw/I");
   t->Branch("hits_biw_x",g_hits_biw.posx,"hits_biw_x[nhits_biw]/D");
   t->Branch("hits_biw_y",g_hits_biw.posy,"hits_biw_y[nhits_biw]/D");
   t->Branch("hits_biw_z",g_hits_biw.posz,"hits_biw_z[nhits_biw]/D");
   t->Branch("hits_biw_px",g_hits_biw.momx,"hits_biw_px[nhits_biw]/D");
   t->Branch("hits_biw_py",g_hits_biw.momy,"hits_biw_py[nhits_biw]/D");
   t->Branch("hits_biw_pz",g_hits_biw.momz,"hits_biw_pz[nhits_biw]/D");
   t->Branch("hits_biw_tofg",g_hits_biw.tofg,"hits_biw_tofg[nhits_biw]/D");
   t->SetAlias("hits_biw_pt","sqrt(hits_biw_px*hits_biw_px+hits_biw_py*hits_biw_py)");
   t->SetAlias("hits_biw_pa","sqrt(hits_biw_pt*hits_biw_pt+hits_biw_pz*hits_biw_pz)");
   t->SetAlias("hits_biw_r", "sqrt(hits_biw_x*hits_biw_x+hits_biw_y*hits_biw_y)");


   // pulls
   t->Branch("invmom",&invmom,"invmom/D");
   t->Branch("sigmasqustate",&sigmasqustate,"sigmasqustate/D");
   t->Branch("sigma_p",&sigma_p,"sigma_p/D");
   t->Branch("momSi",&momSi,"momSi/D");
   t->Branch("momRe",&momRe,"momRe/D");
   t->Branch("momTr",&momTr,"momTr/D");
   t->Branch("momPu",&momPu,"momPu/D");

   t->Branch("qopSi",&qopSi,"qopSi/D");
   t->Branch("qopRe",&qopRe,"qopRe/D");
   t->Branch("qopTr",&qopTr,"qopTr/D");
   t->Branch("qopPu",&qopPu,"qopPu/D");

   t->Branch("upSi",&upSi,"upSi/D");
   t->Branch("upRe",&upRe,"upRe/D");
   t->Branch("upTr",&upTr,"upTr/D");
   t->Branch("upPu",&upPu,"upPu/D");
   t->Branch("vpSi",&vpSi,"vpSi/D");
   t->Branch("vpRe",&vpRe,"vpRe/D");
   t->Branch("vpTr",&vpTr,"vpTr/D");
   t->Branch("vpPu",&vpPu,"vpPu/D");

   t->Branch("upSi",&uSi,"uSi/D");
   t->Branch("upRe",&uRe,"uRe/D");
   t->Branch("upTr",&uTr,"uTr/D");
   t->Branch("upPu",&uPu,"uPu/D");
   t->Branch("vpSi",&vSi,"vSi/D");
   t->Branch("vpRe",&vRe,"vRe/D");
   t->Branch("vpTr",&vTr,"vTr/D");
   t->Branch("vpPu",&vPu,"vPu/D");
   return t;

}

static int first_call;
static TF1* f1;
int get_init_pos_cone(TVector3& pos)
{
   /*
    this assumes angle 45 and z = 10 cm.
    Need to modify if test with other version.
    */


   // unifrom in circle (R=10)
   double xpos,ypos,zpos;
   double rpos;
   while (1) {
      xpos = gRandom->Uniform(-g_target_radius,g_target_radius);
      ypos = gRandom->Uniform(-g_target_radius,g_target_radius);
      rpos = sqrt(xpos*xpos+ypos*ypos);
      if (rpos>=g_target_thickness/1.4142135623 && rpos<=g_target_radius) break;
   }
   // determin depth (z-axis)
   double zmin,zmax;
   int disk_number = (int)gRandom->Integer(10) +1; // [1,10]
   if (disk_number%2==1) {
      zmax = rpos;
      zmin = zmax - g_target_thickness/1.4142135623;
   } else {
      zmin = 10.0-rpos;
      zmax = zmin + g_target_thickness/1.4142135623;
   }
   zpos = gRandom->Uniform(zmin,zmax);
   zpos += (disk_number-5)*10.0 - 10.0;
   pos.SetXYZ(xpos,ypos,zpos);
   //printf("##disk_number %d xpos %lf ypos %lf zpos %lf\n",disk_number,xpos,ypos,zpos);

   return disk_number;
}
int get_init_pos_cylinder(TVector3& pos)
{
   double zmin = -40.0; 
   double zmax = +40.0;
   double rmin = 5.0; //5cm
   double thickness_cylinder = g_target_thickness;//0.004; // 40um
   double rmax = rmin + thickness_cylinder;

   if (first_call==0) {
      f1 = new TF1("f1","x*x/2.0",rmin,rmax);
      first_call=1;
   }

   // 10000um = 10mm = 1cm
   // 10um = 0.001cm

   double rpos = f1->GetRandom();
   double theta = gRandom->Uniform(0.0,2.0*TMath::Pi());
   double xpos = rpos*cos(theta);
   double ypos = rpos*sin(theta);
   double zpos = gRandom->Uniform(zmin,zmax);
   pos.SetXYZ(xpos,ypos,zpos);
   //printf("xpos %lf ypos %lf zpos %lf\n",xpos,ypos,zpos);

   int disk_number = -1;
   return disk_number;
}
int get_init_pos(TVector3& pos)
{
   double xpos,ypos,zpos;
#if 0
   // uniform in x [0,10] cm, 10cm = R of Al disk
   double xpos = gRandom->Uniform(0,10); // cm
   double ypos = 0.0;
   double zpos = 0.0;
#endif
   // fist, select disk number

   // generate 0 to imax-1
   // Integer(imax)
   int number_of_disk = g_target_number_of_disk;
   int disk_number = (int)gRandom->Integer(number_of_disk) +1; // [1,17]

   int center_disk_number = g_target_center_disk_number;
   double spacing = g_target_spacing;//5.0; //cm
   double thickness_disk = g_target_thickness; //0.01; // 100um
   double zmin = -thickness_disk/2.0 + (disk_number-center_disk_number)*spacing  + g_target_offset_z;
   double zmax = +thickness_disk/2.0 + (disk_number-center_disk_number)*spacing  + g_target_offset_z;
   zpos = gRandom->Uniform(zmin,zmax); // cm (thickness of disk is 200 um)

   while (1) {
      xpos = gRandom->Uniform(-g_target_radius,g_target_radius); // cm
      ypos = gRandom->Uniform(-g_target_radius,g_target_radius); // cm
      double rpos = sqrt(xpos*xpos+ypos*ypos);
      if (rpos<=g_target_radius) {
         break;
      }
   }
   pos.SetXYZ(xpos,ypos,zpos);
   return disk_number;
}
void get_init_mom_test(TVector3& mom)
{
   double ptotal = 0.105; //GeV/c
   /* fail
      double xmom = 0.0994145;
      double ymom = 0.0117207;
      double zmom = 0.0316918;
      */
   /* fail
      double xmom = 0.0117207;
      double ymom = 0.0994145;
      double zmom = 0.0316918;
      */
   /* same!
      double xmom = 0.05;
      double ymom = 0.05;
      double zmom = sqrt(0.105*0.105-xmom*xmom-ymom*ymom);
      */
   /* same
      double xmom = 0.06;
      double ymom = 0.03;
      double zmom = sqrt(0.105*0.105-xmom*xmom-ymom*ymom);
      */
   /* fail
      double xmom = 0.09;
      double ymom = 0.01;
      double zmom = sqrt(0.105*0.105-xmom*xmom-ymom*ymom);
      */
   /* fail
      double xmom = 0.09;
      double ymom = 0.00;
      double zmom = sqrt(0.105*0.105-xmom*xmom-ymom*ymom);
      */
   /* fail
      double xmom = 0.00;
      double ymom = 0.09;
      double zmom = sqrt(0.105*0.105-xmom*xmom-ymom*ymom);
      */
   /* ok
      double xmom = 0.00;
      double ymom = 0.08;
      double zmom = sqrt(0.105*0.105-xmom*xmom-ymom*ymom);
      */
   double ymom = 0.0;
   //double xmom=0.1049;
   double xmom=0.1043;
   //double zmom = 0.075; // OK
   //double zmom = 0.073;
   double zmom = sqrt(0.105*0.105-xmom*xmom-zmom*zmom);
   mom.SetXYZ(xmom,ymom,zmom);
}
/* Spectrum is cited from Physical Review D 84, 013006 (2011) */
void get_init_mom_DIO(TVector3& mom)
{
   //double e_startpoint = 0.103; // GeV/c
   double e_startpoint = g_dio_startpoint / 1000.0; // GeV/c
   double e_endpoint   = 0.104973; // GeV/c

   TF1* f1 = new TF1("watanabe_shanker","[1]*pow([0]-x,5)+[2]*pow([0]-x,6)+[3]*pow([0]-x,7)+[4]*pow([0]-x,8)",e_startpoint,e_endpoint);
   f1->SetParameter(0, e_endpoint);
   f1->SetParameter(1, 8.6434e-17);
   f1->SetParameter(2,1.16874e-17);
   f1->SetParameter(3,-1.87828e-19);
   f1->SetParameter(4,9.16327e-20);

   double e_mass = 0.00051099892; // GeV/c
   double e_energy;
   /* Skip energy below e_startpoint to save simulation time */
   while (1) {
      e_energy = f1->GetRandom();
      if (e_energy>=e_startpoint && e_energy<=e_endpoint) {
         break;
      }
   }
   double ptotal = sqrt(e_energy*e_energy-e_mass*e_mass);
   double xmom;
   double ymom;
   double zmom;
#if 0
   // emit 2pi (forword)
   while (1) {
      gRandom->Sphere(xmom,ymom,zmom,ptotal);
      //if (zmom>0) break;
      if (zmom>0 && (xmom!=0&&ymom!=0)) break;
   }
#endif
   // emit 4pi
   gRandom->Sphere(xmom,ymom,zmom,ptotal);
   mom.SetXYZ(xmom,ymom,zmom);
}
/* Proton emission spectrum (Silicon data)
 * Paper: Comments on Proton Emission after Muon Capture (meco034.pdf)
 * */
void get_init_mom_Proton(TVector3& mom)
{
   // x = Momentum(MeV)
   // K = E - M = sqrt(P^2+M^2) - M
   // P = 51.275 (MeV/c) @ K = 1.4 MeV
   TF1* f2 = new TF1("f2","[0]*pow(1-[1]/(sqrt(x*x+[4]*[4])-[4]),[2])*exp(-(sqrt(x*x+[4]*[4])-[4])/[3])",51.3,300);
   f2->SetParameter(0, 0.105); // A (MeV-1)
   f2->SetParameter(1, 1.4); // T_th (MeV)
   f2->SetParameter(2,1.3279); // alpha
   f2->SetParameter(3,3.1); // T_0
   f2->SetParameter(4,938.272); // Proton Mass (MeV)

   double ptotal = f2->GetRandom() * 0.001; // MeV -> GeV
   //   double ptotal = 100; 

   double xmom;
   double ymom;
   double zmom;
#if 0
   while (1) {
      gRandom->Sphere(xmom,ymom,zmom,ptotal);
      //if (zmom>0) break;
      if (zmom>0 && (xmom!=0&&ymom!=0)) break;
   }
#endif
   gRandom->Sphere(xmom,ymom,zmom,ptotal);
   mom.SetXYZ(xmom,ymom,zmom);
}
void get_init_mom(TVector3& mom)
{
   // magnetide is 105 MeV/c
   // direction is isotropically, but only forward direction
   //double ptotal = 0.105; //GeV/c
   double ptotal = 0.104973; //GeV/c
   double xmom;
   double ymom;
   double zmom;
#if 0
   while (1) {
      gRandom->Sphere(xmom,ymom,zmom,ptotal);
      //if (zmom>0) break;
      if (zmom>0 && (xmom!=0&&ymom!=0)) break;
   }
#endif
   gRandom->Sphere(xmom,ymom,zmom,ptotal);
   mom.SetXYZ(xmom,ymom,zmom);
}
int root_opened=0;
double root_ini_posx;
double root_ini_posy;
double root_ini_posz;
double root_ini_momx;
double root_ini_momy;
double root_ini_momz;
int root_nhits;
double root_hits_posx[MAX_HIT];
double root_hits_posy[MAX_HIT];
double root_hits_posz[MAX_HIT];
double root_hits_momx[MAX_HIT];
double root_hits_momy[MAX_HIT];
double root_hits_momz[MAX_HIT];
TFile* root_f;
TTree* root_t;
int read_hit(int iev, char* root_file, TVector3& posini, TVector3& momini)
{
   if (root_opened==0) {
      if (root_file==NULL)
         exit(1);
      root_f = new TFile(root_file);
      root_t = (TTree*)root_f->Get("t");
      root_t->SetBranchAddress("nhits",&root_nhits);
      root_t->SetBranchAddress("hit_x",&root_ini_posx);
      root_t->SetBranchAddress("hit_y",&root_ini_posy);
      root_t->SetBranchAddress("hit_z",&root_ini_posz);
      root_t->SetBranchAddress("hit_px",&root_ini_momx);
      root_t->SetBranchAddress("hit_py",&root_ini_momy);
      root_t->SetBranchAddress("hit_pz",&root_ini_momz);
      root_t->SetBranchAddress("hits_posx",root_hits_posx);
      root_t->SetBranchAddress("hits_posy",root_hits_posy);
      root_t->SetBranchAddress("hits_posz",root_hits_posz);
      root_t->SetBranchAddress("hits_momx",root_hits_momx);
      root_t->SetBranchAddress("hits_momy",root_hits_momy);
      root_t->SetBranchAddress("hits_momz",root_hits_momz);
      root_opened = 1;
   }

   int ret = root_t->GetEntry(iev);
   if (ret == 0) { // no more events
      return -1;
   }
   posini.SetX(root_ini_posx);
   posini.SetY(root_ini_posy);
   posini.SetZ(root_ini_posz);
   momini.SetX(root_ini_momx);
   momini.SetY(root_ini_momy);
   momini.SetZ(root_ini_momz);

   /* use RKTrackRep just to calc stMCT,covMCT */
   GFAbsTrackRep* rephits = new RKTrackRep( posini, momini, PDGcode);
   stMCT->ResizeTo(rephits->getState());
   *stMCT = rephits->getState();
   covMCT->ResizeTo(rephits->getCov());
   *covMCT = rephits->getCov();

   printf("root_nhits %d\n",root_nhits);
   g_hits_det.nhits = 0;
   for (int i=0; i<root_nhits; i++) {
      g_hits_det.posx[g_hits_det.nhits] = root_hits_posx[i];
      g_hits_det.posy[g_hits_det.nhits] = root_hits_posy[i];
      g_hits_det.posz[g_hits_det.nhits] = root_hits_posz[i];
      g_hits_det.momx[g_hits_det.nhits] = root_hits_momx[i];
      g_hits_det.momy[g_hits_det.nhits] = root_hits_momy[i];
      g_hits_det.momz[g_hits_det.nhits] = root_hits_momz[i];
      g_hits_det.nhits++;
   }
   return 0;
}
void read_init_param(char* fname)
{
   char line[128];
   FILE* fp = fopen(fname,"r");
   if (fp==NULL) {
      fprintf(stderr,"ERROR: failed to open init_param file %s\n",fname);
      exit(1);
   }
   while (fgets(line,sizeof(line),fp)) {
      //fprintf(stderr,"failed to read initial parameter, no line?\n");
      //exit(1);
      if (line[0]!='#') { // skip this lines
         break;
      }
   } 
   sscanf(line,"%lf %lf %lf %lf %lf %lf",
         &arg_posini_x,
         &arg_posini_y,
         &arg_posini_z,
         &arg_momini_x,
         &arg_momini_y,
         &arg_momini_z);
   printf("init_para: x=%lf y=%lf z=%lf px=%lf py=%lf pz=%lf (pt=%lf)\n",
         arg_posini_x,
         arg_posini_y,
         arg_posini_z,
         arg_momini_x,
         arg_momini_y,
         arg_momini_z,
         sqrt2( arg_momini_x, arg_momini_y)
         );
   fclose(fp);
}
int get_layer_number(const char* volname)
{
   int id;
   int ret = sscanf(volname,"layer%d_7",&id);
   if (ret!=1)
      return -1;

   return id+1;
#if 0
   //printf("volname %s\n",volname);
   if      (strcmp(volname,"layer0_7")==0) { return 1; } 
   else if (strcmp(volname,"layer1_7")==0) { return 2; } 
   else if (strcmp(volname,"layer2_7")==0) { return 3; } 
   else if (strcmp(volname,"layer3_7")==0) { return 4; } 
   else if (strcmp(volname,"layer4_7")==0) { return 5; } 
   else if (strcmp(volname,"layer5_7")==0) { return 6; }
   else { return -1; }
   return -1;
#endif
}
int get_hit_at_target(const char* volname, const char* pre_volname, int* disk_number)
{
   // return 1 if hit at target, otherwise return 0

   *disk_number = 0;
   int disk_idx = -1; // [0..xx]

   //printf("get_hits_at_target pre_volname %s volname %s\n",pre_volname,volname);
   if (strcmp(pre_volname,"TOPPER")==0) {
      if (strcmp(volname,"TOPPER")==0) {
         return 0;
      }
      if (strstr(volname,"disk")!=NULL) {
         sscanf(volname,"disk%d",&disk_idx);
      } else if (strcmp(volname,"cylinder")==0) {
         disk_idx = 0;
      } else if (strstr(volname,"cone")!=NULL) {
         sscanf(volname,"cone%d",&disk_idx);
      } else {
         return 0;
      }
      *disk_number = disk_idx+1; // [1.xx]
      //printf("volname_pre %s volanme %s *disk_number %d disk_idx %d\n",volname,pre_volname,*disk_number,disk_idx);
      return 1;
   }
   //printf("get_hits_at_target return 0\n");
   return 0;
}
int get_hit_at_scinti(const char* volname, const char* pre_volname)
{
   //printf("get_hits_at_scinti pre_volname %s volname %s\n",pre_volname,volname);
   // volume around scinti_layer should be TOP!!
   if (strcmp(pre_volname,"TOPPER")==0 &&
         strcmp(volname,"scinti_layer")==0) {
      return 1;
   }
   return 0;
}
int get_hit_at_biw(const char* volname, const char* pre_volname)
{
   //printf("get_hits_at_biw pre_volname %s volname %s\n",pre_volname,volname);
   // layer before inner wall
   if (strcmp(pre_volname,"TOPPER")==0 &&
         strcmp(volname,"biw_layer")==0) {
      return 1;
   }
   return 0;
}
int make_true_hit(int iev, struct tree_value* tv, 
      TVector3& posini, 
      TVector3& momini, 
      TVector3 &posErr,
      TVector3 &momErr,
      GFAbsTrackRep* rephits)
{

   /* initial hit is recored in target as well */
   g_hits_tgt.layer[g_hits_tgt.nhits] = tv->ini_disk_number;
   g_hits_tgt.posx[g_hits_tgt.nhits] = posini.X();
   g_hits_tgt.posy[g_hits_tgt.nhits] = posini.Y();
   g_hits_tgt.posz[g_hits_tgt.nhits] = posini.Z();
   g_hits_tgt.momx[g_hits_tgt.nhits] = momini.X();
   g_hits_tgt.momy[g_hits_tgt.nhits] = momini.Y();
   g_hits_tgt.momz[g_hits_tgt.nhits] = momini.Z();
   g_hits_tgt.tofg[g_hits_tgt.nhits] = 0.0;
   g_hits_tgt.nhits++;
   g_nhits_tgt_before_chamber++;

   if (arg_writeout_track_fp) {
      fprintf(arg_writeout_track_fp,"%lf %lf %lf %lf %lf %lf 0\n",posini.X(),posini.Y(),posini.Z(),momini.X(),momini.Y(),momini.Z());
      //fprintf(writeout_track_fp,"%lf %lf %lf %lf %lf %lf %lf 0\n",posini.X(),posini.Y(),posini.Z(),momini.X(),momini.Y(),momini.Z(),momini.Mag());
   }
   //printf("momini: %lf %lf %lf\n",momini.X(),momini.Y(),momini.Z());
   //get_init_mom_test(momini);


   stMCT->ResizeTo(rephits->getState());
   *stMCT = rephits->getState();
   covMCT->ResizeTo(rephits->getCov());
   *covMCT = rephits->getCov();

   double tofg;

   double end_of_solenoid = 154.0/2.0; //cm

   time_t t_start = time(NULL);

   const char* pre_volname = "dummy_volname";

   int ilayer;
   int ilayer_prev = -100;
   double step_size;
   step_size = g_target_thickness/2.0;
   while (1) {
      TVector3 pos = rephits->getPos();
      TVector3 mom = rephits->getMom();
      TVector3 dir = mom.Unit();
      
      if (sqrt2(pos.X(),pos.Y())>g_target_radius) {
         step_size = 0.1;
      }

      //printf("step_size %lf\n",step_size);
      TVector3 next_pos = pos + step_size*dir;

      if (arg_writeout_track_fp) {
         fprintf(arg_writeout_track_fp,"%lf %lf %lf %lf %lf %lf",pos.X(),pos.Y(),pos.Z(),mom.X(),mom.Y(),mom.Z());
         //fprintf(writeout_track_fp,"%lf %lf %lf %lf %lf %lf %lf",pos.X(),pos.Y(),pos.Z(),mom.X(),mom.Y(),mom.Z(),mom.Mag());
      }

      // break when particle go out solenoid
      if (abs(pos.Z())>end_of_solenoid) {
         break;
      }
      if (g_hits_det.nhits>=MAX_HIT) {
         break;
      }
      if (g_hits_abs.nhits>=MAX_HIT) {
         break;
      }
      if (g_hits_tgt.nhits>=MAX_HIT) {
         break;
      }
      if (g_hits_biw.nhits>=MAX_HIT) {
         break;
      }
      // wait 3 second
      double e_time = time(NULL) - t_start;
      tv->elapsed_time = e_time;
      if (e_time>=3) {
         break;
      }

      GFDetPlane d(next_pos,dir);

      try{
         tofg = rephits->extrapolate(d);
      }
      catch(GFException& e){
         e.what();
         std::cerr<<"Exceptions wont be further handled ->exit(1)"<<std::endl;
         return -1;
         //exit(1);
      }

      TVector3 posR(rephits->getPos());
      TVector3 momR(rephits->getMom());

      const char* volname = rephits->getVolName();



      int disk_number; // [1..xx] start from 1, same as ini_disk_number
      if (get_hit_at_target(volname,pre_volname,&disk_number)==1) {
         //printf("### iev %d disk_number %d\n",iev,disk_number);
         //printf("momR.mag %lf \n",momR.Mag());
         g_hits_tgt.layer[g_hits_tgt.nhits] = disk_number;
         g_hits_tgt.posx[g_hits_tgt.nhits] = posR.X();
         g_hits_tgt.posy[g_hits_tgt.nhits] = posR.Y();
         g_hits_tgt.posz[g_hits_tgt.nhits] = posR.Z();
         g_hits_tgt.momx[g_hits_tgt.nhits] = momR.X();
         g_hits_tgt.momy[g_hits_tgt.nhits] = momR.Y();
         g_hits_tgt.momz[g_hits_tgt.nhits] = momR.Z();
         g_hits_tgt.tofg[g_hits_tgt.nhits] = tofg;
         g_hits_tgt.nhits++;
         // only count before hit at chamber
         if (g_hits_det.nhits==0) {
            g_nhits_tgt_before_chamber++;
         }
      }


      if (get_hit_at_biw(volname,pre_volname)==1) {
         g_hits_biw.posx[g_hits_biw.nhits] = posR.X();
         g_hits_biw.posy[g_hits_biw.nhits] = posR.Y();
         g_hits_biw.posz[g_hits_biw.nhits] = posR.Z();
         g_hits_biw.momx[g_hits_biw.nhits] = momR.X();
         g_hits_biw.momy[g_hits_biw.nhits] = momR.Y();
         g_hits_biw.momz[g_hits_biw.nhits] = momR.Z();
         g_hits_biw.tofg[g_hits_biw.nhits] = tofg;
         g_hits_biw.nhits++;
      }

      if (get_hit_at_scinti(volname,pre_volname)==1) {
         //printf("momR.mag %lf \n",momR.Mag());
         g_hits_abs.posx[g_hits_abs.nhits] = posR.X();
         g_hits_abs.posy[g_hits_abs.nhits] = posR.Y();
         g_hits_abs.posz[g_hits_abs.nhits] = posR.Z();
         g_hits_abs.momx[g_hits_abs.nhits] = momR.X();
         g_hits_abs.momy[g_hits_abs.nhits] = momR.Y();
         g_hits_abs.momz[g_hits_abs.nhits] = momR.Z();
         g_hits_abs.tofg[g_hits_abs.nhits] = tofg;
         g_hits_abs.nhits++;
      }

      pre_volname = volname;

      int ilayer = get_layer_number(volname);
      //printf("volname %s ilayer %d\n",volname,ilayer);
      // skip if no hit on layer (ilayer==-1) or same hit in the same layer (not to include continues hits)
      if (ilayer==-1 || (ilayer==ilayer_prev)) {
         ilayer_prev = ilayer;
         if (arg_writeout_track_fp) {
            fprintf(arg_writeout_track_fp," 0\n");
         }
         continue;
      }
      if (arg_writeout_track_fp) {
         fprintf(arg_writeout_track_fp," 1\n");
      }
      //fprintf(stderr,"make_true_hit: volname %s\n",volname);
      ilayer_prev = ilayer;

      g_hits_det.layer[g_hits_det.nhits] = ilayer;
      g_hits_det.posx[g_hits_det.nhits] = posR.X();
      g_hits_det.posy[g_hits_det.nhits] = posR.Y();
      g_hits_det.posz[g_hits_det.nhits] = posR.Z();
      g_hits_det.momx[g_hits_det.nhits] = momR.X();
      g_hits_det.momy[g_hits_det.nhits] = momR.Y();
      g_hits_det.momz[g_hits_det.nhits] = momR.Z();
      g_hits_det.tofg[g_hits_det.nhits] = tofg;

      g_hits_det.nhits++;

      //printf("volname %s pre_volname %s nhits (tgt:%d abs:%d det:%d)\n",volname,pre_volname, g_hits_tgt.nhits, g_hits_abs.nhits,g_hits_det.nhits);

   }

   return 0;
}
int calc_pulls(GFAbsTrackRep *rep)
{
   stREC->ResizeTo(rep->getState());
   *stREC = rep->getState();
   covREC->ResizeTo(rep->getCov());
   *covREC = rep->getCov();

   const int iqop=0;
   const int iup=1;
   const int ivp=2;
   const int iu=3;
   const int iv=4;
   const int iqopT=0;
   const int iupT=1;
   const int ivpT=2;
   const int iuT=3;
   const int ivT=4;

   invmom = (*stREC)[iqop][0];
   sigmasqustate = (*covREC)[iqop][iqop];
   //if(sigmasqustate<1.E-16) continue; 
   if(sigmasqustate<1.E-16) return -1;
   sigma_p = 1/pow(invmom,4.) * sigmasqustate;
   momSi=TMath::Sqrt(sigma_p);
   momRe=fabs(1./((*(stREC))[iqop][0]));
   momTr=fabs(1./((*(stMCT))[iqopT][0]));
   momPu=(momRe-momTr)/momSi;

   qopSi=TMath::Sqrt((*covREC)[iqop][iqop]);
   qopRe=(*stREC)[iqop][0];
   qopTr=(*stMCT)[iqopT][0];
   qopPu=(qopRe-qopTr)/qopSi;

   upSi=TMath::Sqrt((*covREC)[iup][iup]);
   upRe=(*stREC)[iup][0];
   upTr=(*stMCT)[iupT][0];
   upPu=(upRe-upTr)/upSi;
   vpSi=TMath::Sqrt((*covREC)[ivp][ivp]);
   vpRe=(*stREC)[ivp][0];
   vpTr=(*stMCT)[ivpT][0];
   vpPu=(vpRe-vpTr)/vpSi;

   uSi=TMath::Sqrt((*covREC)[iu][iu]);
   uRe=(*stREC)[iu][0];
   uTr=(*stMCT)[iuT][0];
   uPu=(uRe-uTr)/uSi;
   vSi=TMath::Sqrt((*covREC)[iv][iv]);
   vRe=(*stREC)[iv][0];
   vTr=(*stMCT)[ivT][0];
   vPu=(vRe-vTr)/vSi;

   return 0;
}

int fitting(int iev, struct tree_value* tv, GFTrack& fitTrack, GFAbsTrackRep* rep, GFDetPlane& plane)
{

   GFKalman k;
   tv->error = 0;
   try {
      k.processTrack(&fitTrack);
   } catch(GFException& e) {
      std::cerr << e.what();
      std::cerr<< "Exceptions wont be further handled ->exit(1)  line " << __LINE__<<std::endl;
      tv->error = 100;
   }

   tv->ndf = rep->getNDF();

   if (rep->getStatusFlag()==0) {
   //fprintf(stdout,"OK! PDGcode %d\n",PDGcode);
      try {
         //fprintf(stderr,"##GFDetPlane: pos.x %lf pos.y %lf pos.z %lf\n",pos.X(),pos.Y(),pos.Z());
         //GFDetPlane plane(pos,TVector3(1.,0.,0.),TVector3(0.,1.,0.));
         
         rep->extrapolate(plane);

      } catch(GFException& e) {
         e.what();
         std::cerr<<"Exceptions wont be further handled ->exit(1)  line " << __LINE__<<std::endl;
         tv->error =2;
      }

			//tv->fit.x = rep->getPos().X();
			//tv->fit.y = rep->getPos().Y();
			//tv->fit.z = rep->getPos().Z();
			tv->fit.x = fitTrack.getPos().X();
			tv->fit.y = fitTrack.getPos().Y();
			tv->fit.z = fitTrack.getPos().Z();

			tv->chi2 = rep->getRedChiSqu();
			//tv->chi2 = fitTrack.getRedChiSqu();
			tv->nfail = fitTrack.getFailedHits();
	 } else {
		 tv->error=3;
	 }
	 /* Do not set chi2 if error occured */
	 if (tv->error!=0) {
		 tv->chi2 = 10E10;
	 }
	 // chi2-probability
	 //tv->prob = TMath::Prob(tv->ndf*tv->chi2,tv->ndf);
	 double chisq = tv->chi2*tv->ndf;
	 tv->prob = TMath::Prob(chisq,tv->ndf);

	 tv->fit.px = fitTrack.getMom().X();
	 tv->fit.py = fitTrack.getMom().Y();
	 tv->fit.pz = fitTrack.getMom().Z();
	 //tv->fit.px = rep->getMom().X();
	 //tv->fit.py = rep->getMom().Y();
	 //tv->fit.pz = rep->getMom().Z();

	 double retval =  calc_pulls(rep);
	 //delete rep; ~GFTrack delete rep;
	 return retval;
}
int get_nhits_for_fitting2()
{
	int nhits_fit;
	int nhits = g_hits_det.nhits;//layer_nums.size();
	std::vector<int> layer_nums;
	for (int i=0; i<nhits; i++) {
		layer_nums.push_back(g_hits_det.layer[i]);
	}

	int count_layer[26];
	int idx_layer[26];
	for (int i=0; i<26; i++) {
		//printf("i %d layer_num %d\n",i,layer_nums[i]);
		count_layer[i] = 0;
		idx_layer[i] = 0;
	}

	int count_first_layer=0;
	for (int i=0; i<nhits; i++) {
		for (int id=1; id<=26; id++) {
			if (layer_nums[i] == id) {
				count_layer[id-1]++;
				if (count_layer[id-1]==2) {
					idx_layer[id-1] = i;
				}
			}
		}
	}
	//for (int i=0; i<26; i++) {
	//   printf("i %d count_layer %d idx_layer %d\n",i,count_layer[i],idx_layer[i]);
	//}
	nhits_fit = nhits;
	for (int i=0; i<26; i++) {
		if (count_layer[i]>=2) { nhits_fit = idx_layer[i]+1;  break; } 
	}

	return nhits_fit;
}
void sscanf_with_unit(char* line, const char* format, double* value_in_cm)
{
	double value;
	char unit[32];
	sscanf(line,format,&value,unit);
	if (strcmp(unit,"cm")==0) { *value_in_cm = value * 1.0; }
	else if (strcmp(unit,"mm")==0) { *value_in_cm = value * 0.1; }
	else if (strcmp(unit,"um")==0) { *value_in_cm = value * 0.0001; }
	else { fprintf(stderr,"unknown unit %s",unit); exit(1); }
}
int read_config(char* fname)
{
	FILE*fp = fopen(fname,"r");
	if (fp==NULL) {
		fprintf(stderr,"read_config: failed to open %s\n",fname);
		return -1;
	}
	char line[128];
	enum { CMD_SOLENOID, CMD_INWALL, CMD_SCINTI, CMD_TARGET, CMD_CHAMBER, CMD_SIZE } cmd;
	cmd = CMD_SIZE;
	while (fgets(line,sizeof(line),fp)) {
		line[strlen(line)-1] = '\0';
		if      (line[0]=='#') continue;
		else if (strcmp(line,"[inwall]")==0) { cmd = CMD_INWALL; }
		else if (strcmp(line,"[scinti]")==0) { cmd = CMD_SCINTI; }
		else if (strcmp(line,"[target]")==0) { cmd = CMD_TARGET; }
		else if (strcmp(line,"[chamber]")==0) { cmd = CMD_CHAMBER; }
		else if (strcmp(line,"[solenoid]")==0) { cmd = CMD_SOLENOID; }

      //printf("cmd -> %d\n",cmd);
      if (cmd==CMD_TARGET) {
         if (strstr(line,"type"))           { sscanf(line,"type %s",g_target_type); }
         else if (strstr(line,"material"))  { sscanf(line,"material %s",g_target_material); }
         else if (strstr(line,"thickness")) { sscanf_with_unit(line,"thickness %lf %s",&g_target_thickness); }
         else if (strstr(line,"offset_z"))  { sscanf_with_unit(line,"offset_z %lf %s",&g_target_offset_z); }
         else if (strstr(line,"radius"))    { sscanf_with_unit(line,"radius %lf %s",&g_target_radius); }
         else if (strstr(line,"spacing"))   { sscanf_with_unit(line,"spacing %lf %s",&g_target_spacing); }
         else if (strstr(line,"numbers"))   { sscanf(line,"numbers %d",&g_target_number_of_disk); }
         else if (strstr(line,"center"))    { sscanf(line,"center %d",&g_target_center_disk_number); }
      } else if (cmd==CMD_SCINTI) {
         if (strstr(line,"type"))           { sscanf(line,"type %s",g_scinti_type); }
         else if (strstr(line,"material"))  { sscanf(line,"material %s",g_scinti_material); }
         else if (strstr(line,"thickness")) { sscanf_with_unit(line,"thickness %lf %s",&g_scinti_thickness); }
         else if (strstr(line,"first_pos")) { sscanf_with_unit(line,"first_pos %lf %s",&g_scinti_first_pos); }
      } else if (cmd==CMD_SOLENOID) {
         if (strstr(line,"type"))           { sscanf(line,"type %s",g_solenoid_bfld_type); }
         else if (strstr(line,"tesla"))     { sscanf(line,"tesla %lf",&g_solenoid_bfld_tesla); }
         else if (strstr(line,"material"))  { sscanf(line,"material %s",g_solenoid_material); }
      } else if (cmd==CMD_INWALL) {
         if (strstr(line,"material"))       { sscanf(line,"material %s",g_inwall_material); }
         else if (strstr(line,"thickness")) { sscanf_with_unit(line,"thickness %lf %s",&g_inwall_thickness); }
         else if (strstr(line,"first_pos")) { sscanf_with_unit(line,"first_pos %lf %s",&g_inwall_first_pos); }
      } else if (cmd==CMD_CHAMBER) {
         if (strstr(line,"type"))           { sscanf(line,"type %s",g_chamber_type); }
         else if (strstr(line,"first_pos")) { sscanf(line,"first_pos %lf",&g_chamber_first_pos); }
         else if (strstr(line,"last_pos"))  { sscanf(line,"last_pos %lf",&g_chamber_last_pos); }
         else if (strstr(line,"spacing"))   { sscanf(line,"spacing %lf",&g_chamber_spacing); }
         else if (strstr(line,"gas"))       { sscanf(line,"gas %s",g_chamber_gas); }
      }
   }
   printf("[inwall]\n");
   printf("first_pos %lf\n",g_inwall_first_pos);
   printf("thickness %lf\n",g_inwall_thickness);
   printf("material %s\n",g_inwall_material);

   printf("[chamber]\n");
   printf("type %s\n",g_chamber_type);
   printf("first_pos %lf\n",g_chamber_first_pos);
   printf("last_pos %lf\n",g_chamber_last_pos);
   printf("spacing %lf\n",g_chamber_spacing);
   g_chamber_num_layer = (g_chamber_last_pos - g_chamber_first_pos)/g_chamber_spacing;
   printf("==> num_layer %d\n",g_chamber_num_layer);
   printf("gas %s\n",g_chamber_gas);

   printf("[scinti]\n");
   printf("type %s\n",g_scinti_type);
   printf("material %s\n",g_scinti_material);
   printf("thickness %lf\n",g_scinti_thickness);
   printf("first_pos %lf\n",g_scinti_first_pos);

   printf("[target]\n");
   printf("type %s\n",g_target_type);
   printf("material %s\n",g_target_material);
   printf("thickness %lf\n",g_target_thickness);
   printf("offset_z %lf\n",g_target_offset_z);
   printf("radius %lf\n",g_target_radius);
   printf("spacing %lf\n",g_target_spacing);
   printf("numbers %d\n",g_target_number_of_disk);
   printf("center %d\n",g_target_center_disk_number);

   printf("[solenoid]\n");
   printf("type %s\n",g_solenoid_bfld_type);
   printf("tesla %lf\n",g_solenoid_bfld_tesla);
   printf("material %s\n",g_solenoid_material);

   fclose(fp);
   return 0;
}
void construct_geom()
{
   double um = 0.0001; // 1 um = 1000 mm = 10000 cm
   double mm = 0.1;

   //--- Definition of a simple geometry
   //   gSystem->Load("libGeom");
   new TGeoManager("genfitGeom", "GENFIT geometry");
   gROOT->Macro("../../geometry/media.C");

   TGeoMedium *vacuum = gGeoManager->GetMedium("vacuum"); assert(vacuum!=NULL);
   TGeoMedium *mylar = gGeoManager->GetMedium("mylar"); assert(mylar!=NULL);
   TGeoMedium *alumi = gGeoManager->GetMedium("alumi"); assert(alumi!=NULL);
   TGeoMedium *scinti = gGeoManager->GetMedium(g_scinti_material); assert(scinti!=NULL);
   TGeoMedium *strgas = gGeoManager->GetMedium(g_chamber_gas); assert(strgas!=NULL);
   TGeoMedium *mat_wall = gGeoManager->GetMedium(g_inwall_material); assert(mat_wall!=NULL);
   TGeoMedium *mat_sol = gGeoManager->GetMedium(g_solenoid_material); assert(mat_sol!=NULL);
   TGeoMedium *mat_target = gGeoManager->GetMedium(g_target_material); assert(mat_sol!=NULL);

   TGeoVolume *top = gGeoManager->MakeBox("TOPPER", vacuum, 500., 500., 500.);
   gGeoManager->SetTopVolume(top); // mandatory !

   // AMY Solenoid
   TGeoVolume *solenoid = gGeoManager->MakeTube("solenoid", mat_sol, 238.6/2.0,258.4/2.0,154.0/2.0);
   solenoid->SetLineColor(kBlack);
   top->AddNode(solenoid, 1, gGeoIdentity);

   //double fFirstPos = g_inwall_first_pos; // 55cm  <= SECOND to reject 70 MeV/c protons
   double fFirstPos = g_chamber_first_pos; // 55cm  <= SECOND to reject 70 MeV/c protons
   //     double fFirstPos = 50; // 50cm  <= to reject 60 MeV/c protons
   double fSpacing = g_chamber_spacing;//1; // 5cm


   double wall_thick = g_inwall_thickness; //400.0*um;
   //double rin_inwall = 54.0;
   //double rin_inwall = fFirstPos - 1;
   double rin_inwall = g_inwall_first_pos;
   double rout_inwall = rin_inwall+wall_thick;
   double rin_outwall = 81.0;
   double rout_outwall = 81.0+wall_thick;
   // Inner wall of drift chamber
   TGeoVolume *inwall = gGeoManager->MakeTube("inwall", mat_wall, rin_inwall, rout_inwall,154.0/2.0);
   //  TGeoVolume *inwall = gGeoManager->MakeTube("inwall", vacuum, rin_inwall, rout_inwall,154.0/2.0);
   inwall->SetLineColor(kRed);
   top->AddNode(inwall, 1, gGeoIdentity);

   // Outer wall of drift chamber
   TGeoVolume *outwall = gGeoManager->MakeTube("outwall", mat_wall, rin_outwall, rout_outwall,154.0/2.0);
   outwall->SetLineColor(kRed);
   top->AddNode(outwall, 1, gGeoIdentity);


   if (strcmp(g_target_type,"disk")==0) {
      int center_disk_number = g_target_center_disk_number;
      int center_disk_idx = center_disk_number - 1;
      int number_of_disk = g_target_number_of_disk;
      if (number_of_disk>100) {
         fprintf(stderr,"too many disks %d\n",number_of_disk);
         exit(1);
      }
      TGeoVolume *disks[100];
      for (int i=0; i<number_of_disk; i++) {
         //disks[i] = gGeoManager->MakeTube(Form("disk%d",i), mat_target, 0.,10.0, g_target_thickness/2.0);
         disks[i] = gGeoManager->MakeTube(Form("disk%d",i), mat_target, 0.,g_target_radius, g_target_thickness/2.0);
         top->AddNode(disks[i], 1, new TGeoTranslation(0,0,g_target_offset_z + g_target_spacing*(i-center_disk_idx)));
      }
   } else if (strcmp(g_target_type,"cylinder")==0) {

      double thickness_cylinder = g_target_thickness; //40.0 *um; // 40um
      double rmin_target = 5.0; //cm
      double rmax_target = rmin_target + thickness_cylinder;
      //double hz_target = 80.0; // cm
      double hz_target = 80.0/2.0; // cm
      TGeoVolume *target_cylinder = gGeoManager->MakeTube("cylinder", mat_target, rmin_target,rmax_target,hz_target);
      top->AddNode(target_cylinder, 1, gGeoIdentity);

   } else if (strcmp(g_target_type,"cone")==0) {
      TGeoVolume *cones[17];
      //TGeoCone(Double_t dz,Double_t rmin1,Double_t rmax1,Double_t rmin2, Double_t rmax2);
      double dz = 5.0;
      double rmin1 = 0.0;
      double rmax1 = g_target_thickness/1.4142135623;
      double rmin2 = g_target_radius;
      double rmax2 = rmin2+rmax1;
      for (int i=0; i<10; i++) {
         if (i%2==0) {
            cones[i] = gGeoManager->MakeCone(Form("cone%d",i), mat_target, dz, rmin1, rmax1, rmin2, rmax2);
         } else {
            cones[i] = gGeoManager->MakeCone(Form("cone%d",i), mat_target, dz, rmin2, rmax2, rmin1, rmax1);
         }
         top->AddNode(cones[i], 1, new TGeoTranslation(0,0,g_target_offset_z + 10.0*(i-4) - 5));
         cones[i]->SetLineColor(kGreen);
      }
   } else {
      fprintf(stderr,"unkonw %s\n",g_target_type);
      exit(1);
   }



   /* Chamber */
   if (strcmp(g_chamber_type,"drift")==0) {

      // Drift wire layers
      TGeoVolume *layers6[26];
      TGeoVolume *layers7[26];
      TGeoVolume *layers8[26];

      // Thickness
      double t_Vt = 10.0* mm; // thickness of wire layer gas region
      double t_Vt1 = 5.0*mm; // thickness of center layer
      double t_Vt2 = 2.5*mm; // thickness of virtual layer around center wire layer

      double din_before = rout_inwall;
      double dout_before = fFirstPos - t_Vt/2.0;
      TGeoVolume* layers_before =  gGeoManager->MakeTube("layer_before", strgas, din_before, dout_before, 154.0/2.0);

      // Before InnerWall Layer
      double biw_rpos = rin_inwall -1.0;
      double biw_thick = 1.0; //
      double biw_rmin = biw_rpos - biw_thick/2.0;
      double biw_rmax = biw_rmin + biw_thick/2.0;
      TGeoVolume* biw_layer = gGeoManager->MakeTube("biw_layer", vacuum, biw_rmin, biw_rmax, 154.0/2.0);
      biw_layer->SetLineColor(kGreen);
      top->AddNode(biw_layer, 1, gGeoIdentity);


      for (int i=0; i<g_chamber_num_layer; i++) {
         double r = fFirstPos+ i*fSpacing;
         // == 
         // == 
         // == t_Vt2 (2.5mm) 
         // == t_Vt1 (5mm)
         // == t_Vt2 (2.5mm)
         double din7  = r-t_Vt1/2.0; // 5mm
         double dout7 = r+t_Vt1/2.0;

         double din6  = din7 - t_Vt2;
         double dout6 = din7;

         double din8  = dout7;
         double dout8 = dout7 + t_Vt2;

         // make virutal layers with thin thickness so that hit position R is near r, such as 55, 60, ...
         layers6[i] =  gGeoManager->MakeTube(Form("layer%d_6",i), strgas, din6, dout6, 154.0/2.0);
         layers7[i] =  gGeoManager->MakeTube(Form("layer%d_7",i), strgas, din7, dout8, 154.0/2.0);
         layers8[i] =  gGeoManager->MakeTube(Form("layer%d_8",i), strgas, din8, dout8, 154.0/2.0);

         if (layers7[i]) layers7[i]->SetLineColor(kGreen);
         if (layers7[i]) top->AddNode(layers7[i], 1, gGeoIdentity);
         //printf("i %d r %lf din6 %lf dout6 %lf\n",i,r,din6,dout6);
         //printf("i %d r %lf din7 %lf dout7 %lf\n",i,r,din7,dout7);
         //printf("i %d r %lf din8 %lf dout8 %lf\n",i,r,din8,dout8);
      }

   } else if (strcmp(g_chamber_type,"layer")==0) {

      // Straw tube layers
      TGeoVolume *layers1[26];
      TGeoVolume *layers2[26];
      TGeoVolume *layers3[26];
      TGeoVolume *layers4[26];
      TGeoVolume *layers5[26];
      TGeoVolume *layers6[26];
      TGeoVolume *layers7[26];
      TGeoVolume *layers8[26];
      for (int i=0; i<g_chamber_num_layer; i++) {
         double r = fFirstPos+ i*fSpacing;

         // Thickness
         double t_Al = 0.1 *um;
         //double t_Al = 1.1 *um;
         //double t_My = 25.0 *um;
         double t_My = 50.0 *um; // to consider double layers
         //double t_Vt = 0.1 *um; <= miss!!
         double t_Vt = 5.0* mm; // gas region


         // din1 == (r+t_Vt/2.0+t_Al+t_My) 
         // din2 == (r+t_Vt/2.0+t_Al)
         // din3 == (r+t_Vt/2.0)
         // din  == (r-t_Vt/2.0)
         // din4 == (r-t_Vt/2.0-t_Al)
         // din5 == (r-t_Vt/2.0-t_Al-t_My)
         // din6 == (r-t_Vt/2.0)
         double din1 = r + t_Vt/2.0+t_Al+t_My;
         double din2 = r + t_Vt/2.0+t_Al;
         double din3 = r + t_Vt/2.0;
         double din  = r - t_Vt/2.0;
         double din4 = r - t_Vt/2.0-t_Al;
         double din5 = r - t_Vt/2.0-t_Al-t_My;
         double din6 = r - t_Vt/2.0-t_Al-t_My-t_Al;

         layers1[i] =  gGeoManager->MakeTube(Form("layer%d_1",i), alumi, din1, din1+t_Al, 154.0/2.0);
         layers2[i] =  gGeoManager->MakeTube(Form("layer%d_2",i), mylar, din2, din2+t_My, 154.0/2.0);
         layers3[i] =  gGeoManager->MakeTube(Form("layer%d_3",i), alumi, din3, din3+t_Al, 154.0/2.0);
         layers4[i] =  gGeoManager->MakeTube(Form("layer%d_4",i), alumi, din4, din4+t_Al, 154.0/2.0);
         layers5[i] =  gGeoManager->MakeTube(Form("layer%d_5",i), mylar, din5, din5+t_My, 154.0/2.0);
         layers6[i] =  gGeoManager->MakeTube(Form("layer%d_6",i), alumi, din6, din6+t_Al, 154.0/2.0);
         layers7[i] =  gGeoManager->MakeTube(Form("layer%d_7",i), strgas, din, din+t_Vt,  154.0/2.0);

         if (layers1[i]) layers1[i]->SetLineColor(kGreen);
         if (layers2[i]) layers2[i]->SetLineColor(kGreen);
         if (layers3[i]) layers3[i]->SetLineColor(kGreen);
         if (layers4[i]) layers4[i]->SetLineColor(kGreen);
         if (layers5[i]) layers5[i]->SetLineColor(kGreen);
         if (layers6[i]) layers6[i]->SetLineColor(kGreen);
         if (layers7[i]) layers7[i]->SetLineColor(kGreen);

         if (layers1[i]) top->AddNode(layers1[i], 1, gGeoIdentity);
         if (layers2[i]) top->AddNode(layers2[i], 1, gGeoIdentity);
         if (layers3[i]) top->AddNode(layers3[i], 1, gGeoIdentity);
         if (layers4[i]) top->AddNode(layers4[i], 1, gGeoIdentity);
         if (layers5[i]) top->AddNode(layers5[i], 1, gGeoIdentity);
         if (layers6[i]) top->AddNode(layers6[i], 1, gGeoIdentity);
         if (layers7[i]) top->AddNode(layers7[i], 1, gGeoIdentity);
      }





   } else {
      fprintf(stderr,"unknown chamber type %s\n",g_chamber_type);
      exit(1);
   }

   /* Scintillating layer */
   if (strcmp(g_scinti_type,"cylinder")==0) {
      //   double scinti_rpos = 50.0; // cm
      //double scinti_rpos = fFirstPos-5; // cm
      double scinti_rpos = g_scinti_first_pos;
      double scinti_thick = g_scinti_thickness; //0.5; // 5mm
      double scinti_rmin = scinti_rpos - scinti_thick/2.0;
      double scinti_rmax = scinti_rmin + scinti_thick/2.0;
      TGeoVolume* scinti_layer = gGeoManager->MakeTube("scinti_layer", scinti, scinti_rmin, scinti_rmax, 154.0/2.0);
      scinti_layer->SetLineColor(kBlue);
      top->AddNode(scinti_layer, 1, gGeoIdentity);
   } else if (strcmp(g_scinti_type,"end-cylinder")==0) {
      fprintf(stderr,"end-cylinder is not implemented\n");
      exit(1);
   } else {
      fprintf(stderr,"unknown scinti type\n");
      exit(1);
   }

   //--- close the geometry
   gGeoManager->CloseGeometry();

   //--- draw the ROOT box
   //gGeoManager->SetVisLevel(10);
   //top->Draw("ogl");

   if (strcmp(arg_genfitGeom_name,"empty")!=0) {
      TFile *outfile = TFile::Open(arg_genfitGeom_name,"RECREATE");
      gGeoManager->Write();
      outfile->Close();
   }
}
void init_args()
{
   arg_seed = -1;
   arg_total = 1; // one event
   arg_writeout_track_fp = NULL;
   strcpy(arg_genfitGeom_name, "empty");
   strcpy(arg_output_root, "genfit.root");
   strcpy(arg_input_root, "empty");
   strcpy(arg_event_type,"SIG"); // SIG,DIO,Proton,TXT
   strcpy(arg_hit_type,"GeaneTrackRep2"); // RKTrackRep,GeaneTrackRep2,ROOT
   strcpy(arg_track_type,"GeaneTrackRep2"); // RKTrackRep,GeaneTrackRep2
   strcpy(arg_extrap_pos,"FirstHit"); //FirtHit,BeforeInwall
}

void print_usage(char* prog_name)
{
   fprintf(stderr,"Usage %s [inputs] [outputs]\n",prog_name);
   fprintf(stderr,"[inputs]\n");
   fprintf(stderr,"\t -s <seed>\n");
   fprintf(stderr,"\t\t Note that seed=0 is an option to set random seed\n");
   fprintf(stderr,"\t -n <number_of_events>\n");
   fprintf(stderr,"\t -c <config.txt>\n");
   fprintf(stderr,"\t\t Check sample config.txt for the format\n");
   fprintf(stderr,"\t -e <event_type>\n");
   fprintf(stderr,"\t\t SIG(104.973MeV/c), DIO(watanabe_shanker, p > <dio_startpoint> (MeV/c), Proton(Silicon spectrum,thre=50MeV/c), TXT(-f <input_par.txt>)\n");
   fprintf(stderr,"\t -m <dio_startpoint>\n");
   fprintf(stderr,"\t\t Unit is MeV/c, 0.101 etc.\n");
   fprintf(stderr,"\t -j <hit_type>\n");
   fprintf(stderr,"\t\t RKTrackRep, GeaneTrackRep2, ROOT(-i <input.root>)\n");
   fprintf(stderr,"\t -t <track_type>\n");
   fprintf(stderr,"\t\t RKTrackRep, GeaneTrackRep2\n");
   fprintf(stderr,"\t -f <ini_param.txt>\n");
   fprintf(stderr,"\t\t format: x y z px py pz\n");
   fprintf(stderr,"\t\t x y z   : initial position in cm\n");
   fprintf(stderr,"\t\t px py pz: initial momentum in GeV/c\n");
   fprintf(stderr,"[outputs]\n");
   fprintf(stderr,"\t -o <output.root>\n");
   fprintf(stderr,"\t\t Check source code for the description on braches\n");
   fprintf(stderr,"\t -g <TGeometryFile.root>\n");
   fprintf(stderr,"\t\t TGeometryFile.root can be viewed by ROOT as follows.\n");
   fprintf(stderr,"\t\t $ root TGeometryFile.root\n");
   fprintf(stderr,"\t\t root[0] genfitGeom->GetTopVolume()->Draw(\"ogl\")\n");
   fprintf(stderr,"\t -z <hit_pos.txt>\n");
   fprintf(stderr,"\t\t format: x y z px py pz hit\n");
   fprintf(stderr,"\t\t x y z   : track position in cm\n");
   fprintf(stderr,"\t\t px py pz: track momentum in GeV/c\n");
   fprintf(stderr,"\t\t hit     : hit flag (1=hit at detector, 0=no hit at detector)\n");
	 fprintf(stderr,"\n");
	 fprintf(stderr,"[example]\n");
	 fprintf(stderr,"\t(1) Run 1k signal events with config8.txt and records output to root/output.root\n");
	 fprintf(stderr,"\t\t%s -s 1 -n 1000 -c user_config/config8.txt -e SIG -j GeaneTrackRep2 -t GeaneTrackRep2 -o root/output.root\n",prog_name);
	 fprintf(stderr,"\t(2) Record TGeo root file for config8.txt\n");
	 fprintf(stderr,"\t\t%s -s 1 -n 1 -c user_config/config8.txt -g tmp/tgeo.root\n",prog_name);
	 fprintf(stderr,"\t(3) Record hit positions in tmp/hit_pos.txt run with the initial parameters in tmp/ini_param.txt\n");
	 fprintf(stderr,"\t\t%s -s 1 -n 1 -c user_config/config8.txt -f tmp/ini_param.txt -z tmp/hit_pos.txt\n",prog_name);
}

int main(int argc, char** argv)
{
	if (argc==1) {
		print_usage(argv[0]);
		return -1;
	}
   init_args();

   int result;
   while((result=getopt(argc,argv,"s:i:n:e:j:t:p:f:c:o:g:z:m:h"))!=-1){
      switch(result){
         /* INPUTS */
         case 's':
            arg_seed = atoi(optarg);
            gRandom->SetSeed(arg_seed);
            printf("gRandom->SetSeed(%d)\n",arg_seed);
            break;
         case 'i':
            strcpy(arg_input_root,optarg);
            printf("input_root %s\n",arg_input_root);
            break;
         case 'n':
            arg_total = atoi(optarg);
            printf("total %d\n",arg_total);
            break;
         case 'e':
            strcpy(arg_event_type,optarg);
            printf("event_type %s\n",arg_event_type);
            break;
         case 'j':
            strcpy(arg_hit_type,optarg);
            printf("hit_type %s\n",arg_hit_type);
            break;
         case 't':
            strcpy(arg_track_type,optarg);
            printf("track_type %s\n",arg_track_type);
            break;
         case 'p':
            strcpy(arg_extrap_pos,optarg);
            printf("extrap_pos %s\n",arg_extrap_pos);
            break;
         case 'f':
            read_init_param(optarg);
            break;
         case 'c':
            if (read_config(optarg)==-1) exit(1);
            break;
         case 'm':
            g_dio_startpoint = atof(optarg);
            printf("dio_startpoint %lf (MeV/c)\n",g_dio_startpoint);
            break;

            /* OUTPUTS */
         case 'o':
            strcpy(arg_output_root,optarg);
            printf("output_root %s\n",arg_output_root);
            break;
         case 'g':
            strcpy(arg_genfitGeom_name,optarg);
            printf("genfitGeom_name %s\n",arg_genfitGeom_name);
            break;
         case 'z':
            arg_writeout_track_fp = fopen(optarg,"w");
            printf("writeout_track: fname %s\n",optarg);
            break;
         case 'h':
         default:
            print_usage(argv[0]);
            exit(1);
      }
   }

   stMCT  = new TMatrixT<double>;
   covMCT = new TMatrixT<double>;
   stREC  = new TMatrixT<double>;
   covREC = new TMatrixT<double>;

   /* Output Root file */
   TFile* fout = new TFile(arg_output_root,"recreate");
   struct tree_value tv;
   TTree* tout = make_branch("t",&tv);
   tout->SetMarkerStyle(2);

   fill_tree_config("t2");

   construct_geom();
   //TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
   //TGeoManager::Import(genfitGeom_name);

   if (strcmp(g_solenoid_bfld_type,"uniform")==0) {
      GFFieldManager::getInstance()->init(new GFConstField(0.,0.,g_solenoid_bfld_tesla*10.0));// kGaus
   } else {
      fprintf(stderr,"bfld_type %s is not implemented\n",g_solenoid_bfld_type);
      exit(1);
   }

   /* Initialize TGeant3 (need in use of GeaneTrackRep2) */
   gROOT->Macro("config/Geane.C");
   //gROOT->Macro("config/Geane_gaus_scat.C");

   //TVector3 posErr(0.01,0.01,0.2); //cm
   //TVector3 posErr(0.01,0.01,1.5); //cm
   TVector3 posErr(0.01,0.01,0.2); //cm

   //TVector3 posErr(0.01,0.01,0.01); //cm
   //TVector3 posErr(1.01,1.01,1.01); //cm
   //TVector3 posErr(100.1,100.1,100.1); //cm
   TVector3 momErr(0.001,0.001,0.001); //GeV/c = 1MeV/c
   //TVector3 momErr(0.000001,0.000001,0.000001); //GeV/c = 1keV/c <<<===== same as above

   for (int iev=0; iev<arg_total; iev++) {

      tv.iev = iev;
      tv.error = -1;
      tv.nfail = -1;
      tv.ndf = 0;
      tv.prob = 0;
      tv.chi2 = 10e10;
      g_nhits_fit = 0;
      g_hits_biw.nhits = 0;
      g_hits_det.nhits = 0;
      g_hits_abs.nhits = 0;
      g_hits_tgt.nhits = 0;
      g_nhits_tgt_before_chamber = 0.0;
      tv.ini_disk_number = -1;

      TVector3 posini,momini;

      // Determine starting position
      if       (strcmp(g_target_type,"disk")==0)     { tv.ini_disk_number = get_init_pos(posini); } 
      else  if (strcmp(g_target_type,"cylinder")==0) { tv.ini_disk_number = get_init_pos_cylinder(posini); } 
      else  if (strcmp(g_target_type,"cone")==0)     { tv.ini_disk_number = get_init_pos_cone(posini); } 
      else {
         fprintf(stderr,"unkonw target type %s\n",g_target_type);
         exit(1);
      }
      // Determine starting momentum
      if      (strcmp(arg_event_type,"SIG")==0)    { get_init_mom(momini); } 
      else if (strcmp(arg_event_type,"DIO")==0)    { get_init_mom_DIO(momini); } 
      else if (strcmp(arg_event_type,"Proton")==0) { PDGcode = 2212; get_init_mom_Proton(momini); } 
      else if (strcmp(arg_event_type,"TXT")==0)    { 
         posini.SetX(arg_posini_x);
         posini.SetY(arg_posini_y);
         posini.SetZ(arg_posini_z);
         momini.SetX(arg_momini_x);
         momini.SetY(arg_momini_y);
         momini.SetZ(arg_momini_z);
      } else {
         fprintf(stderr,"unkown event type %s\n",arg_event_type);
         exit(1);
      }
      tv.ini.x = posini.X();
      tv.ini.y = posini.Y();
      tv.ini.z = posini.Z();
      tv.ini.px = momini.X();
      tv.ini.py = momini.Y();
      tv.ini.pz = momini.Z();

      // Determine hit type
      if ( strcmp(arg_hit_type,"RKTrackRep")==0 || strcmp(arg_hit_type,"GeaneTrackRep2")==0) {

         GFAbsTrackRep* rephits;
         if (strcmp(arg_hit_type,"RKTrackRep")==0) { 
            rephits = new RKTrackRep( posini, momini, PDGcode);
         } else if (strcmp(arg_hit_type,"GeaneTrackRep2")==0) {
            GFDetPlane plane = GFDetPlane(posini,momini);
            rephits = new GeaneTrackRep2( plane, momini, posErr, momErr, PDGcode);
         }

         make_true_hit(iev,&tv, posini,momini,posErr,momErr,rephits);
         delete rephits;

      } else if (strcmp(arg_hit_type,"ROOT")==0) {

         if (read_hit(iev,arg_input_root, posini,momini)==-1) 
            break; // root_hit returns -1 at the end of events

      } else {
         fprintf(stderr,"unknown hit type %s\n",arg_hit_type);
         exit(1);
      }


      /* Start tracking */

      if (strcmp(arg_track_type,"NoFit")==0) {
         tout->Fill();
         continue;
      }

      g_nhits_fit = get_nhits_for_fitting2();
      printf("iev %d numhits %d numhits_fit %d\n",iev,g_hits_det.nhits,g_nhits_fit);

      /* Fitting will be performed when more than 3 hits in the detector */
      if (g_hits_det.nhits<3) {
         //printf("iev skip fitting %d\n",iev);
         tout->Fill();
         continue;
      }

      // Determine extrapolated position
      if (strcmp(arg_extrap_pos,"FirstHit")==0) {

         tv.hit.x = g_hits_det.posx[0];
         tv.hit.y = g_hits_det.posy[0];
         tv.hit.z = g_hits_det.posz[0];
         tv.hit.px = g_hits_det.momx[0];
         tv.hit.py = g_hits_det.momy[0];
         tv.hit.pz = g_hits_det.momz[0];

      } else if (strcmp(arg_extrap_pos,"BeforeInwall")==0) {

         tv.hit.x = g_hits_biw.posx[0];
         tv.hit.y = g_hits_biw.posy[0];
         tv.hit.z = g_hits_biw.posz[0];
         tv.hit.px = g_hits_biw.momx[0];
         tv.hit.py = g_hits_biw.momy[0];
         tv.hit.pz = g_hits_biw.momz[0];

      } else {
         fprintf(stderr,"unknown extrap_pos %s\n",arg_extrap_pos);
         exit(1);
      }

      TVector3 pos_extrap = TVector3(tv.hit.x,tv.hit.y,tv.hit.z);
      TVector3 mom_extrap = TVector3(tv.hit.px,tv.hit.py,tv.hit.pz);
      GFDetPlane plane_extrap = GFDetPlane(pos_extrap,mom_extrap.Unit());
      //GFDetPlane plane_extrap = GFDetPlane(pos_extrap,TVector3(1,0,0),TVector3(0,1,0));

      // Determine track follower
      GFAbsTrackRep* rep;
      if      (strcmp(arg_track_type,"RKTrackRep") == 0 )     { rep = new RKTrackRep( pos_extrap, mom_extrap, PDGcode); } 
      else if (strcmp(arg_track_type,"GeaneTrackRep2") == 0 ) { rep = new GeaneTrackRep2( plane_extrap, mom_extrap, posErr, momErr, PDGcode); }  
      else {
         fprintf(stderr,"unknown track type %s\n",arg_track_type);
         exit(1);
      }

      // Set hit positions
      GFTrack fitTrack(rep);
      int nhits = g_nhits_fit;
      for (int i=0; i<nhits; i++) {
         TVector3 point( g_hits_det.posx[i], g_hits_det.posy[i], g_hits_det.posz[i]);
         fitTrack.addHit(
               new PointHit(point,posErr.X(),posErr.Z()),
               3,//dummy detector id
               i, // hit Id
               i, // rho
               i); // plane Id
      }

      fitting(iev,&tv,fitTrack,rep,plane_extrap);

      print_tree_value("", stdout, &tv);

      tout->Fill();
   }

   fout->cd();
   tout->Write();
   fout->Close();
   if (arg_writeout_track_fp) {
      fclose(arg_writeout_track_fp);
   }

   return 0;
}
