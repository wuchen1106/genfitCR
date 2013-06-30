#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
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

#include "WirepointHit.h"
#include "WireHit.h"
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

#include "TRandom3.h"

//#define MAX_HIT 2000
#define MAX_HIT 20000

double g_rdrift_err; // R err


/* 
 * Proton data
 *    Check how tracking will be effected by rejecting cellls on which protons hit
 * 
 * */
TTree* g_tree_proton;
// input
int g_proton_ntracks; // number of proton tracks
int g_proton_nallhits;
int g_proton_ilayers[MAX_HIT];
int g_proton_icells[MAX_HIT];
int g_proton_iruns[MAX_HIT]; // [g_proton_nallhits]
int g_proton_ievs[MAX_HIT]; // [g_proton_nallhits]
double g_proton_xwires[MAX_HIT]; // [g_proton_nallhits]
double g_proton_ywires[MAX_HIT]; // [g_proton_nallhits]

// check how many hits on wires(sense/field) in the 1st turn (added 2013/6/24)
int g_nhits_sense;
int g_nhits_field;
#define MAX_HIT_ON_WIRE 100
int g_sense_hit_ilayer[MAX_HIT_ON_WIRE];
int g_sense_hit_icell[MAX_HIT_ON_WIRE];
double g_sense_hit_xpos[MAX_HIT_ON_WIRE];
double g_sense_hit_ypos[MAX_HIT_ON_WIRE];
double g_sense_hit_zpos[MAX_HIT_ON_WIRE];
double g_sense_hit_xmom[MAX_HIT_ON_WIRE];
double g_sense_hit_ymom[MAX_HIT_ON_WIRE];
double g_sense_hit_zmom[MAX_HIT_ON_WIRE];


int g_field_hit_ilayer[MAX_HIT_ON_WIRE];
int g_field_hit_icell[MAX_HIT_ON_WIRE];
double g_field_hit_xpos[MAX_HIT_ON_WIRE];
double g_field_hit_ypos[MAX_HIT_ON_WIRE];
double g_field_hit_zpos[MAX_HIT_ON_WIRE];
double g_field_hit_xmom[MAX_HIT_ON_WIRE];
double g_field_hit_ymom[MAX_HIT_ON_WIRE];
double g_field_hit_zmom[MAX_HIT_ON_WIRE];
double g_end_time_of_1st_turn;
int g_last_idx_of_fitting;

// output
int g_proton_cur_iev;
int g_proton_entries;
int g_proton_ov_ncells_uniq; // number of overwapped cells
int g_proton_ov_iruns_uniq[MAX_HIT]; // [g_proton_nov_cells]
int g_proton_ov_ievs_uniq[MAX_HIT]; // [g_proton_nov_cells]
int g_proton_ov_ilayers_uniq[MAX_HIT]; // [g_proton_nov_cells]
int g_proton_ov_icells_uniq[MAX_HIT]; // [g_proton_nov_cells]
double g_proton_ov_xwires_pro_uniq[MAX_HIT];
double g_proton_ov_ywires_pro_uniq[MAX_HIT];
double g_proton_ov_xwires_sig_uniq[MAX_HIT];
double g_proton_ov_ywires_sig_uniq[MAX_HIT];
void count_overwap_cell(int iev, int last_idx);
#define MAX_THE_HIT 1000

/* arguments */
int arg_seed;
char arg_genfitGeom_name[256];
char arg_output_root[128];
int arg_total;
char arg_event_type[128];
char arg_hit_type[128];
char arg_track_type[128];
char arg_extrap_pos[128];
char arg_input_root[128];
char arg_input_txt[128];
char arg_hitpoint_type[128];
FILE* arg_writeout_track_fp;
double arg_posini_x;
double arg_posini_y;
double arg_posini_z;
double arg_momini_x;
double arg_momini_y;
double arg_momini_z;
char arg_proton_data[1024];
char arg_wire_config_fname[1024];
int arg_fitting_ifirst;
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
double g_scinti_length;

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
TH1F* g_hist_RPC;
/* end configure */

/*
 * Configure for DriftChamber (delete after test of Wirepointhit)
 */

//#include "DriftChamber_new6/config.h"
#include "../../geant4_vmc/examples/B01/DriftChamber_new7/config.h"
struct config* g_config;

#include "../../geant4_vmc/examples/C01/calc_ini_track_par/new.h" // IniPar

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


// for debug (2013/05/09)
double g_fitini_theta1;
double g_fitini_theta2;
double g_fitini_dtheta;
double g_fitini_zpos1;
double g_fitini_zpos2;
double g_fitini_dzpos;
// added (2013/5/28)
double g_fitini_1st_wire_posx;
double g_fitini_1st_wire_posy;
double g_fitini_1st_wire_rdrift;

// added 2013/5/21
// Theset are output from IniPar, while tv.fitinit.x are values used for initail fitting
double g_inipar_x; // 1st hit
double g_inipar_y; // 1st hit
double g_inipar_z; // 1st hit
double g_inipar_px; // 1st hit
double g_inipar_py; // 1st hit
double g_inipar_pz; // 1st hit
double g_inipar_rad_1st_hit; // 1st hit
int g_fit_ix; // search range
int g_fit_iy; // search range
int g_fit_iz; // search range
int g_fit_ideg; // search range
int g_fit_min_ipattern;

double g_diff_fitpa; // [for tail cut] difference of fitted pa between that done by ifirst=0 and ifrist=15

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
struct hits g_hits_dum;
int g_nhits_tgt_before_chamber;



// wire infomation

#define MAX_WIRE_HIT 5000
// cells after some cuts(proton hit, rdrift, etc..)
int    g_wire_nhits;
int    g_wire_ilayer[MAX_WIRE_HIT];
int    g_wire_icell[MAX_WIRE_HIT];
double g_wire_rdrift[MAX_WIRE_HIT];
double g_wire_rdrift_smeared[MAX_WIRE_HIT]; // added (2013/3/25)
double g_wire_zreco[MAX_WIRE_HIT];
double g_wire_hitx[MAX_WIRE_HIT];
double g_wire_hity[MAX_WIRE_HIT];
double g_wire_track_length;  // track length during track is inside CDC layer
double g_wire_posx[MAX_WIRE_HIT];
double g_wire_posy[MAX_WIRE_HIT];
double g_wire_momx[MAX_WIRE_HIT];
double g_wire_momy[MAX_WIRE_HIT];
double g_wire_momz[MAX_WIRE_HIT];
double g_wire_tofg[MAX_WIRE_HIT]; // added 2013/6/2

// all of the wire hit cell
int    g_allwire_nhits;
int    g_allwire_cut_type[MAX_WIRE_HIT]; // 0: no cut, 1: proton hit, 2: rdrift_min ....
int    g_allwire_ilayer[MAX_WIRE_HIT];
int    g_allwire_icell[MAX_WIRE_HIT];
double g_allwire_rdrift[MAX_WIRE_HIT];
double g_allwire_rdrift_smeared[MAX_WIRE_HIT];
double g_allwire_zreco[MAX_WIRE_HIT];
double g_allwire_hitx[MAX_WIRE_HIT];
double g_allwire_hity[MAX_WIRE_HIT];
double g_allwire_posx[MAX_WIRE_HIT];
double g_allwire_posy[MAX_WIRE_HIT];
double g_allwire_momx[MAX_WIRE_HIT];
double g_allwire_momy[MAX_WIRE_HIT];
double g_allwire_momz[MAX_WIRE_HIT];
double g_allwire_tofg[MAX_WIRE_HIT];

// for initial value
double g_wire_ini_posx;
double g_wire_ini_posy;
double g_wire_ini_posz;
double g_wire_ini_momx;
double g_wire_ini_momy;
double g_wire_ini_momz;

// Check momentum difference added (2013/5/31)
double g_max_diff_pz;
int g_max_diff_pz_iter;
int g_max_diff_pz_ilayer;

double g_max_diff_pt;
int g_max_diff_pt_iter;
int g_max_diff_pt_ilayer;

double g_max_diff_pa;
int g_max_diff_pa_iter;
int g_max_diff_pa_ilayer;

#define MAX_ENDPLATE_HIT 50
int g_hit_endplate_nhits; // added (2013/6/2)
double g_hit_endplate_tofg[MAX_ENDPLATE_HIT]; // added (2013/6/2)


// For uniq detid list (copied from rpos_zpos)
int g_uniq_npoint;
int g_uniq_detid[MAX_HIT];
double g_uniq_tofg[MAX_HIT];


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

// Diff of pa
double diff_pa_max;

// smeard pos
double g_hits_posx_smeared[MAX_HIT];
double g_hits_posy_smeared[MAX_HIT];
double g_hits_posz_smeared[MAX_HIT];

// smeard pos
double g_hits_exit_posx_smeared;
double g_hits_exit_posy_smeared;
double g_hits_exit_posz_smeared;


// first hit
double f_hit_x;
double f_hit_y;
double f_hit_z;
double f_hit_px;
double f_hit_py;
double f_hit_pz;

double f_fit_x;
double f_fit_y;
double f_fit_z;
double f_fit_px;
double f_fit_py;
double f_fit_pz;

// last hit
double l_hit_x;
double l_hit_y;
double l_hit_z;
double l_hit_px;
double l_hit_py;
double l_hit_pz;

double l_fit_x;
double l_fit_y;
double l_fit_z;
double l_fit_px;
double l_fit_py;
double l_fit_pz;

// sum of diff at each hit points and divide by #/points
double g_diffx;
double g_diffy;
double g_max_diffx;
double g_max_diffy;

void get_position_on_circle(double x0, double y0, double r, double rad, double& x, double& y)
{
   x = TMath::Cos(rad) + x0;
   y = TMath::Sin(rad) + y0;
}

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
   void copy(struct trackpar* other)
   {
      x = other->x;
      y = other->y;
      z = other->z;
      px = other->px;
      py = other->py;
      pz = other->pz;
   };
};


// target
double g_edep_in_target; // GeV
int g_max_ilayer;
int g_nhits_stereo;
int g_nhits_axial;
double g_max_dist;
// overflow added (2013/3/25)
int g_wire_overflow;
int g_wire_seg_overflow;
int g_wire_hit_overflow;

int g_chi2_nhits;
double g_chi2_hit[6][1000]; // ipass,ihit
double g_chi2_hit_max[6];
double g_zfit_hit[6][1000]; // ipass,ihit
double g_zfit_hit_max[6];
double g_track_length_hit[6][1000]; // ipass,ihit
double g_track_length_hit_max[6];
double g_pa_hit[6][1000]; // ipass,ihit
double g_pa_hit_rms[6];
double g_pa_hit_min[6];
double g_pa_hit_max[6];
double g_pz_hit[6][1000]; // ipass,ihit
double g_pz_hit_rms[6];


struct tree_value
{
   int iev;
   int error;
   int nfail;
   int ndf;
   int ifirst;
   int ilast;
   double prob;
   double chi2;
   int ini_disk_number;
   double elapsed_time; // second
   struct trackpar ini;  // hit at start point
   struct trackpar hit;  // first hit
   struct trackpar fit;  // extrapolated at first hit
   struct trackpar fitini;  // initial value for fitting 
   void copy(struct tree_value* other)
   {
      iev              = other->iev;
      error            = other->error;
      nfail            = other->nfail;
      ndf              = other->ndf;
      ifirst           = other->ifirst;
      ilast            = other->ilast;
      prob             = other->prob;
      chi2             = other->chi2;
      ini_disk_number  = other->ini_disk_number;
      ini.copy(&other->ini);
      hit.copy(&other->hit);
      fit.copy(&other->fit);
      fitini.copy(&other->fitini);

   };
   void init(int a_iev)
   {
      iev = a_iev;
      error = 1000;
      nfail = 1000;
      ndf = 0;
      ifirst = 0;
      ilast = 0;
      prob = -100.0;
      chi2 = -100.0;
      ini_disk_number = -100;
      elapsed_time = -100.0;
   };
};

void print_tree_value(const char* prefix, FILE*fp, struct tree_value* tv)
{
   fprintf(fp,"=============\n");
   fprintf(fp,"%s\n",prefix);
   fprintf(fp,"=============\n");
   fprintf(fp,"iev %d\n",tv->iev);
   fprintf(fp,"edep_in_target %f (keV)\n",g_edep_in_target*1e6);
   fprintf(fp,"error %d\n",tv->error);
   fprintf(fp,"nfail %d\n",tv->nfail);
   fprintf(fp,"ndf %d\n",tv->ndf);
   fprintf(fp,"ifirst %d\n",tv->ifirst);
   fprintf(fp,"ilast %d\n",tv->ilast);
   fprintf(fp,"max_ilayer %d\n",g_max_ilayer);

   //fprintf(fp,"nhits %d\n",g_nhits);
   fprintf(fp,"nhits %d\n",g_hits_det.nhits);
   fprintf(fp,"nwrirehits %d\n",g_wire_nhits);
   fprintf(fp,"nhits_tgt %d\n",g_hits_tgt.nhits);
   fprintf(fp,"numhits_tgt_before_chamber %d\n",g_nhits_tgt_before_chamber);

   fprintf(fp,"nhits_fit %d\n",g_nhits_fit);
   fprintf(fp,"chi2 %lf\n",tv->chi2);
   fprintf(fp,"prob %lf\n",tv->prob);
   fprintf(fp,"elapsed_time %lf\n",tv->elapsed_time);
   fprintf(fp,"ini_disk_num %d\n",tv->ini_disk_number);
   double ipt = sqrt2(tv->ini.px,tv->ini.py);
   double ipa = sqrt3(tv->ini.px,tv->ini.py,tv->ini.pz);
   double fitipt = sqrt2(tv->fitini.px,tv->fitini.py);
   double fitipa = sqrt3(tv->fitini.px,tv->fitini.py,tv->fitini.pz);
   fprintf(fp,"start (x,y,z)    %10.5f %10.5f %10.5f\n",tv->ini.x,tv->ini.y,tv->ini.z);
   fprintf(fp,"start (px,py,pz) %10.5f %10.5f %10.5f : pt %10.5lf pa %10.5lf\n",tv->ini.px,tv->ini.py,tv->ini.pz, ipt, ipa);
   fprintf(fp,"IniPar: initial (x,y,z)    %10.5f %10.5f %10.5f\n",tv->fitini.x,tv->fitini.y,tv->fitini.z);
   fprintf(fp,"IniPar: initial (px,py,pz) %10.5f %10.5f %10.5f : pt %10.5lf pa %10.5lf\n",tv->fitini.px,tv->fitini.py,tv->fitini.pz, fitipt, fitipa);
   fprintf(fp,"x(cm)    :  %10.5f => %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.x,  tv->fitini.x, tv->fit.x  , tv->hit.x - tv->fit.x);
   fprintf(fp,"y(cm)    :  %10.5f => %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.y,  tv->fitini.y, tv->fit.y  , tv->hit.y - tv->fit.y);
   fprintf(fp,"z(cm)    :  %10.5f => %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.z,  tv->fitini.z, tv->fit.z  , tv->hit.z - tv->fit.z);
   fprintf(fp,"px(MeV/c):  %10.5f => %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.px*1000.0, tv->fitini.px*1000.0, tv->fit.px*1000.0,  (tv->hit.px- tv->fit.px)*1000.0);
   fprintf(fp,"py(MeV/c):  %10.5f => %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.py*1000.0, tv->fitini.py*1000.0, tv->fit.py*1000.0,  (tv->hit.py- tv->fit.py)*1000.0);
   fprintf(fp,"pz(MeV/c):  %10.5f => %10.5f => %10.5f : diff(%10.5f)\n",tv->hit.pz*1000.0, tv->fitini.pz*1000.0, tv->fit.pz*1000.0,  (tv->hit.pz- tv->fit.pz)*1000.0);
   double hpt = sqrt2(tv->hit.px,tv->hit.py);
   double hpa = sqrt3(tv->hit.px,tv->hit.py,tv->hit.pz);
   double fpt = sqrt2(tv->fit.px,tv->fit.py);
   double fpa = sqrt3(tv->fit.px,tv->fit.py,tv->fit.pz);
   fprintf(fp,"pt(MeV/c):  %10.5f => %10.5f : diff(%10.5f)\n",hpt*1000.0, fpt*1000.0, ( hpt- fpt)*1000.0);
   fprintf(fp,"pa(MeV/c):  %10.5f => %10.5f : diff(%10.5f)\n",hpa*1000.0, fpa*1000.0, ( hpa- fpa)*1000.0);
   fprintf(fp,"diff_pa(MeV/c) %10.5f\n", g_diff_fitpa*1000.0);
   fprintf(fp,"=============\n");
   if (abs(hpa-fpa)*1000.0 >2 ) {
      printf("hoge iev %d abs(hpa-fpa)>2 MeV/c\n",tv->iev);
   }
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
   t->Branch("hitpoint_type",arg_hitpoint_type,"hitpoint_type/C");
   t->Branch("fitting_ifirst",&arg_fitting_ifirst,"fitting_ifirst/I");

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
   t->Branch("scinti_length",&g_scinti_length,"g_scinti_length/D");

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
   t->Branch("ifirst",&tv->ifirst,"ifirst/I");
   t->Branch("ilast",&tv->ilast,"ilast/I");
   t->Branch("chi2",&tv->chi2,"chi2/D");
   t->Branch("prob",&tv->prob,"prob/D"); // chi2-probability


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

   // wire infos
   t->Branch("wire_nhits",&g_wire_nhits,"wire_nhits/I");
   t->Branch("wire_ilayer",g_wire_ilayer,"wire_ilayer[wire_nhits]/I");
   t->Branch("wire_icell" ,g_wire_icell ,"wire_icell[wire_nhits]/I");
   t->Branch("wire_rdrift",g_wire_rdrift ,"wire_rdrift[wire_nhits]/D");
   t->Branch("wire_rdrift_smeared",g_wire_rdrift_smeared ,"wire_rdrift_smeared[wire_nhits]/D"); // added 2013/3/25
   t->Branch("wire_zreco", g_wire_zreco ,"wire_zreco[wire_nhits]/D");
   //    additional info
   t->Branch("wire_hitx", g_wire_hitx ,"wire_hitx[wire_nhits]/D");
   t->Branch("wire_hity", g_wire_hity ,"wire_hity[wire_nhits]/D");
   t->Branch("wire_track_length", &g_wire_track_length ,"wire_track_length/D");
   
   t->Branch("wire_posx", g_wire_posx ,"wire_posx[wire_nhits]/D");
   t->Branch("wire_posy", g_wire_posy ,"wire_posy[wire_nhits]/D");

   // added (2013/3/15)
   t->Branch("wire_momx", g_wire_momx ,"wire_momx[wire_nhits]/D");
   t->Branch("wire_momy", g_wire_momy ,"wire_momy[wire_nhits]/D");
   t->Branch("wire_momz", g_wire_momz ,"wire_momz[wire_nhits]/D");
   t->Branch("wire_tofg", g_wire_tofg ,"wire_tofg[wire_nhits]/D");

   // added (2013/5/31)
   t->Branch("max_diff_pz", &g_max_diff_pz ,"max_diff_pz/D");
   t->Branch("max_diff_pz_iter", &g_max_diff_pz_iter ,"max_diff_pz_iter/I");
   t->Branch("max_diff_pz_ilayer", &g_max_diff_pz_ilayer ,"max_diff_pz_ilayer/I");
   t->Branch("max_diff_pt", &g_max_diff_pt ,"max_diff_pt/D");
   t->Branch("max_diff_pt_iter", &g_max_diff_pt_iter ,"max_diff_pt_iter/I");
   t->Branch("max_diff_pt_ilayer", &g_max_diff_pt_ilayer ,"max_diff_pt_ilayer/I");
   t->Branch("max_diff_pa", &g_max_diff_pa ,"max_diff_pa/D");
   t->Branch("max_diff_pa_iter", &g_max_diff_pa_iter ,"max_diff_pa_iter/I");
   t->Branch("max_diff_pa_ilayer", &g_max_diff_pa_ilayer ,"max_diff_pa_ilayer/I");

   // added (2013/6/2)
   t->Branch("hit_endplate_nhits", &g_hit_endplate_nhits,"hit_endplate_nhits/I");
   t->Branch("hit_endplate_tofg", g_hit_endplate_tofg ,"hit_endplate_tofg[hit_endplate_nhits]/D");
   t->Branch("uniq_npoint", &g_uniq_npoint,"uniq_npoint/I");
   t->Branch("uniq_detid", g_uniq_detid ,"uniq_detid[uniq_npoint]/I");
   t->Branch("uniq_tofg", g_uniq_tofg ,"uniq_tofg[uniq_npoint]/D");

   t->Branch("proton_cur_iev", &g_proton_cur_iev ,"proton_cur_iev/I");
   t->Branch("proton_ntracks", &g_proton_ntracks ,"proton_ntracks/I");
   t->Branch("proton_ov_ncells", &g_proton_ov_ncells_uniq ,"proton_ov_ncells/I");
   t->Branch("proton_ov_iruns", g_proton_ov_iruns_uniq ,"proton_ov_iruns[proton_ov_ncells]/I");
   t->Branch("proton_ov_ievs",  g_proton_ov_ievs_uniq ,"proton_ov_ievs[proton_ov_ncells]/I");
   t->Branch("proton_ov_ilayers",  g_proton_ov_ilayers_uniq ,"proton_ov_ilayers[proton_ov_ncells]/I");
   t->Branch("proton_ov_icells",  g_proton_ov_icells_uniq ,"proton_ov_icells[proton_ov_ncells]/I");
   t->Branch("proton_ov_xwires_pro",  g_proton_ov_xwires_pro_uniq ,"proton_ov_xwires_pro[proton_ov_ncells]/D"); // z is for proton
   t->Branch("proton_ov_ywires_pro",  g_proton_ov_ywires_pro_uniq ,"proton_ov_ywires_pro[proton_ov_ncells]/D"); // z is for proton
   t->Branch("proton_ov_xwires_sig",  g_proton_ov_xwires_sig_uniq ,"proton_ov_xwires_sig[proton_ov_ncells]/D"); // z is for signal
   t->Branch("proton_ov_ywires_sig",  g_proton_ov_ywires_sig_uniq ,"proton_ov_ywires_sig[proton_ov_ncells]/D"); // z is for signal


   // all wire infos (added 2013/6/4)
   t->Branch("allwire_nhits",        &g_allwire_nhits,           "allwire_nhits/I");
   t->Branch("allwire_cut_type",      g_allwire_cut_type,        "allwire_cut_type[allwire_nhits]/I");
   t->Branch("allwire_ilayer",        g_allwire_ilayer,          "allwire_ilayer[allwire_nhits]/I");
   t->Branch("allwire_icell" ,        g_allwire_icell ,          "allwire_icell[allwire_nhits]/I");
   t->Branch("allwire_rdrift",        g_allwire_rdrift ,         "allwire_rdrift[allwire_nhits]/D");
   t->Branch("allwire_rdrift_smeared",g_allwire_rdrift_smeared , "allwire_rdrift_smeared[allwire_nhits]/D");
   t->Branch("allwire_zreco",         g_allwire_zreco ,          "allwire_zreco[allwire_nhits]/D");
   t->Branch("allwire_hitx",          g_allwire_hitx ,           "allwire_hitx[allwire_nhits]/D");
   t->Branch("allwire_hity",          g_allwire_hity ,           "allwire_hity[allwire_nhits]/D");
   t->Branch("allwire_posx",          g_allwire_posx ,           "allwire_posx[allwire_nhits]/D");
   t->Branch("allwire_posy",          g_allwire_posy ,           "allwire_posy[allwire_nhits]/D");
   t->Branch("allwire_momx",          g_allwire_momx ,           "allwire_momx[allwire_nhits]/D");
   t->Branch("allwire_momy",          g_allwire_momy ,           "allwire_momy[allwire_nhits]/D");
   t->Branch("allwire_momz",          g_allwire_momz ,           "allwire_momz[allwire_nhits]/D");
   t->Branch("allwire_tofg",          g_allwire_tofg ,           "allwire_tofg[allwire_nhits]/D");



   // Number of hits on wires (2013/6/24)
   t->Branch("wirehit_sense_nhits", &g_nhits_sense ,   "wirehit_sense_nhits/I");
   t->Branch("wirehit_sense_ilayer",g_sense_hit_ilayer,  "wirehit_sense_ilayer[wirehit_sense_nhits]/I");
   t->Branch("wirehit_sense_icell", g_sense_hit_icell ,  "wirehit_sense_icell[wirehit_sense_nhits]/I");
   t->Branch("wirehit_sense_xpos", g_sense_hit_xpos ,  "wirehit_sense_xpos[wirehit_sense_nhits]/D");
   t->Branch("wirehit_sense_ypos", g_sense_hit_ypos ,  "wirehit_sense_ypos[wirehit_sense_nhits]/D");
   t->Branch("wirehit_sense_zpos", g_sense_hit_zpos ,  "wirehit_sense_zpos[wirehit_sense_nhits]/D");
   t->Branch("wirehit_sense_xmom", g_sense_hit_xmom ,  "wirehit_sense_xmom[wirehit_sense_nhits]/D");
   t->Branch("wirehit_sense_ymom", g_sense_hit_ymom ,  "wirehit_sense_ymom[wirehit_sense_nhits]/D");
   t->Branch("wirehit_sense_zmom", g_sense_hit_zmom ,  "wirehit_sense_zmom[wirehit_sense_nhits]/D");

   t->Branch("wirehit_field_nhits", &g_nhits_field ,   "wirehit_field_nhits/I");
   t->Branch("wirehit_field_ilayer",g_field_hit_ilayer,  "wirehit_field_ilayer[wirehit_field_nhits]/I");
   t->Branch("wirehit_field_icell", g_field_hit_icell ,  "wirehit_field_icell[wirehit_field_nhits]/I");
   t->Branch("wirehit_field_xpos", g_field_hit_xpos ,  "wirehit_field_xpos[wirehit_field_nhits]/D");
   t->Branch("wirehit_field_ypos", g_field_hit_ypos ,  "wirehit_field_ypos[wirehit_field_nhits]/D");
   t->Branch("wirehit_field_zpos", g_field_hit_zpos ,  "wirehit_field_zpos[wirehit_field_nhits]/D");
   t->Branch("wirehit_field_xmom", g_field_hit_xmom ,  "wirehit_field_xmom[wirehit_field_nhits]/D");
   t->Branch("wirehit_field_ymom", g_field_hit_ymom ,  "wirehit_field_ymom[wirehit_field_nhits]/D");
   t->Branch("wirehit_field_zmom", g_field_hit_zmom ,  "wirehit_field_zmom[wirehit_field_nhits]/D");

   t->Branch("endtime_of_1st_turn", &g_end_time_of_1st_turn , "endtime_of_1st_turn/D");
   t->Branch("last_idx_of_fitting", &g_last_idx_of_fitting ,  "last_idx_of_fitting/I"); // <= this shoudl be wire_nhits -1

   // Deposit energy in target
   t->Branch("edep_in_target", &g_edep_in_target ,"edep_in_target/D");
   
   
   // check how many stereo layers does track pass through (2013/3/12)
   t->Branch("max_ilayer",  &g_max_ilayer, "max_ilayer/I");
   
   // added chi2 at each hit points (2013/3/23)
   t->Branch("chi2_nhits",  &g_chi2_nhits, "chi2_nhits/I");
   t->Branch("chi2_hit1",   &g_chi2_hit[0], "chi2_hit1[chi2_nhits]/D");
   t->Branch("chi2_hit2",   &g_chi2_hit[1], "chi2_hit2[chi2_nhits]/D");
   t->Branch("chi2_hit3",   &g_chi2_hit[2], "chi2_hit3[chi2_nhits]/D");
   t->Branch("chi2_hit4",   &g_chi2_hit[3], "chi2_hit4[chi2_nhits]/D");
   t->Branch("chi2_hit5",   &g_chi2_hit[4], "chi2_hit5[chi2_nhits]/D");
   t->Branch("chi2_hit6",   &g_chi2_hit[5], "chi2_hit6[chi2_nhits]/D");

   t->Branch("chi2_hit1_max",   &g_chi2_hit_max[0], "chi2_hit1_max/D");
   t->Branch("chi2_hit2_max",   &g_chi2_hit_max[1], "chi2_hit2_max/D");
   t->Branch("chi2_hit3_max",   &g_chi2_hit_max[2], "chi2_hit3_max/D");
   t->Branch("chi2_hit4_max",   &g_chi2_hit_max[3], "chi2_hit4_max/D");
   t->Branch("chi2_hit5_max",   &g_chi2_hit_max[4], "chi2_hit5_max/D");
   t->Branch("chi2_hit6_max",   &g_chi2_hit_max[5], "chi2_hit6_max/D");

   // added fitted_zpos at each hit points (2013/4/15)
   t->Branch("zfit_hit1",   &g_zfit_hit[0], "zfit_hit1[chi2_nhits]/D");
   t->Branch("zfit_hit2",   &g_zfit_hit[1], "zfit_hit2[chi2_nhits]/D");
   t->Branch("zfit_hit3",   &g_zfit_hit[2], "zfit_hit3[chi2_nhits]/D");
   t->Branch("zfit_hit4",   &g_zfit_hit[3], "zfit_hit4[chi2_nhits]/D");
   t->Branch("zfit_hit5",   &g_zfit_hit[4], "zfit_hit5[chi2_nhits]/D");
   t->Branch("zfit_hit6",   &g_zfit_hit[5], "zfit_hit6[chi2_nhits]/D");

   t->Branch("zfit_hit1_max",   &g_zfit_hit_max[0], "zfit_hit1_max/D");
   t->Branch("zfit_hit2_max",   &g_zfit_hit_max[1], "zfit_hit2_max/D");
   t->Branch("zfit_hit3_max",   &g_zfit_hit_max[2], "zfit_hit3_max/D");
   t->Branch("zfit_hit4_max",   &g_zfit_hit_max[3], "zfit_hit4_max/D");
   t->Branch("zfit_hit5_max",   &g_zfit_hit_max[4], "zfit_hit5_max/D");
   t->Branch("zfit_hit6_max",   &g_zfit_hit_max[5], "zfit_hit6_max/D");

   // added track_length at each hit points (2013/4/15)
   t->Branch("track_length_hit1",   &g_track_length_hit[0], "track_length_hit1[chi2_nhits]/D");
   t->Branch("track_length_hit2",   &g_track_length_hit[1], "track_length_hit2[chi2_nhits]/D");
   t->Branch("track_length_hit3",   &g_track_length_hit[2], "track_length_hit3[chi2_nhits]/D");
   t->Branch("track_length_hit4",   &g_track_length_hit[3], "track_length_hit4[chi2_nhits]/D");
   t->Branch("track_length_hit5",   &g_track_length_hit[4], "track_length_hit5[chi2_nhits]/D");
   t->Branch("track_length_hit6",   &g_track_length_hit[5], "track_length_hit6[chi2_nhits]/D");

   t->Branch("track_length_hit1_max",   &g_track_length_hit_max[0], "track_length_hit1_max/D");
   t->Branch("track_length_hit2_max",   &g_track_length_hit_max[1], "track_length_hit2_max/D");
   t->Branch("track_length_hit3_max",   &g_track_length_hit_max[2], "track_length_hit3_max/D");
   t->Branch("track_length_hit4_max",   &g_track_length_hit_max[3], "track_length_hit4_max/D");
   t->Branch("track_length_hit5_max",   &g_track_length_hit_max[4], "track_length_hit5_max/D");
   t->Branch("track_length_hit6_max",   &g_track_length_hit_max[5], "track_length_hit6_max/D");

   // added pa at each hit points (2013/4/16)
   t->Branch("pa_hit1",   &g_pa_hit[0], "pa_hit1[chi2_nhits]/D");
   t->Branch("pa_hit2",   &g_pa_hit[1], "pa_hit2[chi2_nhits]/D");
   t->Branch("pa_hit3",   &g_pa_hit[2], "pa_hit3[chi2_nhits]/D");
   t->Branch("pa_hit4",   &g_pa_hit[3], "pa_hit4[chi2_nhits]/D");
   t->Branch("pa_hit5",   &g_pa_hit[4], "pa_hit5[chi2_nhits]/D");
   t->Branch("pa_hit6",   &g_pa_hit[5], "pa_hit6[chi2_nhits]/D");

   // added pz at each hit points (2013/5/31)
   t->Branch("pz_hit1",   &g_pz_hit[0], "pz_hit1[chi2_nhits]/D");
   t->Branch("pz_hit2",   &g_pz_hit[1], "pz_hit2[chi2_nhits]/D");
   t->Branch("pz_hit3",   &g_pz_hit[2], "pz_hit3[chi2_nhits]/D");
   t->Branch("pz_hit4",   &g_pz_hit[3], "pz_hit4[chi2_nhits]/D");
   t->Branch("pz_hit5",   &g_pz_hit[4], "pz_hit5[chi2_nhits]/D");
   t->Branch("pz_hit6",   &g_pz_hit[5], "pz_hit6[chi2_nhits]/D");

   t->Branch("pz_hit1_rms",   &g_pz_hit_rms[0], "pz_hit1_rms/D");
   t->Branch("pz_hit2_rms",   &g_pz_hit_rms[1], "pz_hit2_rms/D");
   t->Branch("pz_hit3_rms",   &g_pz_hit_rms[2], "pz_hit3_rms/D");
   t->Branch("pz_hit4_rms",   &g_pz_hit_rms[3], "pz_hit4_rms/D");
   t->Branch("pz_hit5_rms",   &g_pz_hit_rms[4], "pz_hit5_rms/D");
   t->Branch("pz_hit6_rms",   &g_pz_hit_rms[5], "pz_hit6_rms/D");

   t->Branch("pa_hit1_rms",   &g_pa_hit_rms[0], "pa_hit1_rms/D");
   t->Branch("pa_hit2_rms",   &g_pa_hit_rms[1], "pa_hit2_rms/D");
   t->Branch("pa_hit3_rms",   &g_pa_hit_rms[2], "pa_hit3_rms/D");
   t->Branch("pa_hit4_rms",   &g_pa_hit_rms[3], "pa_hit4_rms/D");
   t->Branch("pa_hit5_rms",   &g_pa_hit_rms[4], "pa_hit5_rms/D");
   t->Branch("pa_hit6_rms",   &g_pa_hit_rms[5], "pa_hit6_rms/D");

   // rms is not good ? (2013/05/16)
   t->Branch("pa_hit1_min",   &g_pa_hit_min[0], "pa_hit1_min/D");
   t->Branch("pa_hit2_min",   &g_pa_hit_min[1], "pa_hit2_min/D");
   t->Branch("pa_hit3_min",   &g_pa_hit_min[2], "pa_hit3_min/D");
   t->Branch("pa_hit4_min",   &g_pa_hit_min[3], "pa_hit4_min/D");
   t->Branch("pa_hit5_min",   &g_pa_hit_min[4], "pa_hit5_min/D");
   t->Branch("pa_hit6_min",   &g_pa_hit_min[5], "pa_hit6_min/D");

   t->Branch("pa_hit1_max",   &g_pa_hit_max[0], "pa_hit1_max/D");
   t->Branch("pa_hit2_max",   &g_pa_hit_max[1], "pa_hit2_max/D");
   t->Branch("pa_hit3_max",   &g_pa_hit_max[2], "pa_hit3_max/D");
   t->Branch("pa_hit4_max",   &g_pa_hit_max[3], "pa_hit4_max/D");
   t->Branch("pa_hit5_max",   &g_pa_hit_max[4], "pa_hit5_max/D");
   t->Branch("pa_hit6_max",   &g_pa_hit_max[5], "pa_hit6_max/D");





   // check (summed) number of hits in stereo layers
   t->Branch("nhits_stereo",  &g_nhits_stereo, "nhits_stereo/I");
   t->Branch("nhits_axial",  &g_nhits_axial, "nhits_axial/I");
   
   // maximum track length between neighboring hit points (2013/3/21)
   t->Branch("max_dist",  &g_max_dist, "max_dist/D");

   // overflow (set at B01)
   t->Branch("wire_overflow",      &g_wire_overflow, "wire_overflow/I");
   t->Branch("wire_seg_overflow",  &g_wire_seg_overflow, "wire_seg_overflow/I");
   t->Branch("wire_hit_overflow",  &g_wire_hit_overflow, "wire_hit_overflow/I");



   // first hit
   t->Branch("hif_x",&tv->hit.x,"hit_x/D");
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

   // Added fitting initial vlue (2013/05/09)
   // Initial value for fitting 
   t->Branch("fitini_x",&tv->fitini.x,"fitini_x/D");
   t->Branch("fitini_y",&tv->fitini.y,"fitini_y/D");
   t->Branch("fitini_z",&tv->fitini.z,"fitini_z/D");
   t->Branch("fitini_px",&tv->fitini.px,"fitini_px/D");
   t->Branch("fitini_py",&tv->fitini.py,"fitini_py/D");
   t->Branch("fitini_pz",&tv->fitini.pz,"fitini_pz/D");
   t->SetAlias("fitini_pt","sqrt(fitini_px*fitini_px+fitini_py*fitini_py)");
   t->SetAlias("fitini_pa","sqrt(fitini_pt*fitini_pt+fitini_pz*fitini_pz)");
   t->SetAlias("fitini_r","sqrt(fitini_x*fitini_x+fitini_y*fitini_y)");
   t->Branch("fitini_theta1",&g_fitini_theta1,"fitini_theta1/D");
   t->Branch("fitini_theta2",&g_fitini_theta2,"fitini_theta2/D");
   t->Branch("fitini_dtheta",&g_fitini_dtheta,"fitini_dtheta/D");
   t->Branch("fitini_zpos1",&g_fitini_zpos1,"fitini_zpos1/D");
   t->Branch("fitini_zpos2",&g_fitini_zpos2,"fitini_zpos2/D");
   t->Branch("fitini_dzpos",&g_fitini_dzpos,"fitini_dzpos/D");

   // Added 2013/5/21
   t->Branch("inipar_x", &g_inipar_x, "inipar_x/D");
   t->Branch("inipar_y", &g_inipar_y, "inipar_y/D");
   t->Branch("inipar_z", &g_inipar_z, "inipar_z/D");
   t->Branch("inipar_px",&g_inipar_px,"inipar_px/D");
   t->Branch("inipar_py",&g_inipar_py,"inipar_py/D");
   t->Branch("inipar_pz",&g_inipar_pz,"inipar_pz/D");
   t->Branch("inipar_rad_1st_hit",&g_inipar_rad_1st_hit,"inipar_rad_1st_hit/D");
   t->Branch("fit_ix",&g_fit_ix,"fit_ix/I");
   t->Branch("fit_iy",&g_fit_iy,"fit_iy/I");
   t->Branch("fit_iz",&g_fit_iz,"fit_iz/I");
   t->Branch("fit_ideg",&g_fit_ideg,"fit_ideg/I");
   t->Branch("fit_min_ipattern",&g_fit_min_ipattern,"min_ipattern/I");

   // added 2013/6/7
   t->Branch("diff_fitpa",&g_diff_fitpa,"diff_fitpa/D");


   // comment out 2013/05/09
#if 0

   // last hit
   t->Branch("l_hit_x",&l_hit_x,"l_hit_x/D");
   t->Branch("l_hit_y",&l_hit_y,"l_hit_y/D");
   t->Branch("l_hit_z",&l_hit_z,"l_hit_z/D");
   t->Branch("l_hit_px",&l_hit_px,"l_hit_px/D");
   t->Branch("l_hit_py",&l_hit_py,"l_hit_py/D");
   t->Branch("l_hit_pz",&l_hit_pz,"l_hit_pz/D");
   t->SetAlias("l_hit_pt","sqrt(l_hit_px*l_hit_px+l_hit_py*l_hit_py)");
   t->SetAlias("l_hit_pa","sqrt(l_hit_pt*l_hit_pt+l_hit_pz*l_hit_pz)");
   t->SetAlias("l_hit_r","sqrt(l_hit_x*l_hit_x+l_hit_y*l_hit_y)");

   // extraplated at last hit
   t->Branch("l_fit_x",&l_fit_x,"l_fit_x/D");
   t->Branch("l_fit_y",&l_fit_y,"l_fit_y/D");
   t->Branch("l_fit_z",&l_fit_z,"l_fit_z/D");
   t->Branch("l_fit_px",&l_fit_px,"l_fit_px/D");
   t->Branch("l_fit_py",&l_fit_py,"l_fit_py/D");
   t->Branch("l_fit_pz",&l_fit_pz,"l_fit_pz/D");
   t->SetAlias("l_fit_pt","sqrt(l_fit_px*l_fit_px+l_fit_py*l_fit_py)");
   t->SetAlias("l_fit_pa","sqrt(l_fit_pt*l_fit_pt+l_fit_pz*l_fit_pz)");
   t->SetAlias("l_fit_r","sqrt(l_fit_x*l_fit_x+l_fit_y*l_fit_y)");

   // Residuals at last hit
   t->SetAlias("l_res_x","(l_hit_x-l_fit_x)");
   t->SetAlias("l_res_y","(l_hit_y-l_fit_y)");
   t->SetAlias("l_res_r","(l_hit_r-l_fit_r)");
   t->SetAlias("l_res_px","(l_hit_px-l_fit_px)");
   t->SetAlias("l_res_py","(l_hit_py-l_fit_py)");
   t->SetAlias("l_res_pz","(l_hit_pz-l_fit_pz)");
   t->SetAlias("l_res_pt","(l_hit_pt-l_fit_pt)");
   t->SetAlias("l_res_pa","(l_hit_pa-l_fit_pa)");


   // sum diff per
   t->Branch("diffx",&g_diffx,"diffx/D");
   t->Branch("diffy",&g_diffy,"diffy/D");
   t->Branch("max_diffx",&g_max_diffx,"max_diffx/D");
   t->Branch("max_diffy",&g_max_diffy,"max_diffy/D");
#endif


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
   //t->Branch("hits_abs_layer",g_hits_abs_hits_abs_abs.layer,"hits_abs_layer[nhits_abs]/I");
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
   // these are hits points actually used for tracking
   t->Branch("nhits_fit",&g_nhits_fit,"nhits_fit/I"); // number of (first series of) hits used for fitting_
   t->Branch("hits_x_smeared",g_hits_posx_smeared,"hits_x_smeared[nhits]/D");
   t->Branch("hits_y_smeared",g_hits_posy_smeared,"hits_y_smeared[nhits]/D");
   t->Branch("hits_z_smeared",g_hits_posz_smeared,"hits_z_smeared[nhits]/D");

   t->Branch("hits_x_exit_smeared",&g_hits_exit_posx_smeared,"hits_x_exit_smeared/D");
   t->Branch("hits_y_exit_smeared",&g_hits_exit_posy_smeared,"hits_y_exit_smeared/D");
   t->Branch("hits_z_exit_smeared",&g_hits_exit_posz_smeared,"hits_z_exit_smeared/D");



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


#if 1 //deleted (2013/3/23) => again 
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
#endif

   // Diff of pa
   t->Branch("diff_pa_max",&diff_pa_max,"diff_pa_max/D");
   return t;

}
// DIO spectrum

struct DIO
{
    TRandom3 *rndm;
    double xpos[86];
    double ypos[86];
    DIO() { rndm = new TRandom3; };
    int Read(char* txt)
    {
        // Data points
        FILE* fp = fopen(txt,"r");
        if (fp==0)
            return -1;

        char line[1024];
        int i=0;
        while(fgets(line,sizeof(line),fp)) {
            if (i==86) break;
            sscanf(line,"%lf %lf",&xpos[i],&ypos[i]);
            i++;
        }
        fclose(fp);

        return 0;
    };
    void Print()
    {
        for (int i=0; i<86; i++) {
            printf("%d %lf %lf\n",i,xpos[i],ypos[i]);
        }
    };
    void GetRandom(double& x)
    {
        while(1) {
            int ix = (int)rndm->Uniform(0,85);
            double yval = rndm->Uniform(0,0.04);
            if (ypos[ix]>yval) {
                x = xpos[ix];
                return;
            }
        }
    };
    void Plot(int num)
    {
        TH1F*h1 = new TH1F("h1","DIO",200,0,105);
        double x;
        for (int i=0; i<num; i++) {
            GetRandom(x);
            h1->GetXaxis()->FindBin(x);
            h1->Fill(x);
        }
        h1->Draw();
    };
};
struct DIO *dio_spec=NULL;

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
int get_init_pos_test(TVector3& pos)
{
   pos.SetXYZ(0.,0.,-0.041); //  just before disk target
   return 9; // ceter disk number
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
//   double ymom = 0.0;
//   //double xmom=0.1049;
//   double xmom=0.1043;
//   //double zmom = 0.075; // OK
//   //double zmom = 0.073;
//   double zmom = sqrt(0.105*0.105-xmom*xmom-zmom*zmom);
//   mom.SetXYZ(xmom,ymom,zmom);
   
   //double a = 0.1034;
   //double b = 0.0182;
   //mom.SetXYZ(a,0.,b);
   
   mom.SetXYZ(0.,0.,ptotal);
}
// Spectrum is made from Figure 3 in  Physical Review D 84, 013006 (2011) */
void get_init_mom_DIOA(TVector3& mom)
{
   double xmom,ymom,zmom;
   double ptotal;
   dio_spec->GetRandom(ptotal);

   gRandom->Sphere(xmom,ymom,zmom,ptotal);
   mom.SetXYZ(xmom,ymom,zmom);
}
/* Spectrum is cited from Physical Review D 84, 013006 (2011) */
void get_init_mom_DIO(TVector3& mom)
{
   //double e_startpoint = 0.103; // GeV/c
   double e_startpoint = g_dio_startpoint;//0.101000; // GeV/c
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
   //printf("xmom %lf ymom %lf zmom %lf\n",xmom,ymom,zmom);
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
/* Radiative Pion Capture
 * Paper: Photon Spectra from Radiative Absorption of Pions in Nuclei
 * Physical Review C, Volume 5, Number 6 June 1972
 * Figure 7 (a)
 * */
TH1F* read_RPC_data(char* fname)
{
   FILE* fp = fopen(fname,"r");
   if (fp==NULL) {
      fprintf(stderr,"ERROR: failed to open RPC_data %s\n",fname);
      exit(1);
   }
   TH1F* h1 = new TH1F("hist_RPC","hist_RPC",74,60,133);
   char line[128];
   double mom,cont;
   while (fgets(line,sizeof(line),fp)) {
      sscanf(line,"%lf %lf",&mom,&cont);
      h1->Fill(mom,cont);
   }
   fclose(fp);

   //TCanvas*c1 = new TCanvas("c1","c1");
   //h1->Draw();
   //c1->Print("hist_RPC.png");

   return h1;
}
void get_init_mom_RPC(TVector3& mom)
{
   double ptotal = g_hist_RPC->GetRandom();
   printf("get_init_mom_RPC: ptotal %lf\n",ptotal);

   double xmom,ymom,zmom;
   gRandom->Sphere(xmom,ymom,zmom,ptotal);
   printf("get_init_om_RPC: xmom %lf ymom %lf zmom %lf\n",xmom,ymom,zmom);
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
double root_hits_tofg[MAX_HIT];
double root_hits_posx[MAX_HIT];
double root_hits_posy[MAX_HIT];
double root_hits_posz[MAX_HIT];
double root_hits_momx[MAX_HIT];
double root_hits_momy[MAX_HIT];
double root_hits_momz[MAX_HIT];
TFile* root_f;
TTree* root_t;
int root_detid[MAX_HIT];
int root_trkid[MAX_HIT];
int root_copyno[MAX_HIT];
int root_pdgid[MAX_HIT];

// additional values for wire closest approach
int    root_wire_nhits;
int    root_wire_ilayer[MAX_WIRE_HIT];
int    root_wire_icell[MAX_WIRE_HIT];
double root_wire_rdrift[MAX_WIRE_HIT];
double root_wire_zreco[MAX_WIRE_HIT];
double root_wire_hitx[MAX_WIRE_HIT];
double root_wire_hity[MAX_WIRE_HIT];
double root_wire_length[MAX_WIRE_HIT];
double   root_wire_posx[MAX_WIRE_HIT];
double   root_wire_posy[MAX_WIRE_HIT];
// added (2013/3/14)
double root_wire_momx[MAX_WIRE_HIT];
double root_wire_momy[MAX_WIRE_HIT];
double root_wire_momz[MAX_WIRE_HIT];
// added (2013/4/12)
int root_wire_is_reverse_cell[MAX_WIRE_HIT];
// added (2013/6/2)
double root_wire_tofg[MAX_WIRE_HIT];

double root_exit_xpos_in_target;
double root_exit_ypos_in_target;
double root_exit_zpos_in_target;
double root_exit_xmom_in_target;
double root_exit_ymom_in_target;
double root_exit_zmom_in_target;
double root_exit_amom_in_target;
double root_edep_in_target;

/* 
 * Proton data
 *    Check how tracking will be effected by rejecting cellls on which protons hit
 * 
 * */
void read_proton_data(char* proton_data)
{
   TFile* f = new TFile(proton_data);
   g_tree_proton= (TTree*)f->Get("t");
   if (g_tree_proton==NULL) {
      fprintf(stderr,"ERROR: tree was not found in (%s)\n", proton_data);
      exit(1);
   }
   g_tree_proton->SetBranchAddress("ntracks", &g_proton_ntracks);
   g_tree_proton->SetBranchAddress("nallhits", &g_proton_nallhits);
   g_tree_proton->SetBranchAddress("ilayers", g_proton_ilayers);
   g_tree_proton->SetBranchAddress("icells", g_proton_icells);
   g_tree_proton->SetBranchAddress("iruns", g_proton_iruns);
   g_tree_proton->SetBranchAddress("ievs", g_proton_ievs);
   g_tree_proton->SetBranchAddress("xwires", g_proton_xwires); // x wire position at given z
   g_tree_proton->SetBranchAddress("ywires", g_proton_ywires); // y wire position at given z
   g_proton_cur_iev = -1;
   g_proton_entries = g_tree_proton->GetEntries();
}
void count_overwap_cell(int iev, int last_idx)
{
   // use global varialbels (root_wire_XXX, come from input root file)
   
   int ov_ilayers_dup[MAX_HIT];
   int ov_icells_dup[MAX_HIT];
   int ov_iruns_dup[MAX_HIT];
   int ov_ievs_dup[MAX_HIT];
   double ov_xwires_pro_dup[MAX_HIT];
   double ov_ywires_pro_dup[MAX_HIT];
   double ov_xwires_sig_dup[MAX_HIT];
   double ov_ywires_sig_dup[MAX_HIT];

   g_proton_cur_iev++;
   g_tree_proton->GetEntry(g_proton_cur_iev);

   if (iev==370) {
   for (int isig=0; isig<root_wire_nhits; isig++) {
      printf("isig %d wire_ilayer %d wire_icell %d\n",isig,root_wire_ilayer[isig],root_wire_icell[isig]);
   }
   for (int iproton=0; iproton<g_proton_nallhits; iproton++) {
      printf("iproton %d proton_ilayers %d proton_icells %d\n",iproton,g_proton_ilayers[iproton], g_proton_icells[iproton]);
   }
   }


   int nhits_dup = 0;
   for (int isig=0; isig<=last_idx; isig++) {
      for (int iproton=0; iproton<g_proton_nallhits; iproton++) {
         if (root_wire_ilayer[isig] == g_proton_ilayers[iproton]) {
            if (root_wire_icell[isig] == g_proton_icells[iproton]) {
               //if (nhits_dup>=max) {
               //   fprintf(stderr,"Warning: iev %d count_overwap_cell too many overwaps. (max=%d), just neglect remaining overwaps\n",iev,max);
               //   goto count_ov_uniq;
               //}
               ov_ilayers_dup[nhits_dup] = root_wire_ilayer[isig];
               ov_icells_dup[nhits_dup] = root_wire_icell[isig];
               ov_iruns_dup[nhits_dup] = g_proton_iruns[iproton];
               ov_ievs_dup[nhits_dup]  = g_proton_ievs[iproton];
               ov_xwires_pro_dup[nhits_dup]  = g_proton_xwires[iproton];
               ov_ywires_pro_dup[nhits_dup]  = g_proton_ywires[iproton];
               ov_xwires_sig_dup[nhits_dup]  = root_wire_posx[isig];
               ov_ywires_sig_dup[nhits_dup]  = root_wire_posy[isig];
               //if (iev==40||iev==52) {
               //   printf("===> dup: iev %d iev_in_root %d iproton %5d iruns %d ievs %5d ilayer %d icells %d\n", 
               //         iev, g_proton_cur_iev,iproton, g_proton_iruns[iproton], g_proton_ievs[iproton], root_wire_ilayer[isig], root_wire_icell[isig]);
               //}
               nhits_dup++;
            }
         }
      }
   }

count_ov_uniq:
   // nhits can be duplicated, so uniq

   int nhits_uniq=0;
   for (int i=0; i<nhits_dup; i++) {

      bool this_is_duplicated = false;
      for (int j=i+1; j<nhits_dup; j++) {
         if (ov_ilayers_dup[i] == ov_ilayers_dup[j] && ov_icells_dup[i] == ov_icells_dup[j])  {
            this_is_duplicated = true;
            break;
         }
      }
      if (!this_is_duplicated) {
         //if (iev==40||iev==52) {
         //   printf("===> uniq: iev %d iev_in_root %d i %5d iruns %d ievs %5d ilayer %d icells %d\n", 
         //         iev, g_proton_cur_iev, i, ov_iruns_dup[i], ov_ievs_dup[i], ov_ilayers_dup[i], ov_icells_dup[i]);
         //}
         g_proton_ov_iruns_uniq[nhits_uniq] = ov_iruns_dup[i];
         g_proton_ov_ievs_uniq[nhits_uniq] = ov_ievs_dup[i];
         g_proton_ov_ilayers_uniq[nhits_uniq] = ov_ilayers_dup[i];
         g_proton_ov_icells_uniq[nhits_uniq] = ov_icells_dup[i];
         g_proton_ov_xwires_pro_uniq[nhits_uniq] = ov_xwires_pro_dup[i];
         g_proton_ov_ywires_pro_uniq[nhits_uniq] = ov_ywires_pro_dup[i];
         g_proton_ov_xwires_sig_uniq[nhits_uniq] = ov_xwires_sig_dup[i];
         g_proton_ov_ywires_sig_uniq[nhits_uniq] = ov_ywires_sig_dup[i];
         nhits_uniq++;
      }
   }

count_ov_end:
   printf("count_overwap_cell:  iev %d root_wire_nhits %d #/tracks %d #/allhits %d #/overwap (dup: %d , uniq: %d )\n", 
         iev, root_wire_nhits,g_proton_ntracks,g_proton_nallhits,nhits_dup, nhits_uniq);

   if (g_proton_cur_iev >= g_proton_entries-1) {
      g_proton_cur_iev = -1;
   }
   g_proton_ov_ncells_uniq = nhits_uniq;
}


int get_disk_number(TVector3& posini)
{
   double z = posini.Z();
   return (int)( (z+1)/5 + 9);
}
int find_index_of_90deg_arc(int iev)
{
   int idx=-1;
   for (int i=1; i<root_wire_nhits; i++) {
      if (root_wire_ilayer[i] - root_wire_ilayer[i-1] < 0) {
         idx=i-1;
         break;
      }
   }
   return idx;
}
int find_index_of_180deg_arc(int iev)
{
   int last_idx=-1;
   int ilayer1_idx =-1;
   //printf("===> i root_wire_nhits %d\n", root_wire_nhits);
   for (int i=0; i<root_wire_nhits; i++) {
      //printf("===> i %d root_wire_ilayer %d root_icell %d\n", i,root_wire_ilayer[i], root_wire_icell[i]);
      if (root_wire_ilayer[i]!=0) {
         ilayer1_idx=i;
         break;
      }
   }
   if (ilayer1_idx==-1) {
      if (root_wire_nhits>0) {
         fprintf(stderr,"ERROR: all hits are in the first laye..., skip this event \n");
         ///fprintf(stderr,"ERROR: first hit is not found at ilayer=0... skip this event\n");
      }
      //last_idx=0;
      last_idx=-1; // fix bug (2013//4/7)
   } else {
      last_idx=root_wire_nhits-1; // if track pass away from CDC in middle layer, last_idx should be set to the last index of hit
      for (int i=ilayer1_idx; i<root_wire_nhits; i++) {
         if (iev==370) {
            printf("find_index_of_180deg_arc: i %d root_wire_ilayer %d root_wire_icell %d\n",i,root_wire_ilayer[i],root_wire_icell[i]);
         }
         if (root_wire_ilayer[i]==0) {
            /* Add check  (2013/05/15)
             * if cell number is larger than 10, then the cell is assumed that track hit end plate and come back CDC region 
             * and hit wrong cell in the 1st layer
             */
            if (TMath::Abs(root_wire_icell[i] - root_wire_icell[i-1]) > 10) {
               last_idx = i-1; // last hit should be one hit before
               break;
            }

            last_idx=i;
            break;
         }
      }
   }
   printf("read_hit: iev %d ilayer1_idx %d last_idx %d root_wire_nhits %d\n",iev,ilayer1_idx,last_idx,root_wire_nhits);
   return last_idx;
}
int read_hit(int iev, char* root_file, TVector3& posini, TVector3& momini)
{
   // ROOT file made by geant4_vmc
   //

   if (root_opened==0) {
      if (root_file==NULL)
         exit(1);
      root_f = new TFile(root_file);
      printf("root_file %s\n",root_file);
      root_t = (TTree*)root_f->Get("t");
      root_t->SetBranchAddress("ini_x_cm",&root_ini_posx);
      root_t->SetBranchAddress("ini_y_cm",&root_ini_posy);
      root_t->SetBranchAddress("ini_z_cm",&root_ini_posz);
      root_t->SetBranchAddress("ini_px_GeV",&root_ini_momx);
      root_t->SetBranchAddress("ini_py_GeV",&root_ini_momy);
      root_t->SetBranchAddress("ini_pz_GeV",&root_ini_momz);
      root_t->SetBranchAddress("nhits",&root_nhits);
      //root_t->SetBranchAddress("tof_ns",root_hits_tofg);
      root_t->SetBranchAddress("tofg_ns",root_hits_tofg);
      root_t->SetBranchAddress("x_cm",root_hits_posx);
      root_t->SetBranchAddress("y_cm",root_hits_posy);
      root_t->SetBranchAddress("z_cm",root_hits_posz);
      root_t->SetBranchAddress("px_GeV",root_hits_momx);
      root_t->SetBranchAddress("py_GeV",root_hits_momy);
      root_t->SetBranchAddress("pz_GeV",root_hits_momz);
      root_t->SetBranchAddress("detid",root_detid);
      root_t->SetBranchAddress("trackid",root_trkid);
      root_t->SetBranchAddress("copyno",root_copyno);
      root_t->SetBranchAddress("pdgid",root_pdgid);

      // values for hitpoint_type=Wire
      root_t->SetBranchAddress("nwirehit", &root_wire_nhits);
      root_t->SetBranchAddress("ilayer", root_wire_ilayer);
      root_t->SetBranchAddress("icell",  root_wire_icell);
      root_t->SetBranchAddress("dist", root_wire_rdrift);
      root_t->SetBranchAddress("minhit_z",  root_wire_zreco);
      //     additional infos
      root_t->SetBranchAddress("minhit_x",  root_wire_hitx);
      root_t->SetBranchAddress("minhit_y",  root_wire_hity);
      root_t->SetBranchAddress("length",  root_wire_length);
      // wirepos added (2013/3/3)
      root_t->SetBranchAddress("wire_x",  root_wire_posx);
      root_t->SetBranchAddress("wire_y",  root_wire_posy);
      // wiremom added (2013/3/14)
      root_t->SetBranchAddress("minhit_px",  root_wire_momx);
      root_t->SetBranchAddress("minhit_py",  root_wire_momy);
      root_t->SetBranchAddress("minhit_pz",  root_wire_momz);
      // added (2013/4/12)
      root_t->SetBranchAddress("is_reverse_cell",  root_wire_is_reverse_cell);
      // added (2013/6/2)
      root_t->SetBranchAddress("time",  root_wire_tofg);

      root_t->SetBranchAddress("wire_overflow",      &g_wire_overflow);
      root_t->SetBranchAddress("wire_seg_overflow",  &g_wire_seg_overflow);
      root_t->SetBranchAddress("wire_hit_overflow",  &g_wire_hit_overflow);

      // Hit position at exiting target
      root_t->SetBranchAddress("exit_xpos_in_target",&root_exit_xpos_in_target);
      root_t->SetBranchAddress("exit_ypos_in_target",&root_exit_ypos_in_target);
      root_t->SetBranchAddress("exit_zpos_in_target",&root_exit_zpos_in_target);
      root_t->SetBranchAddress("exit_xmom_in_target",&root_exit_xmom_in_target);
      root_t->SetBranchAddress("exit_ymom_in_target",&root_exit_ymom_in_target);
      root_t->SetBranchAddress("exit_zmom_in_target",&root_exit_zmom_in_target);
      root_t->SetBranchAddress("edep_in_target",&root_edep_in_target);
      //printf("root_edep_in_target %lf\n",root_edep_in_target);
      //exit(1);




      root_opened = 1;
   }

   int ret = root_t->GetEntry(iev);
   //printf("read_hit: ret %d\n",ret);
   if (ret == 0) { // no more events
      return -1;
   }
   // Set initial value at first hit in CDC
   // (2013/4/12)
#if 1
   posini.SetX(root_ini_posx);
   posini.SetY(root_ini_posy);
   posini.SetZ(root_ini_posz);
   momini.SetX(root_ini_momx);
   momini.SetY(root_ini_momy);
   momini.SetZ(root_ini_momz);
#endif

   if (strcmp(arg_extrap_pos,"TargetExit")==0) {
      posini.SetX(root_exit_xpos_in_target);
      posini.SetY(root_exit_ypos_in_target);
      posini.SetZ(root_exit_zpos_in_target);
      momini.SetX(root_exit_xmom_in_target);
      momini.SetY(root_exit_ymom_in_target);
      momini.SetZ(root_exit_zmom_in_target);
   }



   //printf("root_ini_posx %lf\n",root_ini_posx);
   //printf("root_ini_posy %lf\n",root_ini_posy);
   //printf("root_ini_posz %lf\n",root_ini_posz);
   //printf("root_ini_momx %lf\n",root_ini_momx);
   //printf("root_ini_momy %lf\n",root_ini_momy);
   //printf("root_ini_momz %lf\n",root_ini_momz);

#if 1 // comment-out (2013/05/20)
   /* use RKTrackRep just to calc stMCT,covMCT */
   GFAbsTrackRep* rephits = new RKTrackRep( posini, momini, PDGcode);
   stMCT->ResizeTo(rephits->getState());
   *stMCT = rephits->getState();
   covMCT->ResizeTo(rephits->getCov());
   *covMCT = rephits->getCov();

#endif






      //printf("root_nhits %d\n",root_nhits);
      g_hits_det.nhits = 0;
      g_hits_tgt.nhits = 0;
      g_hits_dum.nhits = 0;
      g_wire_nhits=0; // added 2013/4/7
      g_allwire_nhits=0; // added 2013/6/4 
      g_nhits_sense = 0; // added 2013/6/24
      g_nhits_field = 0; // added 2013/6/24

      int nhit_at_first_layer=0;
      int nhit_at_second_layer=0;
      bool finish_loop=false;


      printf("iev %d root_nhits      %d\n",iev,root_nhits);
      printf("iev %d root_wire_nhits %d\n",iev,root_wire_nhits);

      bool already_enter_chamber=false;
      bool already_gone_chamber=false;
      // Get hits only in first turn
      for (int i=0; i<root_nhits; i++) {
         //if (iev==2) {
         //printf("%d root_trkid %d root_detid %d\n",i, root_trkid[i], root_detid[i]);
         //}

         // Hits at target
         //if (root_trkid[i]==0 && root_detid[i]==5) {
         //   g_hits_tgt.posx[g_hits_tgt.nhits] = root_hits_posx[i];
         //   g_hits_tgt.posy[g_hits_tgt.nhits] = root_hits_posy[i];
         //   g_hits_tgt.posz[g_hits_tgt.nhits] = root_hits_posz[i];
         //   g_hits_tgt.momx[g_hits_tgt.nhits] = root_hits_momx[i];
         //   g_hits_tgt.momy[g_hits_tgt.nhits] = root_hits_momy[i];
         //   g_hits_tgt.momz[g_hits_tgt.nhits] = root_hits_momz[i];
         //   g_hits_tgt.tofg[g_hits_tgt.nhits] = root_hits_tofg[i];
         //   g_hits_tgt.nhits++;
         //}

         // Hits at scinti
         //if (root_trkid[i]==0 && root_detid[i]==12) {
         if (root_trkid[i]==0 && root_detid[i]==8) { // change detid for trig_scingi (2013/4/7)
            g_hits_abs.posx[g_hits_abs.nhits] = root_hits_posx[i];
            g_hits_abs.posy[g_hits_abs.nhits] = root_hits_posy[i];
            g_hits_abs.posz[g_hits_abs.nhits] = root_hits_posz[i];
            g_hits_abs.momx[g_hits_abs.nhits] = root_hits_momx[i];
            g_hits_abs.momy[g_hits_abs.nhits] = root_hits_momy[i];
            g_hits_abs.momz[g_hits_abs.nhits] = root_hits_momz[i];
            g_hits_abs.tofg[g_hits_abs.nhits] = root_hits_tofg[i];
            g_hits_abs.nhits++;
         }

         //printf("read_hit: root_nhits %d i %d root_trkid %d root_detid %d\n",root_nhits,i, root_trkid[i], root_detid[i]);

         if (strcmp(arg_hitpoint_type,"Layer")==0 || strcmp(arg_hitpoint_type,"Wire")==0 || strcmp(arg_hitpoint_type,"Wirepoint")==0) {

            // Hits at layers (only first turn is stored)
            //if (root_trkid[i]==0 && root_detid[i]==8 && nhit_at_first_layer<=1)  <<== 2012/10/10

            // In case of hitpoint_type is 'Wire', g_hits_det is used __only__ for initial guesses
            // so only first hit is ok to record, but it's useful for debug to see all of hits in chamber
            if (strcmp(arg_hitpoint_type,"Wire")==0||strcmp(arg_hitpoint_type,"Wirepoint")==0) {
               //exit(1);
               //if (root_trkid[i]==0 && root_detid[i]==8) { // change detid=8 => 2 (2013/4/7)
               if (root_trkid[i]==0 && root_detid[i]==2) {
                  //if (root_trkid[i]==0 && root_detid[i]==8) {
                  g_hits_det.posx[g_hits_det.nhits] = root_hits_posx[i];
                  g_hits_det.posy[g_hits_det.nhits] = root_hits_posy[i];
                  g_hits_det.posz[g_hits_det.nhits] = root_hits_posz[i];
                  g_hits_det.momx[g_hits_det.nhits] = root_hits_momx[i];
                  g_hits_det.momy[g_hits_det.nhits] = root_hits_momy[i];
                  g_hits_det.momz[g_hits_det.nhits] = root_hits_momz[i];
                  g_hits_det.tofg[g_hits_det.nhits] = root_hits_tofg[i];
                  //if (iev==2) {
                  //   printf("read_hit: [%d] posx %lf",g_hits_det.nhits,g_hits_det.posx[g_hits_det.nhits]);
                  //   printf(" posy %lf",g_hits_det.posy[g_hits_det.nhits]);
                  //   printf(" posz %lf\n",g_hits_det.posz[g_hits_det.nhits]);
                  //}

               g_hits_det.nhits++;
            }

            if (strcmp(arg_extrap_pos,"TargetExit")==0) {
               if (root_trkid[i]==0 && root_detid[i]==13) { // dummy volume for genfit to fit normally

                  if (root_detid[i]==8) {
                     already_enter_chamber=true;
                  }

                  // record until track enter drift chamber
                  if (!already_enter_chamber) {
                     g_hits_dum.posx[g_hits_dum.nhits] = root_hits_posx[i];
                     g_hits_dum.posy[g_hits_dum.nhits] = root_hits_posy[i];
                     g_hits_dum.posz[g_hits_dum.nhits] = root_hits_posz[i];
                     g_hits_dum.momx[g_hits_dum.nhits] = root_hits_momx[i];
                     g_hits_dum.momy[g_hits_dum.nhits] = root_hits_momy[i];
                     g_hits_dum.momz[g_hits_dum.nhits] = root_hits_momz[i];
                     g_hits_dum.tofg[g_hits_dum.nhits] = root_hits_tofg[i];
                     //printf("Hits in dummy volume x=%lf y=%lf z=%lf\n",root_hits_posx[i],root_hits_posy[i],root_hits_posz[i]);
                     g_hits_dum.nhits++;
                  }

               }
            }


         } else if (strcmp(arg_hitpoint_type,"Layer")==0 || strcmp(arg_extrap_pos,"TargetExit")==0) {
            if (root_trkid[i]==0 && root_detid[i]==8) {

               if (root_copyno[i]==1) {
                  if (already_enter_chamber) {
                     already_gone_chamber=true;
                  } else if (already_gone_chamber){
                     // Break because only record first turn
                     break;
                  }
               }

               g_hits_det.posx[g_hits_det.nhits] = root_hits_posx[i];
               g_hits_det.posy[g_hits_det.nhits] = root_hits_posy[i];
               g_hits_det.posz[g_hits_det.nhits] = root_hits_posz[i];
               g_hits_det.momx[g_hits_det.nhits] = root_hits_momx[i];
               g_hits_det.momy[g_hits_det.nhits] = root_hits_momy[i];
               g_hits_det.momz[g_hits_det.nhits] = root_hits_momz[i];
               g_hits_det.tofg[g_hits_det.nhits] = root_hits_tofg[i];
               g_hits_det.nhits++;

               already_enter_chamber=true;
            }
         }

         // The follong code make probmle... g_hits_abs is needed to be filled and those will be after layer hits,
         //  so don't break!! 2012/11/06
         //
         //printf("read_hit: A root_nhits %d i %d root_trkid %d root_detid %d\n",root_nhits,i, root_trkid[i], root_detid[i]);
         //if (already_enter_chamber && root_trkid[i]==0 && root_detid[i]!=8) {
         //   printf("read_hit: already_entr_chamber break at i %d g_hits_det.nhits %d\n", i, g_hits_det.nhits);
         //   break;
         //}

      }
   }

   //printf("g_hits_abs.nhits %d\n",g_hits_abs.nhits);

   int prev_ilayer=0;
   int diff_ilayer = 0;
   int prev_diff_ilayer = -1; // do not change!
   int is_incre_ilayer = 0;
   int is_dicre_ilayer = 0;

   /*
    *  Check number of hits in the first turn
    */
#if 0
   for (int i=1; i<root_wire_nhits; i++) {

      diff_ilayer = root_wire_ilayer[i] - root_wire_ilayer[i-1];
      if (diff_ilayer>0) {
         is_incre_ilayer = 1;
         if (is_dicre_ilayer==1) {
            last_idx = i;
            break;
            //printf("== break here %d\n",i);
         }
      }
      if (diff_ilayer<0) is_dicre_ilayer = 1;
      //printf("i %d diff %d incre %d ciscre %d\n",i,diff_ilayer,is_incre_ilayer,is_dicre_ilayer);
   }
   printf("~=== hoge1 iev %d last_idx %d \n",iev,last_idx);
   // check zpos 
   double diff_zpos = root_wire_zreco[1] - root_wire_zreco[0];
   if (diff_zpos<0) diff_zpos = -diff_zpos;
   for (int i=1; i<last_idx-1; i++)  {
      double diff = root_wire_zreco[i+1] - root_wire_zreco[i];
      if (diff<0) diff = -diff;
      printf("== iev %d i %d diff_zpos %lf diff %lf\n",iev,i,diff_zpos,diff);
      if (diff>diff_zpos*30) {
         last_idx = i-1;
      }
   }
   printf("~=== hoge2 iev %d last_idx %d \n",iev,last_idx);
#endif

/* BUG!! deleted (2013/3/24)
   int hit_more_than_one_layer=0;
   for (int i=1; i<root_wire_nhits; i++) {
      if (root_wire_ilayer[i] > root_wire_ilayer[0]) {
         hit_more_than_one_layer=1;
         continue;
      }
      if (root_wire_ilayer[i] == root_wire_ilayer[0]) {
         if (hit_more_than_one_layer==1) {
            last_idx = i;
            break;
         }
      }
   }
printf("read_hit:  iev %d root_wire_ilayer[0] %d hit_more_than_one_layer %d last_idx %d \n",iev,root_wire_ilayer[0],hit_more_than_one_layer,last_idx);
*/


   /*
    * Make uniq detid list (2013/6/2) Primary track only!!!
    * Assume secondaries will come after primary track hit
    */
   //for (int i=0; i<10000; i++) index_list[i] = 0;

   int last_primary_idx = root_nhits-1;
   for (int i=0; i<root_nhits; i++) {
      if (root_trkid[i]!=0) {
         last_primary_idx = i;
         break;
      }

   }
   g_uniq_npoint=0;
   for (int i=0; i<=last_primary_idx; i++) {
      if (g_uniq_npoint==0) {
         g_uniq_detid[g_uniq_npoint] = root_detid[i];
         g_uniq_tofg[g_uniq_npoint] = root_hits_tofg[i];
         g_uniq_npoint++;
      } else if (root_detid[i] != root_detid[i-1]) {
         g_uniq_detid[g_uniq_npoint] = root_detid[i];
         g_uniq_tofg[g_uniq_npoint] = root_hits_tofg[i];
         g_uniq_npoint++;
      }
   }


   /*
    * Determine last point used for fitting
    */
   g_last_idx_of_fitting =  find_index_of_180deg_arc(iev);
   //   int last_idx  = find_index_of_90deg_arc(iev);


   // calculate number of hits on wires
   g_end_time_of_1st_turn = root_wire_tofg[g_last_idx_of_fitting];
   printf("g_last_idx_of_fitting %d end_time_of_1st_turn %lf (ns)\n", g_last_idx_of_fitting,g_end_time_of_1st_turn);
   // record last_idx and end_time_of_1st_turn!!
   //
   for (int i=1; i<=g_end_time_of_1st_turn; i++) {
      // only for primary track (electron or proton?)
      if (root_trkid[i]!=0)
         continue;
      // suppress double count
      if (root_detid[i]==9  && root_detid[i-1]!=9)  {  // sense wires

         if (g_nhits_sense>=MAX_HIT_ON_WIRE) 
            continue;
         
         g_sense_hit_xpos[g_nhits_sense] = root_hits_posx[i];
         g_sense_hit_ypos[g_nhits_sense] = root_hits_posy[i];
         g_sense_hit_zpos[g_nhits_sense] = root_hits_posz[i];
         g_sense_hit_xmom[g_nhits_sense] = root_hits_momx[i];
         g_sense_hit_ymom[g_nhits_sense] = root_hits_momy[i];
         g_sense_hit_zmom[g_nhits_sense] = root_hits_momz[i];
         g_nhits_sense++; 
      }

      if (root_detid[i]==10 && root_detid[i-1]!=10) { // field wires

         if (g_nhits_field>=MAX_HIT_ON_WIRE) 
            continue;

         g_field_hit_xpos[g_nhits_field] = root_hits_posx[i];
         g_field_hit_ypos[g_nhits_field] = root_hits_posy[i];
         g_field_hit_zpos[g_nhits_field] = root_hits_posz[i];
         g_field_hit_xmom[g_nhits_field] = root_hits_momx[i];
         g_field_hit_ymom[g_nhits_field] = root_hits_momy[i];
         g_field_hit_zmom[g_nhits_field] = root_hits_momz[i];
         g_nhits_field++; 
      }

      //printf("iev %d i %d root_hits_tofg %lf\n", iev,i,root_hits_tofg[i]);


   }
   printf("iev %d g_nhits_(sense= %d, field= %d)\n", iev,g_nhits_sense, g_nhits_field);


   /*
    * If proton_data is set, reject cells where protons hit
    */

   g_proton_ov_ncells_uniq = 0;
   if (strcmp(arg_proton_data,"empty")!=0) {
      count_overwap_cell(iev, g_last_idx_of_fitting);
   }



   double tofg_at_end_of_1st_turn = root_wire_tofg[g_last_idx_of_fitting];
   printf("iev %d root_nhits %d g_last_idx_of_fitting %d\n", iev, root_nhits, g_last_idx_of_fitting);

   // Set Track length in CDC
   g_wire_track_length = -100;
   if (g_last_idx_of_fitting>0) {
      g_wire_track_length = root_wire_length[g_last_idx_of_fitting] - root_wire_length[0];
   }

   double max_abs_diff_wire_momz = -100;
   double max_abs_diff_wire_momt = -100;
   double max_abs_diff_wire_moma = -100;

   //double diff_zpos;
   double prev_diff_zpos;
   // store wire information for hitpoint_type=Wire
   // Only hit information at the 1st turn is recored.
   //
      int nskip_by_rdrift=0;
      int nskip_by_proton=0;
   if (strcmp(arg_hitpoint_type,"Wire")==0 || strcmp(arg_hitpoint_type,"Wirepoint")==0) {
      printf("root_wire_nhits %d (arg_hitpoint_type %s)\n",root_wire_nhits,arg_hitpoint_type);


      // Reject first hit at the chamber, since rdist is wrong, => debug later 
      //int n_max = (root_wire_nhits>10) ? 21: root_wire_nhits;
      //      for (int i=3; i<n_max; i++) {
      //      for (int i=1; i<n_max; i++) {

      //      for (int i=0; i<n_max; i++) {


      for (int i=0; i<=g_last_idx_of_fitting; i++) {
         // DEBUG
         //if (iev==12) {
         //   if (g_wire_nhits==20 || g_wire_nhits==24) {
         //      continue;
         //   }
         //}

         int ilayer = root_wire_ilayer[i];
         int shift = g_config->shift[ilayer];
         if (shift==0) {
            g_nhits_axial++;
         } else {
            g_nhits_stereo++;
         }
#if 0
         // for 1T
         if (
               (root_wire_ilayer[i]>= 0 && root_wire_ilayer[i]<= 6) || 
               (root_wire_ilayer[i]>=12 && root_wire_ilayer[i]<=16) || 
               (root_wire_ilayer[i]>=22 && root_wire_ilayer[i]<=28) ) {
            g_nhits_stereo++;
         } else {
            g_nhits_axial++;
         }
         // for 1.5T
         //if (
         //      (root_wire_ilayer[i]>= 7 && root_wire_ilayer[i]<=11) || 
         //      (root_wire_ilayer[i]>=17 && root_wire_ilayer[i]<=21) || 
         //   g_nhits_stereo++;
         //}
#endif

         //     for (int i=0; i<root_wire_nhits; i++) {

         /*
          * Modify hit points to study chi2
          *
          *
          */
#if 0
         // large chi2 list
         if (
               (i== 0) ||
               (i== 8) ||
               (i== 9) ||
               (i==13) ||
               (i==33) ||
               (i==41) ||
               (i==46) ||
               (i==66) ||
               (i==86) ||
               (i==98) ||
               (i==108) )
            continue;
#endif



         //       for (int i=0; i<root_wire_nhits; i++) {

         //for (int i=0; i<root_wire_nhits; i++) {

         //printf("read_hit: wirehit: i %d ilayer %d icell %d rdrift %lf zreco %lf\n",i, root_wire_ilayer[i], root_wire_icell[i], root_wire_rdrift[i], root_wire_zreco[i]);

         /*
            if (root_wire_ilayer[i]==0) { // <= hits at ilayer 0 is not single (hit at cell0,1 e.g.), so check hit at 2nd layer
            nhit_at_first_layer++;
            if (nhit_at_second_layer>0) {
            finish_loop=true;
            }
            }
            if (root_wire_ilayer[i]==1) {
            nhit_at_second_layer++;
            }
            */


         /*
          *
          *  2013/3/15
          *   First hit is not always layer=0, since there is case that 
          *   track pass through corner of cell in the first layer, so skipped.
          *   Algorithm is that ilayer is increasing at 2 times
          *
          */ 
         /* original code */
         //   if (prev_ilayer==1 && root_wire_ilayer[i]==0) {
         //   g_wire_track_length = root_wire_length[i] - root_wire_length[0];
         //   break;
         //   }


         /* second version  */

         //diff_ilayer = root_wire_ilayer[i] - prev_ilayer;
         //if (iev==10) printf("==hoge iev %d ilayer %d prev_ilayer %d\n",iev,root_wire_ilayer[i],prev_ilayer);
         //if (prev_diff_ilayer < 0 && diff_ilayer > 0) {
         //   num_incre_ilayer++;
         //   if (num_incre_ilayer==2) {
         //      g_wire_nhits -= 1;
         //      if (g_wire_nhits<0) g_wire_nhits = 0;
         //      g_wire_track_length = root_wire_length[i] - root_wire_length[0];
         //      break;
         //   }
         //}
         //prev_diff_ilayer = diff_ilayer;
         //prev_ilayer = root_wire_ilayer[i];


         //DEBUG
         //if (i==5) continue;

         // added (2013/4/12)
//         if (root_wire_is_reverse_cell[i]==1) {
//            continue;
//         }

         // Before skipping cells, fill all of cell
         g_allwire_cut_type[g_allwire_nhits] = 0;
         g_allwire_ilayer[g_allwire_nhits] = root_wire_ilayer[i];
         g_allwire_icell [g_allwire_nhits] = root_wire_icell[i];
         g_allwire_rdrift[g_allwire_nhits] = root_wire_rdrift[i];
         g_allwire_zreco [g_allwire_nhits] = root_wire_zreco[i];
         g_allwire_hitx  [g_allwire_nhits] = root_wire_hitx[i];
         g_allwire_hity  [g_allwire_nhits] = root_wire_hity[i];
         g_allwire_posx  [g_allwire_nhits] = root_wire_posx[i];
         g_allwire_posy  [g_allwire_nhits] = root_wire_posy[i];
         g_allwire_momx  [g_allwire_nhits] = root_wire_momx[i];
         g_allwire_momy  [g_allwire_nhits] = root_wire_momy[i];
         g_allwire_momz  [g_allwire_nhits] = root_wire_momz[i];
         g_allwire_tofg  [g_allwire_nhits] = root_wire_tofg[i];
         g_allwire_nhits ++;

         /*
          * Skip hit (cell) here, before filling g_wire_XXX (2013/6/3)
          *
          * 1) proton hit cell
          * 2) small rdrift
          */

         /* Proton */
         //printf("==> A) Skip due to overwap by proton hit (iev %d) ov_nhits %d\n", iev, ov_nhits);
         bool skip_this_cell_by_proton=false;
         for (int j=0; j<g_proton_ov_ncells_uniq; j++) {
            if (root_wire_ilayer[i] == g_proton_ov_ilayers_uniq[j] && root_wire_icell[i] == g_proton_ov_icells_uniq[j]) {
               //printf("====> Skip due to overwap by proton hit (iev %d ilayer %d  icell %d) g_wire_nhits %d\n", 
               //      iev, g_proton_ov_ilayers_uniq[j], g_proton_ov_icells_uniq[j], g_wire_nhits);
               skip_this_cell_by_proton=true;
               break;
            }
         }

         if (root_wire_rdrift[i]<0) {
            fprintf(stderr,"Warning!! root_wire_rdrift < 0, so skip this hit...\n");
            continue;
         }


         /* R smearing */
         g_rdrift_err = 0.02; // 200um 
//         g_rdrift_err = 0.0001; // debug (2013/6/25)



//         g_rdrift_err = 0.01; // 100um 
         double rdrift_smeared;
         while (1) {
            rdrift_smeared = gRandom->Gaus(root_wire_rdrift[i], g_rdrift_err);
            if (rdrift_smeared>0) break;
         }
         g_wire_rdrift_smeared[g_wire_nhits] = rdrift_smeared;
         g_allwire_rdrift_smeared[g_allwire_nhits-1] = rdrift_smeared;

         if (strcmp(arg_proton_data,"empty")!=0) {
            if (skip_this_cell_by_proton) {
               g_allwire_cut_type[g_allwire_nhits-1] = 1;
               nskip_by_proton++;
               continue;
            }
         }

#if 0
//         double rdrift_min = 0.003;

         double rdrift_min = 0.06; // 2103/6/15
//         if (root_wire_rdrift[i]<rdrift_min) { // comment-out (2013/6/4) after run513
         if (rdrift_smeared < rdrift_min) {
            g_allwire_cut_type[g_allwire_nhits-1] = 2;
            nskip_by_rdrift++;
            continue;
         }
#endif
#if 0
         double rdrift_max = 0.7;
         if (rdrift_smeared > rdrift_max) {
            g_allwire_cut_type[g_allwire_nhits-1] = 2;
            nskip_by_rdrift++;
            continue;
         }
#endif

         /*
         fprintf(stderr,"i %d root_wire_ilayer %d root_wire_icell %d root_wire_rdrift %lf g_wire_nhits %d\n", 
               i, root_wire_ilayer[i], root_wire_icell[i], root_wire_rdrift[i], g_wire_nhits);
               */

         // used for wire hits
         g_wire_ilayer[g_wire_nhits] = root_wire_ilayer[i];
         g_wire_icell [g_wire_nhits] = root_wire_icell[i];
         g_wire_rdrift[g_wire_nhits] = root_wire_rdrift[i];
         g_wire_zreco [g_wire_nhits] = root_wire_zreco[i];
         g_wire_hitx  [g_wire_nhits] = root_wire_hitx[i];
         g_wire_hity  [g_wire_nhits] = root_wire_hity[i];
         g_wire_posx  [g_wire_nhits] = root_wire_posx[i];
         g_wire_posy  [g_wire_nhits] = root_wire_posy[i];
         g_wire_momx  [g_wire_nhits] = root_wire_momx[i];
         g_wire_momy  [g_wire_nhits] = root_wire_momy[i];
         g_wire_momz  [g_wire_nhits] = root_wire_momz[i];
         g_wire_tofg  [g_wire_nhits] = root_wire_tofg[i];

         // Check Momemtum difference (Scattering)
         if (g_wire_nhits>0) {
            // pz
            double momz1 = g_wire_momz[g_wire_nhits];
            double momz2 = g_wire_momz[g_wire_nhits-1];
            double cur_abs_diff_wire_momz = TMath::Abs(momz1-momz2);
            if (max_abs_diff_wire_momz < cur_abs_diff_wire_momz) {
               max_abs_diff_wire_momz  =  cur_abs_diff_wire_momz;
               g_max_diff_pz = max_abs_diff_wire_momz;
               g_max_diff_pz_iter = g_wire_nhits;
               g_max_diff_pz_ilayer = g_wire_ilayer[g_wire_nhits];
            }

            // pt
            double momt1 = sqrt2(g_wire_momx[g_wire_nhits], g_wire_momy[g_wire_nhits]);
            double momt2 = sqrt2(g_wire_momx[g_wire_nhits-1], g_wire_momy[g_wire_nhits-1]);
            double cur_abs_diff_wire_momt = TMath::Abs(momt1-momt2);
            if (max_abs_diff_wire_momt < cur_abs_diff_wire_momt) {
               max_abs_diff_wire_momt  =  cur_abs_diff_wire_momt;
               g_max_diff_pt = max_abs_diff_wire_momt;
               g_max_diff_pt_iter = g_wire_nhits;
               g_max_diff_pt_ilayer = g_wire_ilayer[g_wire_nhits];
            }

            // pa
            double moma1 = sqrt3(g_wire_momx[g_wire_nhits], g_wire_momy[g_wire_nhits], g_wire_momz[g_wire_nhits]);
            double moma2 = sqrt3(g_wire_momx[g_wire_nhits-1], g_wire_momy[g_wire_nhits-1], g_wire_momz[g_wire_nhits-1]);
            double cur_abs_diff_wire_moma = TMath::Abs(moma1-moma2);
            if (max_abs_diff_wire_moma < cur_abs_diff_wire_moma) {
               max_abs_diff_wire_moma  =  cur_abs_diff_wire_moma;
               g_max_diff_pa = max_abs_diff_wire_moma;
               g_max_diff_pa_iter = g_wire_nhits;
               g_max_diff_pa_ilayer = g_wire_ilayer[g_wire_nhits];
            }

         }


         if (g_wire_nhits==0) {
//         if (g_wire_nhits==1) {
            g_wire_ini_posx = root_wire_hitx[i];
            g_wire_ini_posy = root_wire_hity[i];
            g_wire_ini_posz = root_wire_zreco[i];
            g_wire_ini_momx = root_wire_momx[i];
            g_wire_ini_momy = root_wire_momy[i];
            g_wire_ini_momz = root_wire_momz[i];
         }
         if (iev==22) {
            printf("read_hit: [%d] wire_hitx %lf",g_wire_nhits,root_wire_hitx[g_wire_nhits]);
            printf(" wire_hity %lf",root_wire_hity[g_wire_nhits]);
            printf(" wire_zreco %lf",root_wire_zreco[g_wire_nhits]);
            printf(" wire_momx %lf",root_wire_momx[g_wire_nhits]);
            printf(" wire_momy %lf",root_wire_momy[g_wire_nhits]);
            printf(" wire_momz %lf\n",root_wire_momz[g_wire_nhits]);
         }
         g_wire_nhits++;

         //printf("=====>>>> iev %d i %d g_wire_nhits %d\n",iev, i, g_wire_nhits);
         //printf("wirehit: i %d ilayer %d icell %d rdrift %lf zreco %lf\n",i, root_wire_ilayer[i], root_wire_icell[i], root_wire_rdrift[i], root_wire_zreco[i]);
         //prev_ilayer = root_wire_ilayer[i];
      }
      }

      /*
       *
       * 2012/11/06
       * Note that if track pass through only 1st layer of inside layer (tube6?, not recognized by detid==8), 
       * and g_wire_nhits>0. Since initial momentum is set by g_hits_det.mom(x,y,z)[0], which is 0 at initail,
       * this case causes error at GeaneTrackRep2.
       * This event is not need to be fiited. (check after read_hit and before fitting) 
       */



#if 0 // comment out (2013/3/12) to see the difference
      /*
       * 2012/11/07
       * It seems that if there is a gap between hits (because hit is corner .. rhit > min_radius=0.5cm)
       * genfit fails.
       * So, only use (collect) hits until those gap is appear.
       */
      if (g_wire_nhits>0) {
         double prev_zreco = g_wire_zreco[0];
         for (int i=1; i<g_wire_nhits; i++) {
            if (abs(g_wire_zreco[i]-prev_zreco) > 1 /* cm*/ ) {
               fprintf(stderr,"g_wire_nhits is overwritten  since there is gap!!! (#/hits %d -> %d)\n",g_wire_nhits,i);
               g_wire_nhits = i-1; /* do not include overgapped point */
               break;
            }
            prev_zreco = g_wire_zreco[i];
         }
      }
#endif


      //      // dump hits in dummy volume
      //      for (int i=0; i<g_hits_dum.nhits; i++) {
      //         printf("hits in dummy volume: %d x %lf y %lf z %lf\n",
      //               i,g_hits_dum.posx[i],
      //               i,g_hits_dum.posy[i],
      //               i,g_hits_dum.posz[i]);
      //
      //      }

      // check maximum ilayer (2013/3/13)
      g_max_ilayer = -1;
      for (int i=0; i<g_wire_nhits; i++) {
         if (g_max_ilayer<g_wire_ilayer[i]) {
            g_max_ilayer = g_wire_ilayer[i];
         }
      }

      g_max_dist = -1;
      for (int i=1; i<g_wire_nhits; i++) {
         double dist = root_wire_length[i] - root_wire_length[i-1];
         if (g_max_dist<dist) {
            g_max_dist = dist;
         }
      }

      //debug
//      if (iev==3) {
//         g_wire_nhits = 40;
//      }


#if 0
// Swap hit for fitting test (2013/4/5)
      if (iev==11) {
         for (int i=0; i<g_wire_nhits; i++) {
            //swap ihit 10 <=> 20
            if (i==10) {
               g_wire_ilayer[i] = g_wire_ilayer[20];
               g_wire_icell[i]  = g_wire_icell[20];
               g_wire_rdrift[i] = g_wire_rdrift[20];
            } else if (i==20) {
               g_wire_ilayer[i] = g_wire_ilayer[10];
               g_wire_icell[i]  = g_wire_icell[10];
               g_wire_rdrift[i] = g_wire_rdrift[10];
            }
         }
      }
#endif

//      // Set initial value 
//   posini.SetX(g_wire_ini_posx);
//   posini.SetY(g_wire_ini_posy);
//   posini.SetZ(g_wire_ini_posz);
//   momini.SetX(g_wire_ini_momx);
//   momini.SetY(g_wire_ini_momy);
//   momini.SetZ(g_wire_ini_momz);

      /*
       *  Check end plate hit events (for res_pa>0 tail cut study) 2013/6/2
       *
       */
      g_hit_endplate_nhits=0;
      for (int i=0; i<g_uniq_npoint; i++) {
         if (g_uniq_detid[i]==7) {
            if (g_hit_endplate_nhits>=MAX_ENDPLATE_HIT) {
               fprintf(stderr,"Warning: g_hit_endplate_nhits overflowed (max=%d), break anyway\n", MAX_ENDPLATE_HIT);
               break;
            }
            g_hit_endplate_tofg[g_hit_endplate_nhits] = g_uniq_tofg[i];
            g_hit_endplate_nhits++;
         }
      }

      fprintf(stderr,"END OF read_hit: g_hits_det.nhits %d g_wire_nhits %d\n",g_hits_det.nhits, g_wire_nhits);

      printf("read_hit: allwire_nhits %d wire_nhits %d nskip_by_proton %d nskip_by_rdrift %d\n",g_allwire_nhits, g_wire_nhits, nskip_by_proton,nskip_by_rdrift);
      printf("read_hit: g_hits_det.nhits %d g_wire_nhits %d\n",g_hits_det.nhits, g_wire_nhits);
      printf("read_hit: g_wire_ini_posz %lf\n", g_wire_ini_posz);
      printf("read_hit: g_uniq_npoint %d\n", g_uniq_npoint);
      for (int i=0; i<g_uniq_npoint; i++) {
         printf("read_hit: i %d uniq_detid %d uniq_tofg %lf (ns)\n", i,g_uniq_detid[i], g_uniq_tofg[i]);
      }
      printf("read_hit: g_hit_endplate_nhits %d\n", g_hit_endplate_nhits);
      for (int i=0; i<g_hit_endplate_nhits; i++) {
         printf("read_hit: i %d g_hit_endplate_tofg %lf (ns)\n", i,g_hit_endplate_tofg[i]);
      }
      return 0;
   }

   int read_hit_old(int iev, char* root_file, TVector3& posini, TVector3& momini)
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

      //printf("root_nhits %d\n",root_nhits);
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
         fprintf(stderr,"read_iit_param: %s\n",fname);
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
int is_near_target(TVector3& pos)
{
   double x = pos.X();
   double y = pos.Y();
   double z = pos.Z();
   double r = sqrt2(x,y);
   if (r<g_target_radius) {
      /*
         if ((abs(z+40.0)<g_target_thickness+0.01)
         || (abs(z+35.0)<g_target_thickness+0.01) 
         || (abs(z+30.0)<g_target_thickness+0.01) 
         || (abs(z+25.0)<g_target_thickness+0.01) 
         || (abs(z+20.0)<g_target_thickness+0.01) 
      || (abs(z+15.0)<g_target_thickness+0.01) 
      || (abs(z+10.0)<g_target_thickness+0.01) 
      || (abs(z+05.0)<g_target_thickness+0.01) 
      || (abs(z+00.0)<g_target_thickness+0.01) 
      || (abs(z-40.0)<g_target_thickness+0.01) 
      || (abs(z-35.0)<g_target_thickness+0.01) 
      || (abs(z-30.0)<g_target_thickness+0.01) 
      || (abs(z-25.0)<g_target_thickness+0.01) 
      || (abs(z-20.0)<g_target_thickness+0.01) 
      || (abs(z-15.0)<g_target_thickness+0.01) 
      || (abs(z-10.0)<g_target_thickness+0.01) 
      || (abs(z-05.0)<g_target_thickness+0.01)) { //0.01cm = 0.1mm
      */
      //if (fmod(z/5.0)) {
       //  return 1;
      //}
      double absz = fabs(z);
      int z3 = absz/5.0 + 0.5;
      int z5 = z3*5.0;
      double z4 = fabs(z5-absz);
      //printf("z %lf z2 %lf z3 %d z4 %lf\n",absz, absz/5.0,z3,z4);
      if (z4<0.1) {
         return 1;
      }
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

   //const char* pre_volname=NULL;
   //const char* pre_volname = rephits->getVolName();
   const char* pre_volname = "dummy_volname";
   // get_hit_at_target is checking the hits of
   // TOPPER -> diskxxx pattern




   int ilayer;
   int ilayer_prev = -100;
   //for (int nhits=0; nhits<MAX_HIT; nhits++) {
   //for (int nhits=0; nhits<100; nhits++) {
   double step_size;
   step_size = g_target_thickness/2.0;
   while (1) {
      TVector3 pos = rephits->getPos();
      TVector3 mom = rephits->getMom();
      TVector3 dir = mom.Unit();
      //double step_size = 0.001; // cm
      //double step_size = 0.01; // cm
      
      //if (pos.Z()>0.01) {
      //   fprintf(stderr,"### mom.Mag %lf pos.X %lf pos.Y %lf pos.Z %lf\n",mom.Mag(),pos.X(),pos.Y(),pos.Z());
      //   break;
      //}
      
      //if (is_near_target(pos)==1) {
      if (sqrt2(pos.X(),pos.Y())>g_target_radius) {
         step_size = 0.1;
      }
      //step_size = 0.001; //g_target_thickness/2.0;
      //step_size = 0.0001; //g_target_thickness/2.0;


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
      fprintf(stderr,"cpos %lf %lf %lf\n",pos.X(),pos.Y(),pos.Z());
      fprintf(stderr,"npos %lf %lf %lf\n",next_pos.X(),next_pos.Y(),next_pos.Z());
      fprintf(stderr,"dir %lf %lf %lf\n",dir.X(),dir.Y(),dir.Z());

      try{
         tofg = rephits->extrapolate(d);
         //fprintf(stderr,"A\n");
      }
      catch(GFException& e){
         e.what();
         //fprintf(stderr,"B\n");
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

      //hits_at_straw.push_back(posR);
      //layers_at_straw.push_back(ilayer);

      g_hits_det.layer[g_hits_det.nhits] = ilayer;
      g_hits_det.posx[g_hits_det.nhits] = posR.X();
      g_hits_det.posy[g_hits_det.nhits] = posR.Y();
      g_hits_det.posz[g_hits_det.nhits] = posR.Z();
      g_hits_det.momx[g_hits_det.nhits] = momR.X();
      g_hits_det.momy[g_hits_det.nhits] = momR.Y();
      g_hits_det.momz[g_hits_det.nhits] = momR.Z();
      g_hits_det.tofg[g_hits_det.nhits] = tofg;

      //print_hits(g_hits_det.nhits,&g_hits_det);

      g_hits_det.nhits++;


      //printf("volname %s pre_volname %s nhits (tgt:%d abs:%d det:%d)\n",volname,pre_volname, g_hits_tgt.nhits, g_hits_abs.nhits,g_hits_det.nhits);

   }

   return 0;
}
int calc_pulls(int iev, GFAbsTrackRep *rep)
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

   if (iev>0) {
      printf("iev %d\n",iev);
   printf("invmom %lf\n",invmom);
   printf("sigmasqustate %lf\n",sigmasqustate);
   printf("sigma_p %lf\n",sigma_p);

   printf("momSi %lf\n",momSi);
   printf("momRe %lf\n",momRe);
   printf("momTr %lf\n",momTr);
   printf("momPu %lf\n",momPu);

   printf("qopSi %lf\n",qopSi);
   printf("qopRe %lf\n",qopRe);
   printf("qopTr %lf\n",qopTr);
   printf("qopPu %lf\n",qopPu);

   printf("upSi %lf\n",upSi);
   printf("upRe %lf\n",upRe);
   printf("upTr %lf\n",upTr);
   printf("upPu %lf\n",upPu);
   printf("vpSi %lf\n",vpSi);
   printf("vpRe %lf\n",vpRe);
   printf("vpTr %lf\n",vpTr);
   printf("vpPu %lf\n",vpPu);

   printf("uSi %lf\n",uSi);
   printf("uRe %lf\n",uRe);
   printf("uTr %lf\n",uTr);
   printf("uPu %lf\n",uPu);
   printf("vSi %lf\n",vSi);
   printf("vRe %lf\n",vRe);
   printf("vTr %lf\n",vTr);
   printf("vPu %lf\n",vPu);
   }

   return 0;
}

void print_plane(GFDetPlane& plane, const char* pre)
{
#if 1
   TVector3 orgin = plane.getO();
   TVector3 uaxis = plane.getU();
   TVector3 vaxis = plane.getV();
   fprintf(stderr,"%s: orgin: O=(%lf, %lf, %lf)\n",pre, orgin.X(), orgin.Y(), orgin.Z());
   fprintf(stderr,"%s: uaxis: O=(%lf, %lf, %lf)\n",pre, uaxis.X(), uaxis.Y(), uaxis.Z());
   fprintf(stderr,"%s: vrgin: O=(%lf, %lf, %lf)\n",pre, vaxis.X(), vaxis.Y(), vaxis.Z());
#endif
}

double do_the_fitting(int iev, struct tree_value* tv, GFTrack* fitTrack, GFAbsTrackRep* rep, GFDetPlane& plane)
{
   //
   // This will return chi2 value
   //

   GFKalman k;
   //GFDaf k;
   tv->error = 0;
   try {
      //   fprintf(stderr,"hoge fitting: A\n");
//         print_plane(plane,"before_processTrack: ");
      k.processTrack(fitTrack);
//         print_plane(plane,"after_processTrack: ");
      //   fprintf(stderr,"hoge fitting: B\n");
#if 0
         /*
          * Get chi2 at each hit points
          */
         for (int ihit=0; ihit<fitTrack.getNumHits(); ihit++) {

            GFAbsRecoHit* hit = fitTrack.getHit(ihit);
            double chi2 = k.getChi2Hit(hit,rep);
            fprintf(stderr,"fitting: ihit %d chi2 %lf\n",ihit,chi2);

         }

#endif


   } catch(GFException& e) {
      std::cerr << e.what();
      //fprintf(stderr,"iev %d hoge1\n",iev);
      std::cerr<< "Exceptions wont be further handled ->exit(1)  line " << __LINE__<<std::endl;
      tv->error = 100;
      return 100;
   }

   //   fprintf(stderr,"hoge fitting: C\n");
   tv->ndf = rep->getNDF();

   double no_extrap_fit_px; // fit px set before extrapolation (2013/4/14)
   double no_extrap_fit_py; // fit py set before extrapolation (2013/4/14)
   double no_extrap_fit_pz; // fit pz set before extrapolation (2013/4/14)
   double no_extrap_fit_pa; // fit pa set before extrapolation (2013/4/14)
   if (rep->getStatusFlag()==0) {
      //for (int j=0; j<g_nhits_fit; j++) {


      //TVector3 pos_ext(0.,0., g_hits_posz_smeared[j]);

      //GFDetPlane plane_ext(pos_ext,TVector3(1.,0.,0.),TVector3(0.,1.,0.));

      try {
//         print_plane(plane,"fitting: ");
         //GFDetPlane plane(pos,TVector3(1.,0.,0.),TVector3(0.,1.,0.));

#if 1
         /* 
          * Get Fitted values wihtout extraplateion  (2013/4/14)
           It seems that extraplated position is not correct (root343/hoge3.root, iev=22)
           */

    //     fprintf(stderr,"rep->exptrapolate ==hoge==\n");
    //     print_plane(plane,"");
    
         TMatrixT<double> state = rep->getState();
         double s0 = state[0][0];
         double s1 = state[1][0];
         double s2 = state[2][0];
         double s3 = state[3][0];
         double s4 = state[4][0];
         double s5 = state[5][0];

         double hpt = sqrt2(tv->hit.px,tv->hit.py);
         double hpa = sqrt3(tv->hit.px,tv->hit.py,tv->hit.pz);
         double fit_pa = -1./s0; // q=-1
         double hit_pa = hpa;
         double diff_pa = hit_pa - fit_pa;
         no_extrap_fit_px = fit_pa*s3;
         no_extrap_fit_py = fit_pa*s4;
         no_extrap_fit_pz = sqrt(fit_pa*fit_pa - no_extrap_fit_px*no_extrap_fit_px - no_extrap_fit_py*no_extrap_fit_py);
         no_extrap_fit_pa = fit_pa;

//         fprintf(stdout,"rep->getState() (q/p = %lf , u = %lf v = %lf u/w = %lf v/w = %lf spu = %lf) => hit_pa %lf (MeV/c) fit_pa %lf (MeV/c) diff_pa %lf (MeV/c)\n",
//               s0,s1,s2,s3,s4,s5, hit_pa*1000, fit_pa*1000, diff_pa*1000);
//         fprintf(stdout,"rep->getState() no_extrap: (px,py,pz) (%lf,%lf,%lf) pa %lf\n", 
//               no_extrap_fit_px,
//               no_extrap_fit_py,
//               no_extrap_fit_pz,
//               no_extrap_fit_pa);

#endif
         rep->extrapolate(plane);


      //double fit_x= fitTrack->getPos().X();
      //double fit_y= fitTrack->getPos().Y();
      //double fit_z= fitTrack->getPos().Z();

      //double fit_px= 1.0;
      //double fit_py= 1.0;
      //double fit_pz= 1.0;


      //TVector3 fit_pos =  TVector3(fit_x, fit_y, fit_z);
      //TVector3 fit_mom =  TVector3(fit_px, fit_py, fit_pz);
      //GFDetPlane plane_fitval = GFDetPlane(fit_pos, fit_mom.Unit());

      //   rep->extrapolate(plane_fitval);
//         rep->extrapolate(plane_ext);

      } catch(GFException& e) {
         e.what();
         //fprintf(stderr,"iev %d hoge2\n",iev);
         std::cerr<<"Exceptions wont be further handled ->exit(1)  line " << __LINE__<<std::endl;
         tv->error =2;
      }

      //tv->fit.x = rep->getPos().X();
      //tv->fit.y = rep->getPos().Y();
      //tv->fit.z = rep->getPos().Z();

      tv->fit.x = fitTrack->getPos().X();
      tv->fit.y = fitTrack->getPos().Y();
      tv->fit.z = fitTrack->getPos().Z();
      //tv->fit.x = fitTrack->getPos().X();
      //tv->fit.y = fitTrack->getPos().Y();
      //tv->fit.z = fitTrack->getPos().Z();

      tv->chi2 = rep->getRedChiSqu();
      //tv->chi2 = fitTrack->getRedChiSqu();

      tv->nfail = fitTrack->getFailedHits();
      //tv->nfail = fitTrack->getFailedHits();



      //printf("## hoge iev %d j %d chi2 %lf\n",iev,j,tv->chi2);
      //}
   } else {
      tv->error=3;
   }
   /* Do not set chi2 if error occured */
   if (tv->error!=0) {
      tv->chi2 = 9E10;
//      return -1; // rep is already in error
      return -200; // rep is already in error
   }
   // chi2-probability
   //tv->prob = TMath::Prob(tv->ndf*tv->chi2,tv->ndf);
   double chisq = tv->chi2*tv->ndf;
   tv->prob = TMath::Prob(chisq,tv->ndf);
   printf("ndf %d chisq %lf (chisq/ndf=%lf) prob %lf\n",tv->ndf,chisq,tv->chi2,tv->prob);

   tv->fit.px = fitTrack->getMom().X();
   tv->fit.py = fitTrack->getMom().Y();
   tv->fit.pz = fitTrack->getMom().Z();

//   tv->fit.px = no_extrap_fit_px;
//   tv->fit.py = no_extrap_fit_py;
//   tv->fit.pz = no_extrap_fit_pz;


   printf("fitTrack: hit.x,y,z = %lf %lf %lf\n",tv->hit.x,tv->hit.y,tv->hit.z);
   //tv->fit.px = fitTrack->getMom().X();
   //tv->fit.py = fitTrack->getMom().Y();
   //tv->fit.pz = fitTrack->getMom().Z();
   //tv->fit.px = rep->getMom().X();
   //tv->fit.py = rep->getMom().Y();
   //tv->fit.pz = rep->getMom().Z();

   //printf("loop_over_ifirst: iev %d tv->chi2 %lf\n",iev,retval);
   if (strcmp(arg_hit_type,"TXT")!=0) calc_pulls(iev,rep);
// comment-oub (2013/5/20)  if (strcmp(arg_hit_type,"TXT")!=0) calc_pulls(iev,rep);
   double retval =  tv->chi2;

   // only set if cur_value is smaller than prev_min
   
//double g_chi2_hit[6][1000]; // ipass,ihit
//double g_chi2_hit_max[6];
//double g_zfit_hit[6][1000]; // ipass,ihit
//double g_zfit_hit_max[6];
//double g_track_length_hit[6][1000]; // ipass,ihit
//double g_track_length_hit_max[6];
//double g_pa_hit[6][1000]; // ipass,ihit
//double g_pa_hit_rms[6];
//double g_pa_hit_min[6];
//double g_pa_hit_max[6];

   // Get chi2 at each hit points
   for (int ipass=0; ipass<6; ipass++) {
      g_chi2_hit_max[ipass] = -100;
   }
   g_chi2_nhits = k.chi2_nhits;
   for (int ipass=0; ipass<6; ipass++) {
      for (int ihit=0; ihit<g_chi2_nhits; ihit++) {
         g_chi2_hit[ipass][ihit] = k.chi2_hit[ipass][ihit];
         if (g_chi2_hit_max[ipass] < k.chi2_hit[ipass][ihit]) {
            g_chi2_hit_max[ipass] = k.chi2_hit[ipass][ihit];
         }
         //printf("ipass %d ihit %d chi2 %lf\n",ipass,ihit,g_chi2_hit[ipass][ihit]);
      }
   }

   // Get zfit at each hit points
   for (int ipass=0; ipass<6; ipass++) {
      g_zfit_hit_max[ipass] = -1000;
   }
   for (int ipass=0; ipass<6; ipass++) {
      for (int ihit=0; ihit<g_chi2_nhits; ihit++) {
         // check absolute z-position
         //printf("zfit_hit %lf\n", k.zfit_hit[ipass][ihit]);
         g_zfit_hit[ipass][ihit] = TMath::Abs(k.zfit_hit[ipass][ihit]);
         if (g_zfit_hit_max[ipass] < TMath::Abs(k.zfit_hit[ipass][ihit])) {
            g_zfit_hit_max[ipass] = TMath::Abs(k.zfit_hit[ipass][ihit]);
         }
         //printf("ipass %d ihit %d chi2 %lf\n",ipass,ihit,g_chi2_hit[ipass][ihit]);
      }
   }

   // Get track_length at each hit points
   for (int ipass=0; ipass<6; ipass++) {
      g_track_length_hit_max[ipass] = -1000;
   }
   for (int ipass=0; ipass<6; ipass++) {
      for (int ihit=0; ihit<g_chi2_nhits; ihit++) {
         g_track_length_hit[ipass][ihit] = k.track_length_hit[ipass][ihit];
         if (g_track_length_hit_max[ipass] < k.track_length_hit[ipass][ihit]) {
            g_track_length_hit_max[ipass] = k.track_length_hit[ipass][ihit];
         }
         //printf("ipass %d ihit %d chi2 %lf\n",ipass,ihit,g_chi2_hit[ipass][ihit]);
      }
   }

   // Get momenum at each hit points
   double pa_hit_mean[6];
   double pa_hit_mean2[6];
   double pz_hit_mean[6];
   double pz_hit_mean2[6];
   for (int ipass=0; ipass<6; ipass++) {
      pa_hit_mean[ipass] = 0;
      pa_hit_mean2[ipass] = 0;
      g_pa_hit_min[ipass]  = 1e10;
      g_pa_hit_max[ipass]  = -10;
      g_pa_hit_rms[ipass] = 1e10;

      pz_hit_mean[ipass] = 0;
      pz_hit_mean2[ipass] = 0;
      g_pz_hit_rms[ipass] = 1e10;
   }
   for (int ipass=0; ipass<6; ipass++) {
      for (int ihit=0; ihit<g_chi2_nhits; ihit++) {
         g_pa_hit[ipass][ihit] = k.pa_hit[ipass][ihit];
         g_pz_hit[ipass][ihit] = k.pz_hit[ipass][ihit];
         //if (ipass==0) {
         //   printf("ihit %d k.pa_hit[0][ihit] %lf\n", ihit, k.pa_hit[ipass][ihit]);
         //}
         pa_hit_mean[ipass] += g_pa_hit[ipass][ihit];
         pa_hit_mean2[ipass] += g_pa_hit[ipass][ihit]*g_pa_hit[ipass][ihit];
         pz_hit_mean[ipass] += g_pz_hit[ipass][ihit];
         pz_hit_mean2[ipass] += g_pz_hit[ipass][ihit]*g_pz_hit[ipass][ihit];

         if (g_pa_hit_max[ipass] < g_pa_hit[ipass][ihit]) {
            g_pa_hit_max[ipass]  =  g_pa_hit[ipass][ihit];
         }
         if (g_pa_hit_min[ipass] > g_pa_hit[ipass][ihit]) {
            g_pa_hit_min[ipass]  =  g_pa_hit[ipass][ihit];
         }
      }
      if (g_chi2_nhits>0) {
         pa_hit_mean[ipass] /= g_chi2_nhits;
         pa_hit_mean2[ipass] /= g_chi2_nhits;
         g_pa_hit_rms[ipass] = TMath::Sqrt(pa_hit_mean2[ipass] - pa_hit_mean[ipass]*pa_hit_mean[ipass]);

         pz_hit_mean[ipass] /= g_chi2_nhits;
         pz_hit_mean2[ipass] /= g_chi2_nhits;
         g_pz_hit_rms[ipass] = TMath::Sqrt(pz_hit_mean2[ipass] - pz_hit_mean[ipass]*pz_hit_mean[ipass]);
      }
   }
   //for (int ipass=0; ipass<6; ipass++) {
   //   printf("ipass %d pa_hit_mean2 %lf pa_hit_mean*pa_hit_mean %lf pa_hit_rms %lf g_chi2_nhits %d \n",
   //         ipass,pa_hit_mean2[ipass], pa_hit_mean[ipass]*pa_hit_mean[ipass], g_pa_hit_rms[ipass], g_chi2_nhits);
   //}




   //for (int ipass=0; ipass<6; ipass++) {
   //   printf("ipass %d zfit_max %lf\n",ipass,g_zfit_hit_max[ipass]);
   //}

   // Last plane

#if 0
   std::vector<double> vec_diffx;
   std::vector<double> vec_diffy;
   for (int i=0; i<g_nhits_fit; i++) {
      TVector3 v_pos;
      TVector3 v_mom;
      v_pos.SetXYZ(g_hits_posx_smeared[i],g_hits_posy_smeared[i],g_hits_posz_smeared[i]);
      v_mom.SetXYZ(g_hits_det.momx[i],g_hits_det.momy[i],g_hits_det.momz[i]);
      GFDetPlane p1 = GFDetPlane(v_pos,TVector3(1,0,0),TVector3(0,1,0));
      TVector3 a_pos;
      TVector3 a_mom;
      rep->getPosMom(p1, a_pos,a_mom);
      double diff_x  = v_pos.X()-a_pos.X();
      double diff_y  = v_pos.Y()-a_pos.Y();
      double diff_z  = v_pos.Z()-a_pos.Z();
      g_diffx += abs(diff_x);
      g_diffy += abs(diff_y);
      vec_diffx.push_back(g_diffx);
      vec_diffy.push_back(g_diffy);

      double diff_px  = v_mom.X()-a_mom.X();
      double diff_py  = v_mom.Y()-a_mom.Y();
      double diff_pz  = v_mom.Z()-a_mom.Z();
      /*
         if (iev==18788) printf("i %d a_pos %lf %lf %lf | %lf %lf %lf\n",i, a_pos.X(),a_pos.Y(),a_pos.Z(), diff_x,diff_y,diff_z);
         if (iev==18788) printf("i %d a_mom %lf %lf %lf %lf %lf| %lf %lf %lf %lf %lf\n",i, 
         v_mom.X(),v_mom.Y(),v_mom.Z(), sqrt2(v_mom.X(),v_mom.Y()),sqrt3(v_mom.X(),v_mom.Y(),v_mom.Z()), 
         a_mom.X(),a_mom.Y(),a_mom.Z(), sqrt2(a_mom.X(),a_mom.Y()),sqrt3(a_mom.X(),a_mom.Y(),a_mom.Z()));
         */
   }

   g_diffx /= g_nhits_fit;
   g_diffy /= g_nhits_fit;
   g_max_diffx = *max_element(vec_diffx.begin(),vec_diffx.end());
   g_max_diffy = *max_element(vec_diffy.begin(),vec_diffy.end());

#endif
#if 0
   GFDetPlane p2 = GFDetPlane(pos_last,TVector3(1,0,0),TVector3(0,1,0));
   TVector3 b_pos;
   TVector3 b_mom;
   try {
      rep->getPosMom(p2, b_pos,b_mom);
   } catch(GFException& e) {
      e.what();
      std::cerr<<"Exceptions wont be further handled ->exit(1)  line !! last_plane !!" << __LINE__<<std::endl;
      tv->error=100;
      return retval;
   }
   //   printf("b.pos %lf %lf %lf  b.mom %lf %lf %lf\n",
   //         b_pos.X(),b_pos.Y(),b_pos.Z(),
   //         b_mom.X(),b_mom.Y(),b_mom.Z());
   l_fit_x = b_pos.X();
   l_fit_y = b_pos.Y();
   l_fit_z = b_pos.Z();
   l_fit_px = b_mom.X();
   l_fit_py = b_mom.Y();
   l_fit_pz = b_mom.Z();
#endif

   //try {
   //   rep->extrapolate(last_plane);

   //} catch(GFException& e) {
   //   e.what();
   //   std::cerr<<"Exceptions wont be further handled ->exit(1)  line !! last_plane !!" << __LINE__<<std::endl;
   //   tv->error=100;
   //}
   //l_fit_x = fitTrack.getPos().X();
   //l_fit_y = fitTrack.getPos().Y();
   //l_fit_z = fitTrack.getPos().Z();
   //l_fit_px = fitTrack.getMom().X();
   //l_fit_py = fitTrack.getMom().Y();
   //l_fit_pz = fitTrack.getMom().Z();

   //delete rep; ~GFTrack delete rep;
   return retval;
}
/*
   int get_nhits_for_fitting2_for_wirehit()
   {
   int n = root_wire_nhits;
   int count=0;
   for (int i=0; i<n; i++) {
   int ilayer = root_wire_ilayer[i];
   if (ilayer==0) {
   count++;
   }
   if (count==2) {
   return i+1;
   }
   }
   return 0;
}
*/
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
         else if (strstr(line,"length"))    { sscanf_with_unit(line,"length %lf %s",&g_scinti_length); }
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
   printf("length %lf\n",g_scinti_length);

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
   TGeoMedium *mat_target = gGeoManager->GetMedium(g_target_material); assert(mat_target!=NULL);
   //printf("g_target_material %s\n",g_target_material);

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


   //bool debug_of_dummy=true;
   bool debug_of_dummy=false;

   if (debug_of_dummy) {
      rin_inwall = 12.0;
      rout_inwall = rin_inwall+wall_thick;
      fFirstPos = 13.0;
      g_chamber_num_layer = rin_outwall - fFirstPos - 2;
   }



   if (strcmp(g_target_type,"disk")==0 || strcmp(g_target_type,"test")==0) {
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
         // Fixed 2012/6/14
         //layers7[i] =  gGeoManager->MakeTube(Form("layer%d_7",i), strgas, din7, dout8, 154.0/2.0);
         layers7[i] =  gGeoManager->MakeTube(Form("layer%d_7",i), strgas, din7, dout7, 154.0/2.0);
         layers8[i] =  gGeoManager->MakeTube(Form("layer%d_8",i), strgas, din8, dout8, 154.0/2.0);

         if (layers7[i]) layers7[i]->SetLineColor(kGreen);
         if (layers7[i]) top->AddNode(layers7[i], 1, gGeoIdentity);
         // forgot layer6 and 8? (add 2012/6/14)
         if (layers6[i]) layers6[i]->SetLineColor(kFALSE);
         if (layers8[i]) layers8[i]->SetLineColor(kFALSE);
         if (layers6[i]) top->AddNode(layers6[i], 1, gGeoIdentity);
         if (layers8[i]) top->AddNode(layers8[i], 1, gGeoIdentity);


         printf("i %d r %lf din6 %lf dout6 %lf\n",i,r,din6,dout6);
         printf("i %d r %lf din7 %lf dout7 %lf\n",i,r,din7,dout7);
         printf("i %d r %lf din8 %lf dout8 %lf\n",i,r,din8,dout8);
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
      //fprintf(stderr,"end-cylinder is not implemented\n");
      //exit(1);

      double scinti_rpos = g_scinti_first_pos;
      double scinti_thick = g_scinti_thickness; //0.5; // 5mm
      double scinti_rmin = scinti_rpos - scinti_thick/2.0;
      double scinti_rmax = scinti_rmin + scinti_thick/2.0;
      double scinti_length = g_scinti_length;
      TGeoVolume* scinti_layer = gGeoManager->MakeTube("scinti_layer", scinti, scinti_rmin, scinti_rmax, scinti_length/2.0);
      scinti_layer->SetLineColor(kBlue);

      double off_end = 154.0/2.0 - scinti_length/2.0;
      // upstream
      top->AddNode(scinti_layer, 1, new TGeoTranslation(0,0,off_end));
      // downstream
      top->AddNode(scinti_layer, 1, new TGeoTranslation(0,0,-off_end));

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
      //TFile *outfile = TFile::Open("genfitGeom.root","RECREATE");
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
   strcpy(arg_output_root, "empty");
   strcpy(arg_input_root, "empty");
   strcpy(arg_input_txt, "empty");
   strcpy(arg_event_type,"SIG"); // SIG,DIO,Proton,TXT
   strcpy(arg_hit_type,"GeaneTrackRep2"); // RKTrackRep,GeaneTrackRep2,ROOT
   strcpy(arg_track_type,"GeaneTrackRep2"); // RKTrackRep,GeaneTrackRep2
   strcpy(arg_extrap_pos,"FirstHit"); //FirtHit,BeforeInwall
   strcpy(arg_hitpoint_type,"Layer"); //Layer, Wire
   strcpy(arg_proton_data, "empty");
   strcpy(arg_wire_config_fname, "empty");
   arg_fitting_ifirst = 0;
}

void print_usage(char* prog_name)
{
   fprintf(stderr,"Usage %s [inputs] [outputs]\n",prog_name);
   fprintf(stderr,"[inputs]\n");
   fprintf(stderr,"\t -s <seed>\n");
   fprintf(stderr,"\t -c <config.txt>\n");
   fprintf(stderr,"\t -n <number_of_events>\n");
   fprintf(stderr,"\t -e <event_type> = info at starting\n");
   fprintf(stderr,"\t\t SIG(104.973MeV/c), DIOA(0<p<105), DIO(watanabe_shanker over 101MeV/c), Proton(Silicon spectrum,thre=50MeV/c), RPC, TXT(-f <input_par.txt>)\n");
   fprintf(stderr,"\t -j <hit_type>\n");
   fprintf(stderr,"\t\t RKTrackRep, GeaneTrackRep2, ROOT(-i <input.root>), TXT (-x <input.txt>)\n");
   fprintf(stderr,"\t -t <track_type>\n");
   fprintf(stderr,"\t\t RKTrackRep, GeaneTrackRep2\n");
   fprintf(stderr,"\t -p <extrap_pos>\n");
   fprintf(stderr,"\t\t FirstHit, BeforeInwall\n");
   fprintf(stderr,"\t -r <hist_RPC.txt>\n");
   fprintf(stderr,"\t -l <hitpoint_type>\n");
   fprintf(stderr,"\t\t Layer, Wire\n");
   fprintf(stderr,"\t\t Foramt: energy(MeV/c) # (/MeV/Pion x 10^{-4})\n");
   fprintf(stderr,"\t -g <TGeometryFile.root>\n"); // moved from output
   fprintf(stderr,"\t -d <proton_data.root>\n");
   fprintf(stderr,"\t -w <wire_config.txt>\n");
   fprintf(stderr,"\t -f <fitting_ifirst>\n");
   fprintf(stderr,"[outputs]\n");
   fprintf(stderr,"\t -o <output.root>\n");
   fprintf(stderr,"\t -z <hit_pos.txt>\n");

   //strcpy(arg_event_type,"SIG"); // SIG,DIO,Proton,TXT
   //strcpy(arg_hit_type,  "GeaneTrackRep2"); // RKTrackRep,GeaneTrackRep2,ROOT
   //strcpy(arg_track_type,"GeaneTrackRep2"); // RKTrackRep,GeaneTrackRep2
   //strcpy(arg_extrap_pos,"FirstHit"); //FirtHit,BeforeInwall
}

int check_zpos_ordered(int iev, int size, double* posz, int ifirst, int ilast, std::vector<int>& bad_idx)
{

   //t->Branch("hits_z",g_hits_det.posz,"hits_z[nhits]/D");
   for (int i=0; i<size; i++) {
      fprintf(stderr,"check_zpos_ordered: iev %d i %d posz %lf (%lf)\n",iev,i, posz[i], g_hits_det.posz[i]);
   }
   // Assume size >=3

   //printf("check_zpos_orderd: size %d\n",size);
   double z1 = posz[ifirst];
   double z2 = posz[ifirst+1];
   double dz = z2-z1;
   if (dz>0) {
      for (int i=ifirst+2; i<ilast; i++) {

         double diff = posz[i] - posz[i-1];
      fprintf(stderr,"check_zpos_ordered: iev %d i %d posz[i] %lf posz[i-1] %lf diff %lf\n",iev,i, posz[i], posz[i-1],diff);
         if (diff<0)  {
            bad_idx.push_back(i);
            i++; // skip removed point
            return -1;
         }
      }
   } else {
      for (int i=ifirst+2; i<ilast; i++) {
         double diff = posz[i] - posz[i-1];
         if (diff>0) {
            bad_idx.push_back(i);
            i++;
            return -1;
         }
      }
   }
   return 0;
}

void read_hit_from_root_scan(int iev, char* txt, TVector3& posini, TVector3& momini)
{
   // positin and momentum is required!
   /*
    * Format of ROOT scan is as follows.
    *
    ***********************************************************
    *    Row   * Instance * hits_x_sm * hits_y_sm * hits_z_sm * momentum!!
    ***********************************************************
    *    12539 *        0 * -34.49890 * 42.530617 * 45.934633 *
    *
    *
    *
    ***********************************************************
    ==> 21 selected entries
    */
   printf("read_hit_from_root_scan: txt %s\n",txt);
   FILE* fp = fopen(txt,"r");
   if (fp==NULL) {
      fprintf(stderr,"read_hit_from_root_scan: failed to open %s\n",txt);
      exit(1);
   }
   char line[1280];
   int n=1;
   char a1[2];
   char a2[2];
   char a3[2];
   char a4[2];
   char a5[2];
   char a6[2];
   char a7[2];
   char a8[2];
   char a9[2];
   int row,ins;
   double xpos[MAX_HIT], ypos[MAX_HIT], zpos[MAX_HIT];
   double xmom[MAX_HIT], ymom[MAX_HIT], zmom[MAX_HIT];
   int nhits=0;
   while (fgets(line,sizeof(line),fp)) {
      if (n>3) {
         if (line[1]=='*') break;

         sscanf(line,"%s %d %s %d %s %lf %s %lf %s %lf %s %lf %s %lf %s %lf %s",a1,&row,a2,&ins,a3,
               &xpos[nhits],a4,&ypos[nhits],a5,&zpos[nhits],a6,
               &xmom[nhits],a7,&ymom[nhits],a8,&zmom[nhits],a9);
         //printf("nhits %d xpos %lf\n",nhits,xpos[nhits]);
         nhits++;
      }
      n++;
   }
   fclose(fp);

   // Only hit points are filled
   //     skip filling posini, momini

   for (int i=0; i<nhits; i++) {
      g_hits_det.posx[g_hits_det.nhits] = xpos[i];
      g_hits_det.posy[g_hits_det.nhits] = ypos[i];
      g_hits_det.posz[g_hits_det.nhits] = zpos[i];
      g_hits_det.momx[g_hits_det.nhits] = xmom[i];
      g_hits_det.momy[g_hits_det.nhits] = ymom[i];
      g_hits_det.momz[g_hits_det.nhits] = zmom[i];
      printf("%d %lf %lf %lf - %lf %lf %lf\n",i,xpos[i],ypos[i],zpos[i], xmom[i],ymom[i],zmom[i]);
      g_hits_det.nhits++;
   }
}
int set_smeared_hits(int iev, int ifirst, int ilast, struct tree_value& tv, GFAbsRecoHit** list_of_theHit)
{

   int nhits = ilast-ifirst+1;
   GFAbsRecoHit* theHit;
   printf("loop_over_ifrst: iev %d ifirst %d ilast %d nhits %d\n",iev,ifirst,ilast,nhits);

   //std::vector<double> pos_x;
   //std::vector<double> pos_y;
   //std::vector<double> pos_z;

//   print_plane(plane_extrap,"after_loop_over_ifirst: ");

   //
   // Add Hit point in the exit of target
   //
   if (strcmp(arg_extrap_pos,"TargetExit")==0) {
#if 0
      double xy_err_of_silicon = 0.0043; // 0.01cm = 0.1mm = 100um 43um
      double z_err_of_silicon  = 0.0072; // 0.01cm = 0.1mm = 100um 43um
      //double x_exit = gRandom->Gaus(root_exit_xpos_in_target,xy_err_of_silicon);
      //double y_exit = gRandom->Gaus(root_exit_ypos_in_target,xy_err_of_silicon);
      //double z_exit = gRandom->Gaus(root_exit_zpos_in_target,z_err_of_silicon);
      double x_exit = gRandom->Gaus(root_exit_xpos_in_target,0.0);
      double y_exit = gRandom->Gaus(root_exit_ypos_in_target,0.0);
      double z_exit = gRandom->Gaus(root_exit_zpos_in_target,0.0);
      double r_exit = sqrt(x_exit*x_exit+y_exit*y_exit);
      TVector3 point1_exit (x_exit-0.001,y_exit,z_exit);
      TVector3 point2_exit (x_exit-0.001,y_exit,z_exit+10.0);
      theHit = new WireHit(point1_exit, point2_exit, r_exit, z_exit, 0.01, 0.01);
      fitTrack.addHit(
            theHit,
            3,//dummy detector id
            0, // hit Id
            0, // rho
            0); // plane Id
      //fprintf(stderr,"Hit at TargetExit is added\n");



      // Add hits in dummy volume after hit point after target exit

      for (int i=0; i<g_hits_dum.nhits; i++) {
         double x_dum = gRandom->Gaus(g_hits_dum.posx[i],0.0);
         double y_dum = gRandom->Gaus(g_hits_dum.posy[i],0.0);
         double z_dum = gRandom->Gaus(g_hits_dum.posz[i],0.0);
         double r_dum = sqrt(x_dum*x_dum+y_dum*y_dum);
         TVector3 point1_dum (x_dum-0.001,y_dum,z_dum);
         TVector3 point2_dum (x_dum-0.001,y_dum,z_dum+10.0);
         theHit = new WireHit(point1_dum, point2_dum, r_dum, z_dum, 0.01, 0.01);
         fitTrack.addHit(
               theHit,
               3,//dummy detector id
               0, // hit Id
               0, // rho
               0); // plane Id
      }
#endif
#if 0
      //  <<=== It seems that mixture of wirehit and pinthit does not work... => check later
      TVector3 exit_point( root_exit_xpos_in_target, root_exit_ypos_in_target, root_exit_zpos_in_target);
      double xy_err_of_silicon = 0.0043; // 0.01cm = 0.1mm = 100um 43um
      double z_err_of_silicon  = 0.0072; // 0.01cm = 0.1mm = 100um 43um
      theHit = new PointHit(exit_point,xy_err,z_err);
      PointHit* ph_exit = dynamic_cast<PointHit*>(theHit);
      ph_exit->get_smeared_pos(
            g_hits_exit_posx_smeared,
            g_hits_exit_posy_smeared,
            g_hits_exit_posz_smeared);
      //
      fitTrack.addHit(
            theHit,
            3,//dummy detector id
            ilast+1, // hit Id
            ilast+1, // rho
            ilast+1); // plane Id
      //fprintf(stderr,"Hit at TargetExit is added\n");
#endif
   }

   double xy_err = 0.01;
   double z_err = 0.2;

   int num_of_theHit=0;
   for (int i=ifirst; i<ilast+1; i++) {

      TVector3 point( g_hits_det.posx[i], g_hits_det.posy[i], g_hits_det.posz[i]);

      if (strcmp(arg_hit_type,"TXT")==0) {

         if (strcmp(arg_hitpoint_type,"Wire")==0) {
            fprintf(stderr,"Only hitpoint_type=Layer is implemented\n");
            exit(1);
         }

         //theHit = new PointHit(point,posErr.X(),posErr.Z(), false);
         theHit = new PointHit(point,xy_err,z_err, false);
         g_hits_posx_smeared[i] = g_hits_det.posx[i];
         g_hits_posy_smeared[i] = g_hits_det.posy[i];
         g_hits_posz_smeared[i] = g_hits_det.posz[i];

      } else {
         //theHit = new PointHit(point,posErr.X(),posErr.Z());

         if (strcmp(arg_hitpoint_type,"Wirepoint")==0) {

            int ilayer = g_wire_ilayer[i];
            int icell  = g_wire_icell [i];
            double rdrift = g_wire_rdrift[i];    // closest approach 
            double zreco = g_wire_zreco [i];     // and z position
            //printf("theHit: Wire: i %d ilayer %d icell %d rdrift %lf zreco %lf\n",i,ilayer,icell,rdrift,zreco);

            double px, py, pz, vx, vy, vz;
            config_wire_pos_and_vector(g_config, ilayer, icell, &px, &py, &pz, &vx, &vy, &vz);
            double px2 = px+vx;
            double py2 = py+vy;
            double pz2 = pz+vz;
            //printf("px %lf py %lf pz %lf\n",px,py,pz);
            //printf("px2 %lf py2 %lf pz2 %lf\n",px2,py2,pz2);

            TVector3 point1(px,py,pz); // wire position at the downstream end plate
            TVector3 point2(px2, py2, pz2); // wire position at the upstream end plate

            //double rdrift_err = 0.001; // 100um
            //double rdrift_err = 0.1; // 100um
            double rdrift_err = 0.01; // 100um

         } else if (strcmp(arg_hitpoint_type,"Wire")==0) {

            int ilayer = g_wire_ilayer[i];
            int icell  = g_wire_icell [i];

            double px, py, pz, vx, vy, vz;
            config_wire_pos_and_vector(g_config, ilayer, icell, &px, &py, &pz, &vx, &vy, &vz);
            double px2 = px+vx;
            double py2 = py+vy;
            double pz2 = pz+vz;

            //TVector3 point1(px,py,pz); // wire position at the downstream end plate
            //TVector3 point2(px2, py2, pz2); // wire position at the upstream end plate
            
            // 2013/3/18 <= THIS IS!!!
//            TVector3 point1(px,py,pz-150./2.0);
//            TVector3 point2(px2, py2, pz2-150./2.0);
            // 2013/4/14 <= Again, fix this part...
            double wire_length = g_config->length;
            TVector3 point1(px,py,pz-wire_length/2.0);
            TVector3 point2(px2, py2, pz2-wire_length/2.0);

            //double rdrift_err = 0.01; // 100um
            //double rdrift_err = 0.11; // 100um
// comment-out (2013/5/20)            double rdrift_err = 0.02; // 200um 

            //double rdrift = g_wire_rdrift[i];    // closest approach 
            double rdrift_smeared = g_wire_rdrift_smeared[i];    // closest approach 
            //fprintf(stderr,"====>>>>>> rdrfit_smeard %d %lf\n", i, rdrift_smeared);

            theHit = new WireHit(point1, point2, rdrift_smeared, g_rdrift_err, false); // do not smeare inside 

            // comment-out (2013/6/4), moved to read_hit function
            //WireHit* wp = dynamic_cast<WireHit*>(theHit);
            //g_wire_rdrift_smeared[i] = wp->get_smeared_rdrift();
            //fprintf(stderr,"i %d g_wire_rdrift_smeared %lf\n", i, g_wire_rdrift_smeared[i]);

            // Set initial wire position for initial guess (2013/5/28)
            if (i==0) {
               g_fitini_1st_wire_posx = px;
               g_fitini_1st_wire_posy = py;
               g_fitini_1st_wire_rdrift = g_wire_rdrift_smeared[i];
            }

         } else {

               theHit = new PointHit(point,xy_err,z_err);
               PointHit* ph = dynamic_cast<PointHit*>(theHit);
               ph->get_smeared_pos(
                     g_hits_posx_smeared[i],
                     g_hits_posy_smeared[i],
                     g_hits_posz_smeared[i]);
            }
         }
         //if (iev==21) {
         //  printf("iev %d xy_err %lf z_err %lf\n",iev,xy_err,z_err);
         //   printf("iev %d posx %lf posx_sm %lf\n",iev,g_hits_det.posx[i],g_hits_posx_smeared[i]);
         //   printf("iev %d posy %lf posy_sm %lf\n",iev,g_hits_det.posy[i],g_hits_posy_smeared[i]);
         //   printf("iev %d posz %lf posz_sm %lf\n",iev,g_hits_det.posz[i],g_hits_posz_smeared[i]);
         //}

         //pos_x.push_back(g_hits_posx_smeared[i]);
         //pos_y.push_back(g_hits_posy_smeared[i]);
         //pos_z.push_back(g_hits_posz_smeared[i]);

//         printf("fitTrack.addHit i=%d\n",i);

      // DEBUG
//      if (i==12) continue;
//

          list_of_theHit[num_of_theHit++] = theHit;
          if (num_of_theHit >= MAX_THE_HIT) {
             return num_of_theHit;
          }
//         fitTrack.addHit(
//               theHit,
//               3,//dummy detector id
//               i, // hit Id
//               i, // rho
//               i); // plane Id

      }


#if 0
      // Before fitting, check that z_smeard is ordered.
      // By this, large fit_pa is rejected
      //printf("loop_over_ifirst: A iev %d\n",iev);
      std::vector<int> bad_idx;
      if (check_zpos_ordered(iev, nhits, g_hits_posz_smeared, ifirst, ilast, bad_idx)==-1) {
         tv.error = 444;
         // For checking, these events will be fitted. <=== tv.error will be overwritten! (2012/10/02)
         //return -1;
      }
      //printf("loop_over_ifirst: B iev %d\n",iev);
      //std::vector<int> bad_idx;
      //check_zpos_ordered(nhits,g_hits_posz_smeared, ifirst, ilast, bad_idx);
      // remove bad hit positions
#endif

      // this change hit positions (g_hits_det.xxx)
      //reject_bad_zpos();
      //void reject_bad_zpos()
      //{
      //   int nhits = g_hits_det.nhits;
      //   for (int i=0; i<nhits; i++) {
      //      g_hits_det.posz[i];
      //   }
      // do not use bad hit points check by above
      //std::vector< int >::iterator location;
      //location = std::find( bad_idx.begin(), bad_idx.end(), i );
      //// if find, skip this hit
      //if ( location != bad_idx.end() )  {
      //   printf("iev %d skipped due to woron zpos\n",iev);
      //   continue;
      //}
      //
      //
      //if (i2<3) {
      //   //printf("iev skip fitting %d\n",iev);
      //   printf("iev %d continued since nhits<3 after zpos removing\n",iev);
      //   return -2;
      //}


      //};



      /* Fitting will be performed when more than 3 hits in the detector */

      double tmp_px1 = g_hits_det.momx[0];
      double tmp_py1 = g_hits_det.momy[0];
      double tmp_pz1 = g_hits_det.momz[0];
      double tmp_pa1 = sqrt3(tmp_px1,tmp_py1,tmp_pz1);
      double tmp_px2 = g_hits_det.momx[g_nhits_fit-1];
      double tmp_py2 = g_hits_det.momy[g_nhits_fit-1];
      double tmp_pz2 = g_hits_det.momz[g_nhits_fit-1];
      double tmp_pa2 = sqrt3(tmp_px2,tmp_py2,tmp_pz2);
      diff_pa_max = tmp_pa1-tmp_pa2;

      f_hit_x = g_hits_det.posx[0];
      f_hit_y = g_hits_det.posy[0];
      f_hit_z = g_hits_det.posz[0];
      f_hit_px = g_hits_det.momx[0];
      f_hit_py = g_hits_det.momy[0];
      f_hit_pz = g_hits_det.momz[0];

      l_hit_x = g_hits_det.posx[g_nhits_fit-1];
      l_hit_y = g_hits_det.posy[g_nhits_fit-1];
      l_hit_z = g_hits_det.posz[g_nhits_fit-1];
      l_hit_px = g_hits_det.momx[g_nhits_fit-1];
      l_hit_py = g_hits_det.momy[g_nhits_fit-1];
      l_hit_pz = g_hits_det.momz[g_nhits_fit-1];
      //printf("l_hit_pos %lf %lf %lf\n",l_hit_x,l_hit_y,l_hit_z);
      //printf("l_hit_mom %lf %lf %lf\n",l_hit_px,l_hit_py,l_hit_pz);

      TVector3 pos_last_hit = TVector3(l_hit_x,l_hit_y,l_hit_z);
      TVector3 mom_last_hit = TVector3(l_hit_px,l_hit_py,l_hit_pz);
      GFDetPlane plane_last_hit = GFDetPlane(pos_last_hit,mom_last_hit.Unit());

      TVector3 pos_first_hit = TVector3(f_hit_x,f_hit_y,f_hit_z);
      TVector3 mom_first_hit = TVector3(f_hit_px,f_hit_py,f_hit_pz);
      GFDetPlane plane_first_hit = GFDetPlane(pos_first_hit,mom_first_hit.Unit());
      //printf("f_hit_pos %lf %lf %lf\n",f_hit_x,f_hit_y,f_hit_z);
      //printf("f_hit_mom %lf %lf %lf\n",f_hit_px,f_hit_py,f_hit_pz);



//      print_plane(plane_extrap,"before_fitting: ");
      //if (iev<4964) 
      //{
      //   continue;
      //}
      
      // plane_last_hit is not used in fitting..
      //double chi2 = fitting(iev,&tv,fitTrack,rep, plane_extrap, plane_last_hit, pos_extrap,pos_last_hit);

      /*
       *
double fitting(int iev, struct tree_value* tv, GFTrack& fitTrack, GFAbsTrackRep* rep, 
      GFDetPlane& plane,
      GFDetPlane& last_plane,
      TVector3& pos_ini,
      TVector3& pos_last)
      */

      //printf("loop_over_ifirst: C iev %d pos_extrap (%lf,%lf,%lf)\n",iev,pos_extrap.X(),pos_extrap.Y(),pos_extrap.Z());
      //printf("loop_over_ifirst: C iev %d chi2 %lf\n",iev,chi2);
      return num_of_theHit;
}

int read_smeared_hits()
{
   //
   // Call BEFORE set_semared_hits
   //
   /* overwrite 
    *  g_wire_ilayer
    *  g_wire_icell
    *  g_wire_rdrift_smeared
    *
    * */
   int wire_nhits;
   int wire_ilayer[1000];
   int wire_icell[1000];
   double wire_rdrift_smeared[1000];
   double wire_posx[1000];
   double wire_posy[1000];
   double wire_posz[1000];
   double wire_momx[1000];
   double wire_momy[1000];
   double wire_momz[1000];

   TFile* fin = new TFile("output_root_2/out_root652/hoge1.root");
   TTree*tin = (TTree*)fin->Get("t");
         printf("read_smeared_hits tin %p\n",tin);
   tin->SetBranchAddress("wire_nhits", &wire_nhits);
   tin->SetBranchAddress("wire_ilayer", wire_ilayer);
   tin->SetBranchAddress("wire_icell", wire_icell);
   tin->SetBranchAddress("wire_rdrift_smeared", wire_rdrift_smeared);
   tin->SetBranchAddress("wire_posx", wire_posx);
   tin->SetBranchAddress("wire_posy", wire_posy);
   tin->SetBranchAddress("wire_zreco", wire_posz);
   tin->SetBranchAddress("wire_momx", wire_momx);
   tin->SetBranchAddress("wire_momy", wire_momy);
   tin->SetBranchAddress("wire_momz", wire_momz);
   tin->GetEntry(29);
   printf("read_smeared_hits: wire_nhits %d\n", wire_nhits);
   for (int i=0; i<wire_nhits; i++) {
      g_wire_ilayer[i] = wire_ilayer[i];
      g_wire_icell [i] = wire_icell[i];
      g_wire_rdrift_smeared [i] = wire_rdrift_smeared[i];
      printf("read_smeared_hits: ihit %d ilayer %d icell %d rdrift_smeared %lf\n", 
            i, g_wire_ilayer[i], g_wire_icell[i], g_wire_rdrift_smeared[i]);
   }
   fin->Close();

   g_wire_ini_posx = wire_posx[0];
   g_wire_ini_posy = wire_posy[0];
   g_wire_ini_posz = wire_posz[0];
   g_wire_ini_momx = wire_momx[0];
   g_wire_ini_momy = wire_momy[0];
   g_wire_ini_momz = wire_momz[0];

}

int main(int argc, char** argv)
{
   if (argc==1) {
      print_usage(argv[0]);
      return -1;
   }
   init_args();

   int result;
   while((result=getopt(argc,argv,"s:x:i:n:e:j:t:p:f:c:o:g:z:m:r:l:d:w:h:k:"))!=-1){
      switch(result){
         /* INPUTS */
         case 's':
            arg_seed = atoi(optarg);
            gRandom->SetSeed(arg_seed);
            printf("gRandom->SetSeed(%d)\n",arg_seed);
            break;
         case 'x':
            strcpy(arg_input_txt,optarg);
            printf("input_txt %s\n",arg_input_txt);
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
         case 'r':
            g_hist_RPC = read_RPC_data(optarg);
            printf("RPC data %s\n",optarg);
            break;
         case 'l':
            strcpy(arg_hitpoint_type,optarg);
            printf("hitpoint_type %s\n",arg_hitpoint_type);
            break;
         case 'g':
            strcpy(arg_genfitGeom_name,optarg);
            printf("genfitGeom_name %s\n",arg_genfitGeom_name);
            break;
         case 'd': // added 2013/6/2
            strcpy(arg_proton_data,optarg);
            printf("proton_data %s\n",arg_proton_data);
            break;
         case 'w': // added 2013/6/2
            strcpy(arg_wire_config_fname,optarg);
            printf("wire_config_fname %s\n",arg_wire_config_fname);
            break;
         case 'k': // added 2013/6/17
            arg_fitting_ifirst = atoi(optarg);
            printf("fitting_ifirst %d\n",arg_fitting_ifirst);
            break;

            /* OUTPUTS */
         case 'o':
            strcpy(arg_output_root,optarg);
            printf("output_root %s\n",arg_output_root);
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
      if (strcmp(arg_proton_data,"empty")!=0) {
         read_proton_data(arg_proton_data);
      }

      // Wire config
      g_config = (struct config*)malloc(sizeof(struct config));

      if (strcmp(arg_wire_config_fname,"empty")==0) {
         fprintf(stderr,"ERROR: no wire configutation wire is set\n");
         exit(1);
      }


      // Read configure for wires
//      config_init("DriftChamber_new6/config.txt", g_config);
// B=1.0T
//      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/config.txt", g_config);
//      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/config_allaxial.txt", g_config);
//
// B=1.5T
      //config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/config_B1.5T.txt", g_config);
      //config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/config_B1.5T_new.txt", g_config);
//      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/config_allaxial_b15.txt", g_config);


      // New cofigureation (2013/2/22)

#if 0
      if (g_solenoid_bfld_tesla==1.0) {
      // B=1.0T
      printf("config_init 1.0T\n");
      //config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/config.txt", g_config);
      //config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/config_allaxial.txt", g_config);
      //config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/new_config.txt", g_config);
      //config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/all_stereo_config.txt", g_config); // 2013/3/13
      //config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/stereo_axial.txt", g_config); // 2013/3/13
//     config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/B1T_square_less.txt", g_config);
//      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/B1T_square_less_1m.txt", g_config); // 2013/4/11
//      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/B1T_square_less_1m_shift6.txt", g_config); // 2013/4/11
      //config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/B1T_square_less_1.2m.txt", g_config); // 2013/4/12
//      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/B1T_square_20130507.txt", g_config); // 2013/5/16


//      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/B1T_square_less_6446.txt", g_config); // 2013/4/17
//      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/B1T_square_SSSS.txt", g_config); // 2013/5/22
      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/B1T_square_SASA.txt", g_config); // 2013/5/24
//      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/B1T_square_All_axial.txt", g_config); // 2013/5/28
//      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/B1T_square_A6S4A9.txt", g_config); // 2013/5/29

      } else if (g_solenoid_bfld_tesla==0.7) {
         printf("config_init 0.7T\n");
         config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/B1T_square_Kuno.txt", g_config);
      } else if (g_solenoid_bfld_tesla==1.5) {
      // B=1.5T
      printf("config_init 1.5T\n");
      config_init("../../geant4_vmc/examples/B01/DriftChamber_new7/new_config_B1.5T.txt", g_config);
      } else {
         fprintf(stderr,"ERROR: @config_init g_solenoid_bfld_tesla is not neither 1.0 or 1.5.");
         exit(1);
      }
#endif
      config_init(arg_wire_config_fname, g_config);
      config_print(stderr, g_config);


      stMCT  = new TMatrixT<double>;
      covMCT = new TMatrixT<double>;
      stREC  = new TMatrixT<double>;
      covREC = new TMatrixT<double>;

      /* Output Root file */
      TFile* fout = new TFile(arg_output_root,"recreate");
      struct tree_value tv;
      struct tree_value prev_tv; // for previous fitting result
      //TTree* tout = make_branch("t",&min_tv);
      TTree* tout = make_branch("t",&tv);
      tout->SetMarkerStyle(2);

      fill_tree_config("t2");


      //construct_geom();
      TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
      //TGeoManager::Import(genfitGeom_name);
      //TGeoManager::Import("/work2/sakamoto/mu3e/Genfit/genfit/geant4_vmc/examples/C01/genfit_geom/geom_nowire.root");
      //TGeoManager::Import("/work2/sakamoto/mu3e/Genfit/genfit/geant4_vmc/examples/C01/genfit_geom/geom_nowire_spacing_3cm.root");
      //TGeoManager::Import("/work2/sakamoto/mu3e/Genfit/genfit/geant4_vmc/examples/C01/genfit_geom/geom_nowire_spacing_2cm.root");
      //TGeoManager::Import("/work2/sakamoto/mu3e/Genfit/genfit/geant4_vmc/examples/C01/genfit_geom/geom_nowire_spacing_4cm.root");

      TGeoManager::Import(arg_genfitGeom_name);
      //--- close the geometry
      gGeoManager->CloseGeometry();
      //exit(1);


      if (strcmp(g_solenoid_bfld_type,"uniform")==0) {
         GFFieldManager::getInstance()->init(new GFConstField(0.,0.,g_solenoid_bfld_tesla*10.0));// tesla -> kGaus
      } else {
         fprintf(stderr,"bfld_type %s is not implemented\n",g_solenoid_bfld_type);
         exit(1);
      }

      //gRandom->SetSeed(0);
      /* Initialize TGeant3 (need in use of GeaneTrackRep2) */
      gROOT->Macro("config/Geane.C");
      //gROOT->Macro("config/Geane_gaus_scat.C");

      //TVector3 posErr(0.01,0.01,0.2); //cm
      //TVector3 posErr(0.01,0.01,1.5); //cm

      //TVector3 posErr(0.001,0.001,0.001); //cm

      //      TVector3 posErr(0.00,0.00,0.0); //cm  ==>> chi2 should be infinity but genfit does not stop...
      //
      
      TVector3 posErr(0.01,0.01,0.2); //cm 2013/5/27
//      TVector3 posErr(0.1,0.1,0.1); //cm

//      TVector3 posErr(1,1,1); //cm
//      TVector3 posErr(0.02,0.02,0.4); //cm

      //TVector3 posErr(0.01,0.01,0.01); //cm
      //TVector3 posErr(1.01,1.01,1.01); //cm
      //TVector3 posErr(100.1,100.1,100.1); //cm

//      TVector3 momErr(0.001,0.001,0.001); //GeV/c = 1MeV/c
//      TVector3 momErr(0.002,0.002,0.004); //GeV/c

      TVector3 momErr(0.002,0.002,0.004); //GeV/c 2013/5/27

//      TVector3 momErr(0.02,0.02,0.04); //GeV/c <= comment-out 2013/5/30

      //TVector3 momErr(0.000001,0.000001,0.000001); //GeV/c = 1keV/c <<<===== same as above
      //
      //
      //
      
//      GFAbsTrackRep* rep;
//      if (strcmp(arg_track_type,"GeaneTrackRep2") == 0 ) { 
//         rep = new GeaneTrackRep2;
//      }
//      GFTrack *fitTrack = NULL;
//      GFTrack *fitTrack = new GFTrack(rep);
//      GFTrack *fitTrack = new GFTrack;

      printf("arg_total %d\n",arg_total);
       //  fprintf(stderr,"HOGE!!\n");
       //  exit(1);
     for (int iev=0; iev<arg_total; iev++) {
//for (int iev=2787; iev<arg_total; iev++) {
//        if (iev==2788) break;
         // init
         prev_tv.init(iev);

         fprintf(stderr,"iev %d\n",iev);

         tv.iev = iev;
         tv.error = -1;
         tv.nfail = -1;
         tv.ndf = 0;
         tv.prob = 0;
         tv.chi2 = 10e10;
         g_nhits_fit = 0;
         g_wire_nhits=0;
         g_hits_biw.nhits = 0;
         g_hits_det.nhits = 0;
         g_hits_abs.nhits = 0;
         g_hits_tgt.nhits = 0;
         g_nhits_tgt_before_chamber = 0.0;
         tv.ini_disk_number = -1;
         g_nhits_stereo=0;
         g_nhits_axial=0;
         g_proton_ov_ncells_uniq = 0;

         TVector3 posini,momini;

         // Determine starting position
         if       (strcmp(g_target_type,"test")==0)     { tv.ini_disk_number = get_init_pos_test(posini); } 
         else if  (strcmp(g_target_type,"disk")==0)     { tv.ini_disk_number = get_init_pos(posini); } 
         else  if (strcmp(g_target_type,"cylinder")==0) { tv.ini_disk_number = get_init_pos_cylinder(posini); } 
         else  if (strcmp(g_target_type,"cone")==0)     { tv.ini_disk_number = get_init_pos_cone(posini); } 
         else {
            fprintf(stderr,"unkonw target type %s\n",g_target_type);
            exit(1);
         }

         // Determine starting momentum
         if      (strcmp(arg_event_type,"test")==0)   { get_init_mom_test(momini); } 
         else if (strcmp(arg_event_type,"SIG")==0)    { get_init_mom(momini); } 
         else if (strcmp(arg_event_type,"DIO")==0)    { get_init_mom_DIO(momini); } 
         else if (strcmp(arg_event_type,"Proton")==0) { PDGcode = 2212; get_init_mom_Proton(momini); } 
         else if (strcmp(arg_event_type,"RPC")==0)    { PDGcode = 22/*gamma*/; get_init_mom_RPC(momini); }
         else if (strcmp(arg_event_type,"TXT")==0)    { 
            posini.SetX(arg_posini_x);
            posini.SetY(arg_posini_y);
            posini.SetZ(arg_posini_z);
            momini.SetX(arg_momini_x);
            momini.SetY(arg_momini_y);
            momini.SetZ(arg_momini_z);
         } else if (strcmp(arg_event_type,"DIOA")==0)    {  // added 2012/06/17 => hit points should be made in geant4_vmc ...
            if (dio_spec==NULL) {
               dio_spec = new DIO;
               if (dio_spec->Read("DIO/dio_spectrum.txt")==-1) {
                  exit(1);
               }
               dio_spec->Print();
               printf("dio_spectrum.txt is read\n");
               exit(1);
            }

            get_init_mom_DIOA(momini);

         } else {
            fprintf(stderr,"unkown event type %s\n",arg_event_type);
            exit(1);
         }


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

            if (read_hit(iev,arg_input_root, posini,momini)==-1)  {

               break; // root_hit returns -1 at the end of events
            }

         } else if (strcmp(arg_hit_type,"TXT")==0) {
            read_hit_from_root_scan(iev,arg_input_txt, posini,momini);
            // only one event
            arg_total = 1;
         } else {
            fprintf(stderr,"unknown hit type %s\n",arg_hit_type);
            exit(1);
         }

         
         tv.ini_disk_number = get_disk_number(posini);
         tv.ini.x = posini.X();
         tv.ini.y = posini.Y();
         tv.ini.z = posini.Z();
         tv.ini.px = momini.X();
         tv.ini.py = momini.Y();
         tv.ini.pz = momini.Z();
         

         g_fitini_theta1 = -1000.0;
         g_fitini_theta2 = -1000.0;
         g_fitini_dtheta = -1000.0;
         g_fitini_zpos1 = -1000.0;
         g_fitini_zpos2 = -1000.0;
         g_fitini_dzpos = -1000.0;

         g_fit_ix = -1;
         g_fit_iy = -1;
         g_fit_iz = -1;
         g_fit_ideg = -1;
         g_fit_min_ipattern=-1;

         
         //tv.ini.x += 1.0;
         //tv.ini.y += 1.0;
         //tv.ini.z += 1.0;

         //tv.ini.px =-60./1000.0;
         //tv.ini.py = 50./1000.0;
         //tv.ini.pz = 100.0/1000.0;
         
         
         //printf("iev %d tv.ini.pz %lf\n",iev,tv.ini.pz);


         /* Start tracking */
         PDGcode = 11; // e-
//         PDGcode = -11; // e+

         //printf("#### hogehoge arg_track_type %s\n",arg_track_type);
         // initialize
         g_diffx = 0.0;
         g_diffy = 0.0;
         if (strcmp(arg_track_type,"NoFit")==0) {
            tout->Fill();
            continue;
         }


         /*
          * Set DEBUG
          */
         //if (iev==2934) {
         //   setenv("DEBUG","1",1);
         //}



         //printf("qqqq: arg_hit_type %s\n",arg_hit_type);
         if (strcmp(arg_hit_type,"ROOT")==0 || strcmp(arg_hit_type,"TXT")==0) {
            g_nhits_fit = g_hits_det.nhits;
//               printf("hoge10 g_hits_det.nhits %d\n",g_hits_det.nhits);
            if (strcmp(arg_hitpoint_type,"Wire")==0 || strcmp(arg_hitpoint_type,"Wirepoint")==0) {
//               printf("hoge11\n");
               g_nhits_fit = g_wire_nhits;
            }
         } else {
            g_nhits_fit = get_nhits_for_fitting2();
         }
         //fprintf(stderr,"iev %d g_wire_nhits %d g_nhits_fit %d\n",iev, g_wire_nhits, g_nhits_fit);
         //exit(1);

         //printf("#### hogehoge2 g_wire_nhits %d g_hits_det.nhits %d\n",g_wire_nhits,g_hits_det.nhits);

         /* No need to fitting. Check comments in read_hit function. 2011/11/06 */
         if (g_hits_det.nhits==0) {
            tv.error = 111;
            tout->Fill();
            continue;
         }
         // DEBUG
         //if (iev==2) {
         //   g_nhits_fit  = 20;
         //}
         //if (iev==14) g_nhits_fit=65;
         //if (iev==12) g_nhits_fit=34;
         //if (iev==18673) g_nhits_fit=40;
         //if (iev==7928) g_nhits_fit=50;
         //if (iev==2987) g_nhits_fit=60;


         // Determine extrapolated position
         if (strcmp(arg_extrap_pos,"TargetExit")==0) {

            tv.hit.x = root_exit_xpos_in_target;
            tv.hit.y = root_exit_ypos_in_target;
            tv.hit.z = root_exit_zpos_in_target;
            tv.hit.px = root_exit_xmom_in_target;
            tv.hit.py = root_exit_ymom_in_target;
            tv.hit.pz = root_exit_zmom_in_target;
            g_edep_in_target = root_edep_in_target;

            printf("TagetExit x %lf\n",tv.hit.x);
            printf("TagetExit y %lf\n",tv.hit.y);
            printf("TagetExit z %lf\n",tv.hit.z);
            printf("TagetExit px %lf\n",tv.hit.px);
            printf("TagetExit py %lf\n",tv.hit.py);
            printf("TagetExit pz %lf\n",tv.hit.pz);
            printf("TagetExit edep %lf\n",g_edep_in_target);

         } else if (strcmp(arg_extrap_pos,"FirstHit")==0) {

           // /* This is not correct. (2013/3/14)
//            tv.hit.x = g_hits_det.posx[0];
//            tv.hit.y = g_hits_det.posy[0];
//            tv.hit.z = g_hits_det.posz[0];
//            tv.hit.px = g_hits_det.momx[0];
//            tv.hit.py = g_hits_det.momy[0];
//            tv.hit.pz = g_hits_det.momz[0];
            


            
            tv.hit.x = g_wire_ini_posx;
            tv.hit.y = g_wire_ini_posy;
            tv.hit.z = g_wire_ini_posz;
            tv.hit.px = g_wire_ini_momx;
            tv.hit.py = g_wire_ini_momy;
            tv.hit.pz = g_wire_ini_momz;




//            printf("FirstHit: x %lf\n",tv.hit.x);
//            printf("FirstHit: y %lf\n",tv.hit.y);
//            printf("FirstHit: z %lf\n",tv.hit.z);
//            printf("FirstHit: px %lf\n",tv.hit.px);
//            printf("FirstHit: py %lf\n",tv.hit.py);
//            printf("FirstHit: pz %lf\n",tv.hit.pz);
            //exit(1);

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
         /*
          * Check tv.ini_disk_number
          */

         /* 
          * First check of nhits
         */

         // (Second check is at after removing bad zpositions)
         int nhits = g_nhits_fit;
         //int nhits = 10;
         if (nhits<3) {
            //printf("iev skip fitting %d\n",iev);
            tv.error = 50;
            tout->Fill();
            printf("iev %d continued since nhits<3\n",iev);
            continue;
         }
         //int checking_iev = 14639;
         int checking_iev = -1;
         //int checking_iev = 18788;
         //int ifirst=0;
         ////ifirst=8;
         ////nhits = nhits - ifirst;
         //if (iev<checking_iev) {
         //   continue;
         //}
         //printf("iev %d nhits %d g_nhtis %d\n",iev,g_hits_det.nhits,g_nhits_fit);
         //if (iev==14639) {
         //   ifirst=10;
         //   nhits = 30;
         //   g_nhits_fit=nhits-ifirst+1;
         //}


         // Set Hit point after smearing (Drift circle)
         GFAbsRecoHit* list_of_theHit[MAX_THE_HIT];
         int ifirst = 0;
         int ilast = g_wire_nhits-1;
//         double xy_err = 0.0;
//         double z_err = 0.0;
//         double rdrift_err = 0.02; // 200um 
//         double rdrift_err = 0.0150; // 150um 
//         double rdrift_err = 0.0100; // 100um 
//         int num_of_theHit = set_smear_hits(iev, ifirst, ilast, xy_err, z_err, rdrift_err, tv, list_of_theHit);



         // DEBUG (2013/6/25)
         //read_smeared_hits();

         int num_of_theHit = set_smeared_hits(iev, ifirst, ilast, tv, list_of_theHit);


         /*
          *
          * Initail Value for Fitting
          *
          */

         //for (int ipass=0; ipass<6; ipass++) {
         //   g_pa_hit_min[ipass]  = 1e10;
         //   g_pa_hit_max[ipass]  = -10;
         //   g_pa_hit_rms[ipass] = 1234;
         //}


         //TVector3 pos_extrap = TVector3(tv.hit.x+1,tv.hit.y+1,tv.hit.z+1);
//         TVector3 pos_extrap = TVector3(tv.hit.x,tv.hit.y,tv.hit.z);
//         double z_cand = gRandom->Uniform(tv.hit.z, 5.0); // set z with accuracy of 5cm
//         double z_cand = gRandom->Uniform(tv.hit.z, 2.0); // set z with accuracy of 2cm
//         double z_cand = gRandom->Uniform(tv.hit.z, 0.0); // set z with accuracy of 0cm
//         double z_cand = gRandom->Gaus(tv.hit.z, 0.0); // set z with accuracy of 0cm
// <== MC (2013/5/9)         TVector3 pos_extrap = TVector3(tv.hit.x,tv.hit.y,z_cand);
//         printf("z_cand %lf\n", z_cand);
//
//          TVector3 pos_extrap = TVector3(-100, -100, -100);
         //TVector3 pos_extrap = TVector3(55,0,tv.hit.z);
         //TVector3 pos_extrap = TVector3(0., 0., tv.hit.z);
         
         //double init_x = 3.0*root_wire_posx[0] - 2.0*root_wire_posx[1];
         //double init_y = 3.0*root_wire_posy[0] - 2.0*root_wire_posy[1];
         //double init_z = 3.0*root_wire_zreco[0] - 2.0*root_wire_zreco[1];
         //TVector3 pos_extrap = TVector3(init_x, init_y, init_z);
         
          //printf("tv.hit.pz %lf\n",tv.hit.pz);
//         TVector3 mom_extrap = TVector3(0.,0.,1.);

         // Initialize
         g_inipar_x  = -1000.0;
         g_inipar_y  = -1000.0;
         g_inipar_z  = -1000.0;
         g_inipar_px = -1000.0;
         g_inipar_py = -1000.0;
         g_inipar_pz = -1000.0;
         g_inipar_rad_1st_hit = -1000.0;

         tv.fitini.x= -2000.0;
         tv.fitini.y= -2000.0;
         tv.fitini.z= -2000.0;
         tv.fitini.px= -2000.0;
         tv.fitini.py= -2000.0;
         tv.fitini.pz= -2000.0;

  
         if (g_wire_nhits<10) {
            fprintf(stderr,"g_wire_nhits %d is not enough for fitting (should be more than 10), so skip this event\n", g_wire_nhits);
            fprintf(stdout,"g_wire_nhits %d is not enough for fitting (should be more than 10), so skip this event\n", g_wire_nhits);
            tv.error = 60;
            tout->Fill();
            continue;
         }

//         int fitting_type = 5;
//         int fitting_type = 4;
//         int fitting_type = 3;
         int fitting_type = 2;
//         int fitting_type = 1;
         /*
          * Fitting Types:
          * 1) MC
          * 2) MC + smear
          * 3) IniPar
          * 4) IniPar->FitPar (2 times fitting)
          * 5) IniPar(ifirst=0) -> IniPar(ifirst=15) (2 times fitting)
          *    !! Note that g_wire_nhits will be overwritten at 2nd times
          */


         IniPar* par;
         if (fitting_type>=3) {

            par = new IniPar(g_config, 1.0, -1); // config, bfield, charge
            for (int i=0; i<g_wire_nhits; i++) {
               // BIG BUG!!! (2013/05/14) <= g_wire_nhits should be used with g_wire_XXXX
               //      int ilayer = root_wire_ilayer[i];
               //      int icell = root_wire_icell[i];
               //      double radius = root_wire_rdrift[i]; // radius is not used (2013/05/09)
               int ilayer = g_wire_ilayer[i];
               int icell = g_wire_icell[i];
               double radius = g_wire_rdrift[i]; // <= Is it OK not to use smeared r?
               double radius_smeared = g_wire_rdrift_smeared[i]; // replate with smeared R (2013/5/20)
               /*
                  fprintf(stderr,"g_wire_ilayer %d g_wire_icell %d g_wire_rdrift %lf (smeared-> %lf)\n", ilayer, icell, radius, radius_smeared);
                  */
               par->AddHit(ilayer, icell, radius_smeared);
            }
            //fprintf(stderr,"=0=0=0=0 iev %d\n", iev);
            bool verbose=true;
            int ret=0;
            // SetBeforeFitting goes to in par.Calc
            //ret = par.SetBeforeFitting(iev, verbose);
            //if (ret!=0) {
            //   fprintf(stderr,"main.cxx:: IniPar::Calc return %d, so skip this event\n", ret);
            //   fprintf(stdout,"main.cxx:: IniPar::Calc return %d, so skip this event\n", ret);
            //   tv.error = 60;
            //   tout->Fill();
            //   continue;
            //}
            ret = par->Calc(iev, verbose);
            //par.table.Print();
            //par.DrawFit();


#if 0
            if (ret!=0) {
               fprintf(stderr,"main.cxx:: IniPar::Calc return %d, so skip this event\n", ret);
               fprintf(stdout,"main.cxx:: IniPar::Calc return %d, so skip this event\n", ret);
               //tv.error = 70-ret; // 71 or 72
               tv.error = 70; // 2013/5/27
               tout->Fill();
               continue;
            }
#endif

            // outputs from IniPar
            g_fitini_theta1 = par->crs_theta[0];
            g_fitini_theta2 = par->crs_theta[1];
            g_fitini_dtheta = par->d_theta;
            g_fitini_zpos1 = par->crs_zpos[0];
            g_fitini_zpos2 = par->crs_zpos[1];
            g_fitini_dzpos = par->d_zpos;

            g_inipar_x  = par->extrap_x;
            g_inipar_y  = par->extrap_y;
            g_inipar_z  = par->extrap_z;
            g_inipar_px = par->extrap_px;
            g_inipar_py = par->extrap_py;
            g_inipar_pz = par->extrap_pz;

            g_inipar_rad_1st_hit = par->rad_1st_hit;


            // initial parameter for fitting
            //         tv.fitini.x = par.extrap_x; // cm
            //         tv.fitini.y = par.extrap_y; // cm
            //         tv.fitini.z = par.extrap_z; // after changing IniPar::CalBoarderZpos_(Rough->Precision), extrap_z is already geant coord
            //         tv.fitini.px = par.extrap_px; // GeV
            //         tv.fitini.py = par.extrap_py; // GeV
            //         tv.fitini.pz = par.extrap_pz; // GeV
            //
            fprintf(stderr,"END OF IniPar\n");
            fprintf(stderr,"main.cxx: IniPar: x %lf y %lf z %lf px %lf py %lf pz %lf\n", g_inipar_x, g_inipar_y, g_inipar_z, g_inipar_px, g_inipar_py, g_inipar_pz);
            fprintf(stderr,"main.cxx: MCTrue: x %lf y %lf z %lf px %lf py %lf pz %lf\n", tv.hit.x, tv.hit.y, tv.hit.z, tv.hit.px, tv.hit.py, tv.hit.pz);

         }

         //for (int ipass=0;ipass<6;ipass++) {
         //printf("AA ipass %d g_pa_hit_rms*1000 %lf\n", ipass,g_pa_hit_rms[ipass]*1000);
         //}
         //         for (int ideg=0; ideg<1; ideg++) {
         //         for (int ideg=0; ideg<5; ideg++) {

         int npat = 1;
         if (fitting_type==4||fitting_type==5) npat = 2;

         //         int fitting_ifirst = 0;
         //         int fitting_ifirst = 20;
         int fitting_ifirst = arg_fitting_ifirst;
         printf("fitting_ifirst %d\n", fitting_ifirst);
         int fitting_ilast = num_of_theHit-1;
         for (int ipattern=0; ipattern<npat; ipattern++) {

            if (fitting_type==5) {
               if (ipattern==0) {
                  fitting_ifirst = 0;
                  fitting_ilast = num_of_theHit-1;
               } else if (ipattern==1) {
                  fitting_ifirst = 15;
                  fitting_ilast = num_of_theHit-1;
               }
               // update #/hits
               g_wire_nhits = fitting_ilast - fitting_ifirst + 1;
            }

         // LOOP over ifirst
         double min_chi2 = 10e10;

         // scanning length = > 30*0.1 = 3 cm 
         double xpos_test_pattern[3]; 
         double ypos_test_pattern[3]; 
         double zpos_test_pattern[20]; 
         double xmom_test_pattern[3]; 
         double ymom_test_pattern[3]; 
         double zmom_test_pattern[20]; 

         double rad_test_pattern[20]; 
//         // pattern 1
//         for (int i=0; i<3; i++) {
//            xpos_test_pattern[i] = tv.fitini.x + (i-1)*0.5;
//            ypos_test_pattern[i] = tv.fitini.y + (i-1)*0.5;
//            zpos_test_pattern[i] = tv.fitini.z + (i-1)*0.5;
//         }
//
//         // pattern 2
//         for (int i=0; i<3; i++) {
//            xpos_test_pattern[i] = tv.fitini.x + (i-1)*0.5;
//            ypos_test_pattern[i] = tv.fitini.y + (i-1)*0.5;
//         }
//         for (int i=0; i<3; i++) {
//            zpos_test_pattern[i] = tv.fitini.z + (i-10)*0.1;
//         }

//         // pattern 3
//         for (int i=0; i<3; i++) {
//            xpos_test_pattern[i] = tv.hit.x + (i-1)*0.5;
//            ypos_test_pattern[i] = tv.hit.y + (i-1)*0.5;
//            zpos_test_pattern[i] = tv.hit.z + (i-10)*0.1;
//         }

//         // pattern 4 (MC only)
//         for (int i=0; i<1; i++) {
//            xpos_test_pattern[i] = tv.hit.x;
//            ypos_test_pattern[i] = tv.hit.y;
//            zpos_test_pattern[i] = tv.hit.z;
//         }
//

// pattern 5 (IniPar only)
         for (int i=0; i<1; i++) {
            xpos_test_pattern[i] = g_inipar_x;
            ypos_test_pattern[i] = g_inipar_y;
            zpos_test_pattern[i] = g_inipar_z;
            xmom_test_pattern[i] = g_inipar_px;
            ymom_test_pattern[i] = g_inipar_py;
            zmom_test_pattern[i] = g_inipar_pz;
         }

         if (ipattern==1) {
            // extrapolate to the fitting_ifirst
            //par.Extrap_Pos(iev, fitting_ifirst);
         }


         // +/-40deg @10deg
         for (int i=0; i<5; i++) {
            rad_test_pattern[i] = g_inipar_rad_1st_hit + (i-2)*10./180.*TMath::Pi();
         }

         //         // pattern 6 (Onlyl x, y is smread)
         //         for (int i=0; i<3; i++) {
         //            xpos_test_pattern[i] = g_inipar_x + (i-1)*0.5;
         //            ypos_test_pattern[i] = g_inipar_y + (i-1)*0.5;
         //         }
         //         for (int i=0; i<1; i++) {
         //            zpos_test_pattern[i] = g_inipar_z;
         //            xmom_test_pattern[i] = g_inipar_px;
         //            ymom_test_pattern[i] = g_inipar_py;
         //            zmom_test_pattern[i] = g_inipar_pz;
         //         }

         // pattern 7 (Onlyl x, y is smread)
         //for (int i=0; i<3; i++) {
         //   xpos_test_pattern[i] = g_inipar_x + (i-1)*0.05;
         //   ypos_test_pattern[i] = g_inipar_y + (i-1)*0.05;
         //}
         //for (int i=0; i<1; i++) {
         //   zpos_test_pattern[i] = g_inipar_z;
         //   xmom_test_pattern[i] = g_inipar_px;
         //   ymom_test_pattern[i] = g_inipar_py;
         //   zmom_test_pattern[i] = g_inipar_pz;
         // }

         double new_fitini_x;
         double new_fitini_y;
         double new_fitini_z;
         double new_fitini_px;
         double new_fitini_py;
         double new_fitini_pz;

         TVector3 pos_ini;
         TVector3 mom_ini;
         TVector3 pos_extrap;
         TVector3 mom_extrap;




            //         for (int ix=0; ix<3; ix++) {
            //         for (int iy=0; iy<3; iy++) {
            //         for (int iz=0; iz<1; iz++) {
            //int ipattern=-1;
            //ipattern++;

            // <== IniPar
            // double new_fitini_x;
            // double new_fitini_y;
            // double new_fitini_z = zpos_test_pattern[0];
            // get_position_on_circle(g_fitini_1st_wire_posx, g_fitini_1st_wire_posy, g_fitini_1st_wire_rdrift, rad_test_pattern[ideg], new_fitini_x, new_fitini_y);

             if (fitting_type==1) {

               //<== MC only
               
               pos_extrap = TVector3(tv.hit.x,tv.hit.y,tv.hit.z);
               mom_extrap = TVector3(tv.hit.px,tv.hit.py,tv.hit.pz);
               pos_ini = TVector3(tv.hit.x,tv.hit.y,tv.hit.z);
               mom_ini = TVector3(tv.hit.px,tv.hit.py,tv.hit.pz);

            } else if (fitting_type==2) {

               // <= MC + Smear

               //double sm_mc_xpos = gRandom->Gaus(tv.hit.x, 0.5); // 0.5cm
               //double sm_mc_ypos = gRandom->Gaus(tv.hit.y, 0.5); // 0.5cm
               //double sm_mc_zpos = gRandom->Gaus(tv.hit.z, 1.5); // 1.5cm

               double sm_mc_xpos = gRandom->Gaus(tv.hit.x, 0.1); // 1mm
               double sm_mc_ypos = gRandom->Gaus(tv.hit.y, 0.1); // 1mm
               double sm_mc_zpos = gRandom->Gaus(tv.hit.z, 0.2); // 2mm
               //double sm_mc_zpos = gRandom->Gaus(tv.hit.z, 0.1); // 2mm
               double sm_mc_xmom = gRandom->Gaus(tv.hit.px, 0.002); // 2MeV/c
               double sm_mc_ymom = gRandom->Gaus(tv.hit.py, 0.002); // 2MeV/c
               double sm_mc_zmom = gRandom->Gaus(tv.hit.pz, 0.004); // 4MeV/c
               //double sm_mc_xpos = tv.hit.x;
               //double sm_mc_ypos = tv.hit.y;
               //double sm_mc_zpos = tv.hit.z;
               //double sm_mc_xmom = tv.hit.px;
               //double sm_mc_ymom = tv.hit.py;
               //double sm_mc_zmom = tv.hit.pz;
               //double sm_mc_xmom = g_inipar_px;
               //double sm_mc_ymom = g_inipar_py;
               //double sm_mc_zmom = g_inipar_pz;

               // extrapolate to wire position at 1st hit
               pos_extrap = TVector3(sm_mc_xpos, sm_mc_ypos, sm_mc_zpos);
               mom_extrap = TVector3(sm_mc_xmom, sm_mc_ymom, sm_mc_zmom); 
               //mom_extrap = TVector3(sm_mc_xpos, sm_mc_ypos, 0.0); //<= no change... (run517)

               pos_ini = TVector3(sm_mc_xpos, sm_mc_ypos, sm_mc_zpos);
               mom_ini = TVector3(sm_mc_xmom, sm_mc_ymom, sm_mc_zmom);
               //TVector3 mom_ini = TVector3(sm_mc_xpos, sm_mc_ypos, sm_mc_zpos);

               // initial parameter for fitting
               tv.fitini.x = sm_mc_xpos;
               tv.fitini.y = sm_mc_ypos;
               tv.fitini.z = sm_mc_zpos;
               tv.fitini.px = sm_mc_xmom;
               tv.fitini.py = sm_mc_ymom;
               tv.fitini.pz = sm_mc_zmom;

            } else if (fitting_type==3||fitting_type==4 || fitting_type==5) {

               if (ipattern==0) {

                  new_fitini_x = xpos_test_pattern[0];
                  new_fitini_y = ypos_test_pattern[0];
                  new_fitini_z = zpos_test_pattern[0];

                  //double sm_mc_zpos = gRandom->Gaus(tv.hit.z, 0.2); // 2mm
                  //new_fitini_z = sm_mc_zpos;

                  ////            double new_fitini_x = tv.hit.x;
                  ////            double new_fitini_y = tv.hit.y;
                  ////            double new_fitini_z = tv.hit.z;

                  //            double new_fitini_px = xmom_test_pattern[0];
                  //            double new_fitini_py = ymom_test_pattern[0];
                  //            double new_fitini_pz = zmom_test_pattern[0];

                  // New extrapolate plane! (2013/5/28)
                  new_fitini_px = xmom_test_pattern[0];
                  new_fitini_py = ymom_test_pattern[0];
                  //                double new_fitini_pz = zmom_test_pattern[0];
//                  new_fitini_pz = zmom_test_pattern[0]; //<= #/ev get lower
                  new_fitini_pz = 0;
                  //            TVector3 dir_extrap = TVector3(xmom_test_pattern[0], ymom_test_pattern[0], 0);

               } else if (ipattern>=1) {

                  // initial value for next fitting
                  new_fitini_x = tv.fit.x;
                  new_fitini_y = tv.fit.y;
                  new_fitini_z = tv.fit.z;
                  new_fitini_px = tv.fit.px;
                  new_fitini_py = tv.fit.py;
                  new_fitini_pz = 0;
//                  new_fitini_pz = tv.fit.pz;

                  // Correctly, set initial value at ifirst=15 for the 2nd fitting
                  // But, set ifirst=0 at the moment...
                  if (fitting_type==5) {
                     new_fitini_x = xpos_test_pattern[0];
                     new_fitini_y = ypos_test_pattern[0];
                     new_fitini_z = zpos_test_pattern[0];
                     new_fitini_px = xmom_test_pattern[0];
                     new_fitini_py = ymom_test_pattern[0];
                     new_fitini_pz = 0;
                  }


               }

               // Change the extrapolation plane so that surface is vertical to the axis made with 1st hit point and the center (0,0) (2013/5/28)
               pos_ini = TVector3(new_fitini_x,new_fitini_y, new_fitini_z);
               mom_ini = TVector3(new_fitini_px,new_fitini_py, new_fitini_pz);
               tv.fitini.x = new_fitini_x;
               tv.fitini.y = new_fitini_y;
               tv.fitini.z = new_fitini_z;
               tv.fitini.px = new_fitini_px;
               tv.fitini.py = new_fitini_py;
               tv.fitini.pz = new_fitini_pz;

               pos_extrap = TVector3(new_fitini_x,new_fitini_y, new_fitini_z);
               mom_extrap = TVector3(new_fitini_px,new_fitini_py, new_fitini_pz);
               printf("iev %d ipattern %d FitPar (x,y,z)    %lf %lf %lf (MC %lf %lf %lf)\n",    iev, ipattern, new_fitini_x, new_fitini_y, new_fitini_z, tv.hit.x,tv.hit.y,tv.hit.z);
               printf("iev %d ipattern %d FitPar (px,py,pz) %lf %lf %lf pt %lf (MC %lf %lf %lf) pt %lf\n", 
                     iev, ipattern, new_fitini_px, new_fitini_py, new_fitini_pz, sqrt2(new_fitini_px,new_fitini_py), 
                     tv.hit.px,tv.hit.py,tv.hit.pz, sqrt2(tv.hit.px, tv.hit.py));

            }


            //GFDetPlane plane_extrap = GFDetPlane(pos_extrap2,mom_extrap.Unit());
            GFDetPlane plane_ini;
            GFDetPlane plane_extrap;
            if (strcmp(arg_hitpoint_type,"Wire")==0) {

               //TVector3 mom_extrap_new(0., -1., 0.0);
               //plane_extrap = GFDetPlane(pos_extrap,mom_extrap_new.Unit());
               //plane_extrap = GFDetPlane(pos_extrap,mom_extrap.Unit());
               //TVector3 pos_extrap = TVector3(55,0,0);
               //         TVector3 mom_extrap = TVector3(1,1,1);
               //         TVector3 mom_extrap = TVector3(1,0,0);
               //         TVector3 mom_extrap = TVector3(0,0,1);
               plane_extrap = GFDetPlane(pos_extrap,mom_extrap.Unit());
               //plane_extrap = GFDetPlane(pos_extrap,dir_extrap.Unit());
               //TVector3 u(0,1,0);
               //TVector3 v(1,0,0);
               //plane_extrap = GFDetPlane(pos_extrap,u,v);
               //            printf("I'm Wire! (pos_extrap x,y,z=%lf,%lf,%lf)\n",pos_extrap.X(),pos_extrap.Y(),pos_extrap.Z());

               plane_ini = GFDetPlane(pos_ini,mom_ini.Unit());
            } else if (strcmp(arg_hitpoint_type,"Wirepoint")==0) {

               TVector3 mom_extrap_new(0., 1., 0.0);
               plane_extrap = GFDetPlane(pos_extrap,mom_extrap_new.Unit());

            } else {

               /* plne is perpendicular to the direction of mom_extrap */
               plane_extrap = GFDetPlane(pos_extrap,mom_extrap.Unit());
            }

            // Determine track follower
            GFAbsTrackRep* rep;
            if      (strcmp(arg_track_type,"RKTrackRep") == 0 )     { rep = new RKTrackRep( pos_extrap, mom_extrap, PDGcode); } 
            else if (strcmp(arg_track_type,"GeaneTrackRep2") == 0 ) { 

               rep = new GeaneTrackRep2( plane_ini, mom_ini, posErr, momErr, PDGcode); 
               //               rep = new GeaneTrackRep2( plane_extrap, mom_extrap, posErr, momErr, PDGcode); 
               //rep -> SetPlaneMom( plane_extrap, mom_extrap, posErr, momErr, PDGcode); 


            }  else {
               fprintf(stderr,"unknown track type %s\n",arg_track_type);
               exit(1);
            }
            // Set hit positions

            //            if (fitTrack!=NULL) {
            //               delete fitTrack;
            //            }
            //GFTrack fitTrack(rep);
            GFTrack* fitTrack = new GFTrack(rep);
            // TODO: understand why following clode does not work!
            //            fitTrack->reset();
            //            fitTrack->SetRep(rep);

            printf("=~=====>>>>>> num_of_theHit %d\n", num_of_theHit);

            /*
             * Tail Cut Test (2013/5/31)
             *
             * Star from idx=10
             * Compare chi2_hitX_max
             *
             */

            for (int i=fitting_ifirst; i<=fitting_ilast; i++) {
               fitTrack->addHit(
                     list_of_theHit[i],
                     3,//dummy detector id
                     i, // hit Id
                     i, // rho
                     i); // plane Id

            }
            //GFTrack *fitTrack = new GFTrack(rep);

            //         double chi2 = loop_over_ifirst(iev, ifirst, ilast, posErr.X(), posErr.Z(), fitTrack, rep,  plane_extrap, pos_extrap, tv);
            //         double chi2 = loop_over_ifirst(iev, ifirst, ilast, posErr.X(), posErr.Z(), fitTrack, rep,  plane_extrap2, pos_extrap2, tv);

            double chi2 = do_the_fitting(iev,&tv,fitTrack,rep, plane_extrap);
            if (fitting_type==4) {
               if (tv.error!=0) {
                  goto fill;
               }
            }
            if (fitting_type==4||fitting_type==5) {
               prev_tv.copy(&tv);
            }

            printf("#### iev %d ipattern %d chi2 %lf\n", iev, ipattern, chi2);
#if 0
            // using fit_z, fitiing again
            pos_extrap.SetZ(tv.fit.z);
            plane_extrap = GFDetPlane(pos_extrap,mom_extrap.Unit());
            if      (strcmp(arg_track_type,"RKTrackRep") == 0 )     { rep = new RKTrackRep( pos_extrap, mom_extrap, PDGcode); } 
            else if (strcmp(arg_track_type,"GeaneTrackRep2") == 0 ) { rep = new GeaneTrackRep2( plane_extrap, mom_extrap, posErr, momErr, PDGcode); }  
            else {
               fprintf(stderr,"unknown track type %s\n",arg_track_type);
               exit(1);
            }
            // Set hit positions
            GFTrack fitTrack2(rep);
            //GFTrack *fitTrack = new GFTrack(rep);
            chi2 = loop_over_ifirst(iev, ifirst, ilast, posErr.X(), posErr.Z(), fitTrack2, rep,  plane_extrap, pos_extrap, tv);
#endif

#if 0
            if (chi2==-1) {
               //printf("nhits %d\n",nhits);
               //for (int i=0; i<nhits; i++) {
               //   printf("===iev %d i %d posy %lf posy_sm %lf\n",iev,i,g_hits_det.posy[i],g_hits_posy_smeared[i]);
               //}
               tv.error = 200;
               tout->Fill();
               goto end;
            } else if (chi2==-2) {
               tv.error = 55;
               tout->Fill();
               goto end;
            }
#endif

            printf("## iev %d ipattern %d chi2 %lf pa_hit1_rms %lf (MeV/c) \n",iev,ipattern,chi2, g_pa_hit_rms[0]*1000);

            //if (chi2<min_chi2) {

            //               tv.fitini.x = new_fitini_x;
            //               tv.fitini.y = new_fitini_y;
            //               tv.fitini.z = new_fitini_z;
            //               tv.fitini.px = new_fitini_px;
            //               tv.fitini.py = new_fitini_py;
            //               tv.fitini.pz = new_fitini_pz;
            //

            if (ipattern==0) {
               //min_tv.copy(&tv);
               //min_chi2 = chi2;
               //min_tv.chi2 = min_chi2;
               //min_tv.ifirst = ifirst;
               //min_tv.ilast = ilast;

               //               g_fit_ix = ix;
               //               g_fit_iy = iy;
               //               g_fit_iz = iz;
               //               g_fit_ideg = ideg;

               g_fit_min_ipattern = ipattern;
            } else if (ipattern==1) {

               // compare previous result
               // If g_diff_fitpa is large, then this could make tail!!
               double fitpa1 = sqrt3(prev_tv.fit.px, prev_tv.fit.py, prev_tv.fit.pz); //previous result
               double fitpa2 = sqrt3(tv.fit.px, tv.fit.py, tv.fit.pz); //current result
               g_diff_fitpa = fitpa1 - fitpa2;
            }
            // always last fitting result
            //min_tv.copy(&tv);
            //min_chi2 = chi2;
            //min_tv.chi2 = min_chi2;
            //min_tv.ifirst = fitting_ifirst;
            //min_tv.ilast = fitting_ilast;


            //fprintf(stderr,"==== chi2<min_chi2 iev =%d\n", iev);
            //fprintf(stderr,"IniPar: initial (x,y,z)    %10.5f %10.5f %10.5f\n",min_tv.fitini.x,min_tv.fitini.y,min_tv.fitini.z);
            //fprintf(stderr,"IniPar: initial (px,py,pz) %10.5f %10.5f %10.5f\n",
            //      min_tv.fitini.px,min_tv.fitini.py,min_tv.fitini.pz);
            //} // chi2

            //if (min_chi2<5.0) {

            //double fitx = tv.fit.x;
            //double fity = tv.fit.y;
            //double fitz = tv.fit.z;
            //double fitinix = tv.fitini.x; // MC infor cannot be used in real
            //double fitiniy = tv.fitini.y;
            //double fitiniz = tv.fitini.z;
            //double res_x = fitx - fitinix;
            //double res_y = fity - fitiniy;
            //double res_z = fitz - fitiniz;
            //
            //if (TMath::Abs(res_x)<2 && TMath::Abs(res_y)<2 && TMath::Abs(res_z)<10) {
            //   printf("iev %d ipattern %d zpos_test %lf fitx %lf fitinix %lf res_x %lf\n", iev, ipattern, zpos_test_pattern[ipattern], fitx, fitinix, res_x);
            //   printf("iev %d ipattern %d zpos_test %lf fity %lf fitiniy %lf res_y %lf\n", iev, ipattern, zpos_test_pattern[ipattern], fity, fitiniy, res_y);
            //   printf("iev %d ipattern %d zpos_test %lf fitz %lf fitiniz %lf res_z %lf\n", iev, ipattern, zpos_test_pattern[ipattern], fitz, fitiniz, res_z);
            //   break; // ok, enough, let's break from ipattern
            //}
            //  break;
            // }

            //         } //iz
            //         } //iy
            //         } //ix
         } //ipattern
         //         } //ideg
         if (fitting_type==5) {
            tv.copy(&prev_tv); // set first fitting result
         }
         



#if 0 
         // smear is changed...

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
         //GFTrack *fitTrack = new GFTrack(rep);
         tv.ifirst = min_ifirst;
         tv.ilast = min_ilast;
         min_chi2 = loop_over_ifirst(iev, min_ifirst, min_ilast, posErr.X(), posErr.Z(), fitTrack, rep, plane_extrap, pos_extrap, tv);
         //printf("iev %d min_ifirst %d min_ilast %d chi2 %lf\n",iev, min_ifirst,min_ilast, min_chi2);
#endif



         /*
          *
          */
         if (getenv("THIS_IS_BAD_EVENT")!=NULL) {
            tv.error=555;
            unsetenv("THIS_IS_BAD_EVENT");
            fprintf(stderr,"tv.error=555\n");
         }

fill:
//         print_tree_value("", stdout, &min_tv);
         print_tree_value("", stdout, &tv); // tv is filled in tree, not min_tv
         tout->Fill();
end:
         //delete rep;

         //delete theHit;
         //delete fitTrack;
         if (iev==checking_iev) 
         {
            printf("g_nhits_fit %d\n",g_nhits_fit);
            break;
         }
         //for (int ipass=0;ipass<6;ipass++) {
         //printf("ipass %d g_pa_hit_rms*1000 %lf\n", ipass,g_pa_hit_rms[ipass]*1000);
         //}

next_loop:
         ;
     }
     printf("end\n");

     fout->cd();
     tout->Write();
     fout->Close();
     if (arg_writeout_track_fp) {
        fclose(arg_writeout_track_fp);
     }

     return 0;
     }
