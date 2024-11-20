#ifndef __CA_TREE_DIS_HH__
#define __CA_TREE_DIS_HH__


Int_t nphe; 
Int_t mult_e;
Int_t mult_pip;
Int_t mult_pim;
Int_t mult_g;
Int_t tgt_no;
TLorentzVector beam;
TLorentzVector target;
TLorentzVector Pe;
TVector3 vtx_e;
Float_t theta_e;
Float_t phi_e;
Double_t e_dE_pcal;
Double_t e_dE_ecin;
Double_t e_dE_ecout;
Double_t e_x_pcal;
Double_t e_x_ecin;
Double_t e_x_ecout;
Double_t e_y_pcal;
Double_t e_y_ecin;
Double_t e_y_ecout;
Double_t e_dE_ftof1;
Double_t e_dE_ftof2;
Double_t e_dE_ftof3;
Double_t SF;
Double_t e_Chi2DoF;
Double_t runCharge;
Double_t ft_x;
Double_t ft_y;
Double_t ft_r;
Int_t e_sec_dc;
double e_dc_x[3];
double e_dc_y[3];
double e_dc_z[3];
double e_htcc_x;
double e_htcc_y;
double e_htcc_z;
TVector3 e_edge_dc;
Float_t e_nph_htcc;
Float_t e_nph_ltcc;
TLorentzVector q;
Double_t theta_q;
Bool_t fiducial_cut;
Double_t Q2;
Double_t xB;
Double_t W2;
Double_t Yscale;
Int_t sector_e;
Int_t detector_e;
Double_t e_V;
Double_t e_W;
bool isGoodElectron;
/*
Int_t mult_p_cd;
Int_t mult_p_fd;
Int_t mult_p;
std::vector<Bool_t> p_fd;
std::vector<Bool_t> p_cd;
std::vector<Double_t> tcorr;
std::vector<Double_t> pbeta;
std::vector<Double_t> ptof;
std::vector<Double_t> ptofp;
std::vector<Double_t> ppath;
std::vector<Double_t> p_dE_cnd1;
std::vector<Double_t> p_dE_cnd2;
std::vector<Double_t> p_dE_cnd3;
std::vector<Double_t> p_dE_ctof;
std::vector<Double_t> p_tof_cnd1;
std::vector<Double_t> p_tof_cnd2;
std::vector<Double_t> p_tof_cnd3;
std::vector<Double_t> p_tof_ctof;
std::vector<Double_t> p_path_cnd1;
std::vector<Double_t> p_path_cnd2;
std::vector<Double_t> p_path_cnd3;
std::vector<Double_t> p_path_ctof;
std::vector<TLorentzVector> PpV;
Double_t PpX;
Double_t PpY;
Double_t PpZ;
Double_t PpE;
std::vector<std::vector<Double_t>> Pp;
std::vector<std::vector<Double_t>> Pp_true;
std::vector<double> pp;
std::vector<double> pp_true;
std::vector<std::vector<double> > PN_true;
Double_t PleadX;
Double_t PleadY;
Double_t PleadZ;
Double_t PleadE;
std::vector<std::vector<Double_t>> Plead;
std::vector<Double_t> theta_p;
std::vector<Double_t> phi_p;
std::vector<Int_t> sector_p_cd;
std::vector<Int_t> detector;
TVector3 vtx_p[10];
Double_t vtx_p_x;
Double_t vtx_p_y;
Double_t vtx_p_z;
std::vector<std::vector<Double_t>> vtxp;
std::vector<Bool_t> plead;
Int_t plead_index;
Int_t plead_index_reco;
Int_t plead_no;
Int_t plead_no_cd;
Int_t plead_no_fd;
std::vector<Bool_t> p_flag;
Int_t p_no;
Int_t p_no_cd;
Int_t p_no_fd;
*/

//TVector3 pmiss[10];
/*Double_t pmissP;
double pmissX;
double pmissY;
double pmissZ;
std::vector<std::vector<Double_t>> Pmiss;
std::vector<double> pmiss;
std::vector<double> Mmiss;
std::vector<double> Mmiss_fd;
std::vector<double> Mmiss_cd;
std::vector<double> Mmiss2;
std::vector<double> Emiss;
std::vector<double> Emiss0;
std::vector<double> alpha_miss;
std::vector<Double_t> pq;
std::vector<Double_t> theta_pq;
std::vector<Double_t> theta_pmq;
std::vector<Double_t> pmissq_para;
std::vector<Double_t> pmissq_perp;
*/

Double_t weight;
TLorentzVector Pe_true;
Float_t theta_e_true;
Float_t phi_e_true;
TLorentzVector q_true;
Double_t theta_q_true;
Double_t Q2_true;
Double_t xB_true;
Double_t W2_true;
/*Int_t mult_p_true;
TLorentzVector Plead_true;
//TVector3 pmiss_true;
//Double_t theta_p_true;
//Double_t phi_p_true;
Int_t plead_index_true;
Int_t plead_no_true;
double pmissP_true;
//double Mmiss_true;
//double Mmiss2_true;
//double Emiss_true;
//Double_t pq_true;
//Double_t theta_pq_true;
*/
/*
std::vector<std::vector<Double_t>> Pmiss_true;
std::vector<double> pmiss_true;
std::vector<double> Mmiss_true;
std::vector<double> Emiss_true;
std::vector<Double_t> pq_true;
std::vector<Double_t> theta_pq_true;
std::vector<Double_t> theta_p_true;
std::vector<Double_t> phi_p_true;

std::vector<Double_t> pim_dE_cnd1;
std::vector<Double_t> pim_dE_cnd2;
std::vector<Double_t> pim_dE_cnd3;
std::vector<Double_t> pim_dE_ctof;
std::vector<Double_t> pim_p;
*/

tout.Branch("tgt_no", &tgt_no);
tout.Branch("Pbeam", &beam);
tout.Branch("nphe",&nphe);
//tout.Branch("mult_e", &mult_e);
tout.Branch("Pe", &Pe);
tout.Branch("theta_e", &theta_e);
tout.Branch("phi_e", &phi_e);
tout.Branch("sector_e", &sector_e);
tout.Branch("detector_e", &detector_e);
tout.Branch("SF_e", &SF);
tout.Branch("e_dE_ecout",&e_dE_ecout);
tout.Branch("e_dE_pcal",&e_dE_pcal);
tout.Branch("e_dE_ecin",&e_dE_ecin);
tout.Branch("e_x_ecout",&e_x_ecout);
tout.Branch("e_x_pcal",&e_x_pcal);
tout.Branch("e_x_ecin",&e_x_ecin);
tout.Branch("e_y_ecout",&e_y_ecout);
tout.Branch("e_y_pcal",&e_y_pcal);
tout.Branch("e_y_ecin",&e_y_ecin);
tout.Branch("vtx_e", &vtx_e);
tout.Branch("e_dc_x", &e_dc_x,"e_dc_x[3]/d");
tout.Branch("e_dc_y", &e_dc_y,"e_dc_y[3]/d");
tout.Branch("e_dc_z", &e_dc_z,"e_dc_z[3]/d");
tout.Branch("e_edge_dc", &e_edge_dc);
tout.Branch("e_sec_dc", &e_sec_dc);
tout.Branch("e_V", &e_V);
tout.Branch("e_W", &e_W);
tout.Branch("e_Chi2DoF", &e_Chi2DoF);
tout.Branch("isGoodElectron", &isGoodElectron);
tout.Branch("ft_x", &ft_x);
tout.Branch("ft_y", &ft_y);
tout.Branch("ft_r", &ft_r);

tout.Branch("e_htcc_x", &e_htcc_x);
tout.Branch("e_htcc_y", &e_htcc_y);
tout.Branch("e_htcc_z", &e_htcc_z);
tout.Branch("q", &q);
tout.Branch("theta_q", &theta_q);
tout.Branch("Q2", &Q2);
tout.Branch("xB", &xB);
tout.Branch("W2", &W2);
tout.Branch("Yscale", &Yscale);
tout.Branch("runCharge", &runCharge);

/*
tout.Branch("mult_p_cd", &mult_p_cd);
tout.Branch("mult_p_fd", &mult_p_fd);
tout.Branch("mult_p", &mult_p);
tout.Branch("p_fd", &p_fd);
tout.Branch("p_cd", &p_cd);
tout.Branch("tcorr", &tcorr);
tout.Branch("pbeta", &pbeta);
tout.Branch("ptof", &ptof);
tout.Branch("ptofp", &ptofp);
tout.Branch("ppath", &ppath);
tout.Branch("p_dE_cnd1", &p_dE_cnd1);
tout.Branch("p_dE_cnd2", &p_dE_cnd2);
tout.Branch("p_dE_cnd3", &p_dE_cnd3);
tout.Branch("p_dE_ctof", &p_dE_ctof);
tout.Branch("p_tof_cnd1", &p_tof_cnd1);
tout.Branch("p_tof_cnd2", &p_tof_cnd2);
tout.Branch("p_tof_cnd3", &p_tof_cnd3);
tout.Branch("p_tof_ctof", &p_tof_ctof);
tout.Branch("p_path_cnd1", &p_path_cnd1);
tout.Branch("p_path_cnd2", &p_path_cnd2);
tout.Branch("p_path_cnd3", &p_path_cnd3);
tout.Branch("p_path_ctof", &p_path_ctof);
tout.Branch("Pp", &Pp);
tout.Branch("pp", &pp);
tout.Branch("Plead", &Plead);
tout.Branch("theta_p", &theta_p);
tout.Branch("phi_p", &phi_p);
tout.Branch("sector_p_cd", &sector_p_cd);
tout.Branch("vtxp", &vtxp);
tout.Branch("plead_index", &plead_index);
tout.Branch("plead_index_reco", &plead_index_reco);
tout.Branch("plead", &plead);
tout.Branch("plead_no", &plead_no);
tout.Branch("plead_no_cd", &plead_no_cd);
tout.Branch("plead_no_fd", &plead_no_fd);

tout.Branch("p_no", &p_no);
tout.Branch("p_no_cd", &p_no_cd);
tout.Branch("p_no_fd", &p_no_fd);

tout.Branch("theta_pq", &theta_pq);
tout.Branch("pq", &pq);
tout.Branch("Pmiss", &Pmiss);
tout.Branch("pmiss", &pmiss);
tout.Branch("Mmiss", &Mmiss);
tout.Branch("Mmiss_fd", &Mmiss_fd);
tout.Branch("Mmiss_cd", &Mmiss_cd);
tout.Branch("Mmiss2", &Mmiss2);
tout.Branch("Emiss", &Emiss);
tout.Branch("Emiss0", &Emiss0);
tout.Branch("alpha_miss", &alpha_miss);
tout.Branch("theta_pmq", &theta_pmq);
tout.Branch("pmissq_para", &pmissq_para);
tout.Branch("pmissq_perp", &pmissq_perp);

tout.Branch("mult_pip", &mult_pip);
tout.Branch("mult_pim", &mult_pim);
tout.Branch("mult_g", &mult_g);
tout.Branch("pim_dE_cnd1", &pim_dE_cnd1);
tout.Branch("pim_dE_cnd2", &pim_dE_cnd2);
tout.Branch("pim_dE_cnd3", &pim_dE_cnd3);
tout.Branch("pim_dE_ctof", &pim_dE_ctof);
tout.Branch("pim_p", &pim_p);
*/


tout.Branch("weight", &weight);
/*
tout.Branch("Pe_true", &Pe_true);
tout.Branch("theta_e_true", &theta_e_true);
tout.Branch("phi_e_true", &phi_e_true);
tout.Branch("q_true", &q_true);
tout.Branch("theta_q_true", &theta_q_true);
tout.Branch("Q2_true", &Q2_true);
tout.Branch("xB_true", &xB_true);
tout.Branch("W2_true", &W2_true);
*/
/*
tout.Branch("mult_p_true", &mult_p_true);
tout.Branch("Pp_true", &Pp_true);
tout.Branch("pp_true", &pp_true);
//tout.Branch("PN_true", &PN_true);
tout.Branch("Plead_true", &Plead_true);
tout.Branch("theta_p_true", &theta_p_true);
tout.Branch("phi_p_true", &phi_p_true);
tout.Branch("plead_no_true", &plead_no_true);
tout.Branch("plead_index_true", &plead_index_true);

tout.Branch("theta_pq_true", &theta_pq_true);
tout.Branch("pq_true", &pq_true);
tout.Branch("Pmiss_true", &Pmiss_true);
tout.Branch("pmiss_true", &pmiss_true);
tout.Branch("Mmiss_true", &Mmiss_true);
//tout.Branch("Mmiss2_true", &Mmiss2_true);
tout.Branch("Emiss_true", &Emiss_true);
*/



#endif //__CA_TREE_DIS_HH__
