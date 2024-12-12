#ifndef __NEUTRON_SKIM_HH__
#define __NEUTRON_SKIM_HH__

const int n_allowed_hits = 10;
const int n_allowed_particles = 5;

/// PARTICLE VARIABLES ///

Int_t idx_recon[n_allowed_particles];

Int_t pid[n_allowed_particles]; 
Float_t momentum[n_allowed_particles];
Float_t theta[n_allowed_particles];
Float_t phi[n_allowed_particles];
Int_t sector[n_allowed_particles];
Int_t charge[n_allowed_particles];
Float_t beta_arr[n_allowed_particles];
Float_t path_length[n_allowed_particles];

Int_t pid_mc[n_allowed_particles]; 
Float_t momentum_mc[n_allowed_particles];
Float_t theta_mc[n_allowed_particles];
Float_t phi_mc[n_allowed_particles];

Float_t de_dx[n_allowed_particles];
Float_t layermult[n_allowed_particles];

Bool_t inCND[n_allowed_particles];
Bool_t inCND1[n_allowed_particles];
Bool_t inCND2[n_allowed_particles];
Bool_t inCND3[n_allowed_particles];

/// RAW HIT VARIABLES ///

Int_t cnd_hit_id_arr[n_allowed_hits];
Int_t cnd_hit_status_arr[n_allowed_hits];
Int_t cnd_hit_trkID_arr[n_allowed_hits];
Int_t cnd_hit_sector_arr[n_allowed_hits];
Int_t cnd_hit_layer_arr[n_allowed_hits];
Int_t cnd_hit_component_arr[n_allowed_hits];

Float_t cnd_hit_energy_arr[n_allowed_hits];
Float_t cnd_hit_time_arr[n_allowed_hits];
Float_t cnd_hit_energy_unc_arr[n_allowed_hits];
Float_t cnd_hit_time_unc_arr[n_allowed_hits];
Float_t cnd_hit_x_arr[n_allowed_hits];
Float_t cnd_hit_y_arr[n_allowed_hits];
Float_t cnd_hit_z_arr[n_allowed_hits];
Float_t cnd_hit_x_unc_arr[n_allowed_hits];
Float_t cnd_hit_y_unc_arr[n_allowed_hits];
Float_t cnd_hit_z_unc_arr[n_allowed_hits];
Float_t cnd_hit_tx_arr[n_allowed_hits];
Float_t cnd_hit_ty_arr[n_allowed_hits];
Float_t cnd_hit_tz_arr[n_allowed_hits];
Float_t cnd_hit_tlength_arr[n_allowed_hits];
Float_t cnd_hit_pathlength_arr[n_allowed_hits];

Int_t cnd_hit_indexLadc_arr[n_allowed_hits];
Int_t cnd_hit_indexRadc_arr[n_allowed_hits];
Int_t cnd_hit_indexLtdc_arr[n_allowed_hits];
Int_t cnd_hit_indexRtdc_arr[n_allowed_hits];



#endif 
