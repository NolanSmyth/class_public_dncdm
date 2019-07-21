/** @file class.c 
 * Julien Lesgourgues, 17.04.2011    
 */
 
#include "class.h"

/* this main runs only the background part */

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */


  

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }
  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  /****** here you can output the evolution of any background
	  quanitity you are interested in ******/

  int index_tau;

  FILE *fptr = fopen("output/bg_dncdm_test.dat", "w"); 

  double tau_at_z, Omega_dncdm, rho_c;

  background_tau_of_z(&ba, 0., &tau_at_z);

  /* 
   * The prefactors here come from the fact that we actually assume DNCDM is a Boson, so it follows the B-E distribution.
   * This means that the degrees of freedom is different in calculation of the temperature ratio from Neff, 
   * AND, the degrees of freedom is different in computing the final abundance.
   * The former involves the relativistic energy density, while the latter involves the number density
   */
  // This is only correct if Gamma = 0!
  // Bose-Einstein distribution
  Omega_dncdm = (8./7.)*(3./4.)*ba.background_table[(ba.bt_size-1)*ba.bg_size + ba.index_bg_rho_dncdm]/pow(ba.H0,2); 
  // FD Distribution
  Omega_dncdm = ba.background_table[(ba.bt_size-1)*ba.bg_size + ba.index_bg_rho_dncdm]/pow(ba.H0,2); 

  printf(" -> decaying non-cold dark matter species has m = %e eV (so m / omega =%e eV)\n",
             ba.m_dncdm_in_eV,
             ba.m_dncdm_in_eV/Omega_dncdm/ba.h/ba.h);
  printf(" -> Should be ~ 93 eV if Gamma_dncdm = 0\n");
  printf(" -> H0 = %e\n", ba.H0);

  for (index_tau=0; index_tau<ba.bt_size; index_tau++) {

    /*
    fprintf(stdout,
	    "tau=%e z=%e a=%e H=%e\n",
	    ba.tau_table[index_tau],
	    ba.z_table[index_tau],
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_a],
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_H]);
    */
    /*
    fprintf(fptr,//stdout,
	    "%e     %e     %e     %e     %e\n",
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_a],
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_rho_dncdm],
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_p_dncdm],
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_rho_dr], 
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_rho_ur]      
        );
    */
    fprintf(fptr,//stdout,
	    "%e     %e     %e     %e     %e     %e     ",
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_a],
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_rho_dncdm],
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_p_dncdm],
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_n_dncdm]*ba.M_dncdm,
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_rho_dr], 
	    ba.background_table[index_tau*ba.bg_size+ba.index_bg_rho_ur]      
        );
    for (int index_q = 0; index_q < ba.q_size_dncdm_bg; index_q++)
      fprintf(fptr, "%e     ", exp(ba.background_table[index_tau*ba.bg_size+ba.index_bg_lnf_dncdm + index_q]));
    for (int index_q = 0; index_q < ba.q_size_dncdm_bg; index_q++)
      fprintf(fptr, "%e     ", ba.background_table[index_tau*ba.bg_size+ba.index_bg_dlnf_dlnq_dncdm + index_q]);
    for (int index_q = 0; index_q < ba.q_size_dncdm; index_q++)
      fprintf(fptr, "%e     ", exp(ba.background_table[index_tau*ba.bg_size+ba.index_bg_lnf_pt_dncdm + index_q]));
    for (int index_q = 0; index_q < ba.q_size_dncdm; index_q++)
      fprintf(fptr, "%e     ", ba.background_table[index_tau*ba.bg_size+ba.index_bg_dlnf_dlnq_pt_dncdm + index_q]);
    fprintf(fptr,"\n");
  }
  fclose(fptr);

  /****** all calculations done, now free the structures ******/

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}
