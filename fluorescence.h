/* See XrayLib code (c) T. Schoonjans
 * /usr/include/xraylib/xraylib-parser.h      for compoundData
 * /usr/include/xraylib/xraylib.h             for XS
 */

/* inspired from:
 * https://github.com/golosio/xrmc src/photon/photon.cpp (c) Bruno Golosio    
 */
  
/* XRMC_CrossSections: Compute interaction cross sections in [barn/atom]
 * Return total cross section, given Z and E0:
 *   total_xs = XRMC_CrossSections(Z, E0, xs[3]);
 */
double XRMC_CrossSections(int Z, double E0, double *xs);

/* XRMC_SelectFromDistribution: select a random element from a distribution
 *   index = XRMC_SelectFromDistribution(cum_sum[N], N);
 * index is returned within 0 and N-1
 * The x_arr must be a continuously increasing cumulated sum, which last element is the max
 */
int XRMC_SelectFromDistribution(double x_arr[], int N);

/* XRMC_SelectInteraction: select interaction type Fluo/Compton/Rayleigh
 * Return the interaction type from a random choice within cross sections 'xs'
 *   type = XRMC_SelectInteraction(xs[3]);
 * 'xs' is computed with XRMC_CrossSections.
 * type is one of FLUORESCENCE | RAYLEIGH | COMPTON
 */
int XRMC_SelectInteraction(double *xs);

/* XRMC_SelectFluorescenceEnergy: select outgoing fluo photon energy, when incoming with 'E0'
 *   Ef = XRMC_SelectFluorescenceEnergy(Z, E0, &dE);
 */
double XRMC_SelectFluorescenceEnergy(int Z, double E0, double *dE);

/* Function removing spaces from string */
char * removeSpacesFromStr(char *string);

int fluo_PN_list_compare (void const *a, void const *b);

struct fluo_line_data
    {
      double F2;                  /* Value of structure factor */
      double q;                   /* Qvector */
      int j;                      /* Multiplicity */
      double DWfactor;            /* Debye-Waller factor */
      double w;                   /* Intrinsic line width */
      double Epsilon;             /* Strain=delta_d_d/d shift in ppm */
    };


struct fluo_line_info_struct {
  struct fluo_line_data *list;     /* Reflection array */
  int  count;                  /* Number of reflections */
  double Dd;
  double DWfactor;
  double V_0;
  double pow_density;//TODO this is likely not necessary - so check this
  double at_weight;
  double at_nb;
  double sigma_a;
  double sigma_i;
  char   compname[256];
  double flag_barns;
  int    shape; /* 0 cylinder, 1 box, 2 sphere, 3 OFF file */
  int    column_order[9]; /* column signification */
  int    flag_warning;
  char   type;  /* interaction type of event t=Transmit, i=Incoherent, c=Coherent */
  double dq;    /* wavevector transfer [Angs-1] */
  double Epsilon; /* global strain in ppm */
  double XsectionFactor;
  double my_s_k2_sum;
  double my_a;
  double my_inc;
  double lfree; // store mean free path for the last event;
  double *w,*q, *my_s_k2;
  double radius_i,xwidth_i,yheight_i,zdepth_i;
  double k; /* last wavenumber (cached) */
  double Nq;
  double xs_Nq[CHAR_BUF_LENGTH];
  double xs_sum[CHAR_BUF_LENGTH];
  int    nb_reuses, nb_refl, nb_refl_count;
  double k_min, k_max;
  unsigned int photon_passed;
  long   xs_compute, xs_reuse, xs_calls;
  t_Table mat_table;
  int mat_column_order[5]; /*column signification for the coeff. in material data file*/
};

/* ok = fluo_get_material(material, formula)
 * extracts material atoms from file header
 * the result is concatenated into 'formula'
 */
int fluo_get_material(char *filename, char *formula);

int fluo_calc_xsect(double k, double *q, double *my_s_k2, int count, double *sum,
          struct fluo_line_info_struct *line_info);

int fluo_read_line_data(char *reflections, struct fluo_line_info_struct *info);


int XRMC_SelectPowderLineQ(struct fluo_line_info_struct *line_info, double Ei, double *Q);
