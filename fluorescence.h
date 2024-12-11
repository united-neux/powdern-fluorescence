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

/* ok = fluo_get_material(material, formula)
 * extracts material atoms from file header
 * the result is concatenated into 'formula'
 */
int fluo_get_material(char *filename, char *formula);
