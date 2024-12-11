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
double XRMC_CrossSections(int Z, double E0, double *xs) {
  int    i_line;
  
  if (xs == NULL) return 0;
  
  // loop on possible fluorescence lines
  for (i_line=0; i_line<XRAYLIB_LINES_MAX; i_line++) { 
    // cumulative sum of the line cross sections
    xs[FLUORESCENCE] += CSb_FluorLine(Z, -i_line, E0, NULL); /* XRayLib */
  }

  // coherent and incoherent cross sections
  xs[RAYLEIGH] = CSb_Rayl( Z, E0, NULL);
  xs[COMPTON]  = CSb_Compt(Z, E0, NULL);
  
  // total interaction cross section, should converge to CSb_Total(Z, E0, NULL)
  return xs[FLUORESCENCE] + xs[RAYLEIGH] + xs[COMPTON];
} // XRMC_CrossSections

/* XRMC_SelectFromDistribution: select a random element from a distribution
 *   index = XRMC_SelectFromDistribution(cum_sum[N], N);
 * index is returned within 0 and N-1
 * The x_arr must be a continuously increasing cumulated sum, which last element is the max
 */
int XRMC_SelectFromDistribution(double x_arr[], int N)
{
  double x=rand01()*x_arr[N-1];
  if (x<x_arr[0]) {    // x is smaller than lower limit
    return 0;
  }
  if (x>=x_arr[N-1]) { // x is greater or equal to upper limit
    return N-1;
  }
  int id=0, iu=N-1; // lower and upper index of the subarray to search
  while (iu-id>1) { // search until the size of the subarray to search is >1
    int im = (id + iu)/2; // use the midpoint for equal partition
    // decide which subarray to search
    if (x>=x_arr[im]) id=im; // change min index to search upper subarray
    else iu=im; // change max index to search lower subarray
  }

  return id;
} // XRMC_SelectFromDistribution

/* XRMC_SelectInteraction: select interaction type Fluo/Compton/Rayleigh
 * Return the interaction type from a random choice within cross sections 'xs'
 *   type = XRMC_SelectInteraction(xs[3]);
 * 'xs' is computed with XRMC_CrossSections.
 * type is one of FLUORESCENCE | RAYLEIGH | COMPTON
 */
int XRMC_SelectInteraction(double *xs)
{
  double sum_xs, cum_xs[4];
  int    i;
  
  cum_xs[0]=sum_xs=0;
  for (i=0; i< 3; i++) {
    sum_xs += xs[i];
    cum_xs[i+1]= sum_xs;
  }
  return XRMC_SelectFromDistribution(cum_xs, 4);
} // XRMC_SelectInteraction

/* XRMC_SelectFluorescenceEnergy: select outgoing fluo photon energy, when incoming with 'E0'
 *   Ef = XRMC_SelectFluorescenceEnergy(Z, E0, &dE);
 */
double XRMC_SelectFluorescenceEnergy(int Z, double E0, double *dE)
{
  int i_line;
  double sum_xs, cum_xs_lines[XRAYLIB_LINES_MAX+1];

  // compute cumulated XS for all fluo lines
  cum_xs_lines[0] = sum_xs = 0;
  for (i_line=0; i_line<XRAYLIB_LINES_MAX; i_line++) { // loop on fluorescent lines
    double xs = CSb_FluorLine(Z, -i_line, E0, NULL); /* XRayLib */
    // when a line is inactive: E=xs=0
    sum_xs += xs;
    cum_xs_lines[i_line+1] = sum_xs; // cumulative sum of their cross sections
  }
  // select randomly one of these lines
  i_line = XRMC_SelectFromDistribution(cum_xs_lines, XRAYLIB_LINES_MAX); // extract a line
  // get the K shell line width as approximation of fluorescence line width
  if (dE) *dE = AtomicLevelWidth(Z, K_SHELL, NULL); // keV

  return LineEnergy(Z, -i_line, NULL); // fluorescent line energy
} // XRMC_SelectFluorescenceEnergy

// Function removing spaces from string
char * removeSpacesFromStr(char *string)
{
  // non_space_count to keep the frequency of non space characters
  int non_space_count = 0;

  //Traverse a string and if it is non space character then, place it at index non_space_count
  for (int i = 0; string[i] != '\0'; i++)
  {
    if (isalpha(string[i]) || isdigit(string[i]))
    {
      string[non_space_count] = string[i];
      non_space_count++; //non_space_count incremented
    }    
  }

  //Finally placing final character at the string end
  string[non_space_count] = '\0';
  return string;
}

// ok = fluo_get_material(material, formula)
// extracts material atoms from file header
// the result is concatenated into 'formula'
int fluo_get_material(char *filename, char *formula) {

  int  ret = 0;
  char Line[65535];
  int flag_found_cif=0;
  FILE *file = Open_File(filename, "r", NULL);
  if (!file)
    exit(fprintf(stderr, "%s: ERROR: can not open file %s\n", 
          __FILE__, filename));

  // Read the file, and search tokens in rows
  formula[0]=0;
  while (fgets(Line, sizeof(Line), file) != NULL) {
    char *token = NULL;
    int   flag_exit=0;
    char *first_non_space=NULL;
    char *next_non_space =NULL;
    // CIF: _chemical_formula_structural 'chemical_formulae'
    // CIF: _chemical_formula_sum 'chemical_formulae'
    // LAZ/LAU: # ATOM <at> <trailing>
    // LAZ/LAU: # Atom <at> <trailing>
    // LAZ/LAU: # TITLE <at> <at> ... [ trailing...]
    // CFL: Title <chemical_formulae>
    // CFL: Atom <at> <trailing>

    // search for CIF token
    // single line: search " '\'" delimiter after CIF token, marks reading the formula
    if (!strncasecmp(Line, "_chemical_formula_structural", 28) 
        || !strncasecmp(Line, "_chemical_formula_sum", 21) 
        || !strncasecmp(Line, "_chemical_formula_moiety", 24)
        || flag_found_cif) {
      if  (flag_found_cif) { flag_found_cif=0; /* can not span on more that 2 lines */} 
      else flag_found_cif=1;
      // search for delimiter after the CIF token
      char *first_space_after_token=strpbrk(Line, " \'\n");
      if (first_space_after_token) {
        // search for the characters that may compose the formula
        first_non_space = strpbrk(first_space_after_token, "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()[]");
        if (strchr(first_space_after_token,'?')) {
          // we have an invalid/unknown formula in the CIF. Skip CIF token/line
          flag_found_cif=0;
          continue;
        }
      }
      if (first_non_space) 
        next_non_space  = strpbrk(first_non_space+1, "\'\n\r"); // position of formula end
      if (next_non_space) { 
        flag_exit = 1;
        token = first_non_space;
      }
    } else if (!strncasecmp(Line, "# TITLE", 7) && strchr(Line, '['))
      token = Line+7;
    else if (!strncasecmp(Line, "# Atom", 6))
      token = Line+6;
    else if (!strncasecmp(Line, "Atom", 4))
      token = Line+4;
    else if (!strncasecmp(Line, "Title", 5))
      token = Line+7;
    if (!token) continue;
    if (!strncasecmp(Line, "# TITLE", 7)) {
      first_non_space = Line+7;
      next_non_space  = strchr(Line+7, '[');
      if (next_non_space) flag_exit = 1;
    }
    if (!first_non_space) first_non_space = strtok(token, " \t\n");
    if (!first_non_space) continue;
    if (!next_non_space)  next_non_space  = strtok(NULL,  " \t\n"); // end of formulae
    if (!next_non_space)  next_non_space = Line+strlen(Line);
    // remove spaces
    strncat(formula, first_non_space, next_non_space-first_non_space);
    ret++;
    if (flag_exit) break;
  }
  fclose(file);
  removeSpacesFromStr(formula);
  return(ret);
} // fluo_get_material
