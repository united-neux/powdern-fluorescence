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

int fluo_read_line_data(char *SC_file, struct line_info_struct_union *info) {
  struct line_data_union *list = NULL;
  int    size = 0;
  t_Table sTable; /* sample data table structure from SC_file */
  int    i=0;
  int    mult_count  =0;
  char   flag=0;
  double q_count=0, j_count=0, F2_count=0;
  char **parsing;
  int    list_count=0;
  double sum_F2=0;
  char  *filename=NULL;

  if (!SC_file || !strlen(SC_file) || !strcmp(SC_file, "NULL")) {
    printf("PowderN: %s: Using incoherent elastic scattering only.\n",
        info->compname);
    info->count = 0;
    return(0);
  }
  filename = cif2hkl(SC_file, "--mode XRA");
  long retval = Table_Read(&sTable, filename, 1); /* read 1st block data from SC_file into sTable*/
  if (retval<0) {
    fprintf(stderr,"PowderN: Could not open file %s - exiting!\n", filename);
    exit(-1);
  }

  /* parsing of header */
  parsing = Table_ParseHeader(sTable.header,
      "Vc","V_0",
      "column_j",
      "column_d",
      "column_F2",
      "column_DW",
      "column_Dd",
      "column_inv2d", "column_1/2d", "column_sintheta/lambda",
      "column_q", /* 10 */
      "DW", "Debye_Waller",
      "delta_d/d",
      "column_F ",
      "V_rho",
      "density",
      "weight",
      "nb_atoms","multiplicity",
      NULL);

  if (parsing) {
    if (parsing[0] && !info->V_0)     info->V_0    =atof(parsing[0]);
    if (parsing[1] && !info->V_0)     info->V_0    =atof(parsing[1]);
    if (parsing[2])                   info->column_order[0]=atoi(parsing[2]);
    if (parsing[3])                   info->column_order[1]=atoi(parsing[3]);
    if (parsing[4])                   info->column_order[2]=atoi(parsing[4]);
    if (parsing[5])                   info->column_order[3]=atoi(parsing[5]);
    if (parsing[6])                  info->column_order[4]=atoi(parsing[6]);
    if (parsing[7])                  info->column_order[5]=atoi(parsing[7]);
    if (parsing[8])                  info->column_order[5]=atoi(parsing[8]);
    if (parsing[9])                  info->column_order[5]=atoi(parsing[9]);
    if (parsing[10])                  info->column_order[6]=atoi(parsing[10]);
    if (parsing[11] && info->DWfactor<=0)    info->DWfactor=atof(parsing[11]);
    if (parsing[12] && info->DWfactor<=0)    info->DWfactor=atof(parsing[12]);
    if (parsing[13] && info->Dd <0)          info->Dd      =atof(parsing[13]);
    if (parsing[14])                  info->column_order[7]=atoi(parsing[14]);
    if (parsing[15] && !info->V_0)    info->V_0    =1/atof(parsing[15]);
    if (parsing[16] && !info->rho)    info->rho    =atof(parsing[16]);
    if (parsing[17] && !info->at_weight)     info->at_weight    =atof(parsing[17]);
    if (parsing[18] && info->at_nb <= 1)  info->at_nb    =atof(parsing[18]);
    if (parsing[19] && info->at_nb <= 1)  info->at_nb    =atof(parsing[19]);
    for (i=0; i<=19; i++) if (parsing[i]) free(parsing[i]);
    free(parsing);
  }

  if (!sTable.rows)
    exit(fprintf(stderr, "PowderN: %s: Error: The number of rows in %s "
          "should be at least %d\n", info->compname, SC_file, 1));
  else
    size = sTable.rows;

  Table_Info(sTable);
  printf("PowderN: %s: Reading %d rows from %s\n",
      info->compname, size, SC_file);

  if (info->column_order[0] == 4 && info->flag_barns !=0)
    printf("PowderN: %s: Powder file probably of type Crystallographica/Fullprof (lau)\n"
        "WARNING: but F2 unit is set to barns=1 (barns). Intensity might be 100 times too high.\n",
        info->compname);
  if (info->column_order[0] == 17 && info->flag_barns == 0)
    printf("PowderN: %s: Powder file probably of type Lazy Pulver (laz)\n"
        "WARNING: but F2 unit is set to barns=0 (fm^2). Intensity might be 100 times too low.\n",
        info->compname);
  /* allocate line_data array */
  list = (struct line_data_union*) calloc(size, sizeof(struct line_data_union));

  for (i=0; i<size; i++)
  {
    /*      printf("Reading in line %i\n",i);*/
    double j=0, d=0, w=0, q=0, DWfactor=0, F2=0, Epsilon=0;
    int index;

    if (info->Dd >= 0)      w         = info->Dd;
    if (info->DWfactor > 0) DWfactor  = info->DWfactor;
    if (info->Epsilon)      Epsilon   = info->Epsilon*1e-6;

    /* get data from table using columns {j d F2 DW Dd inv2d q F} */
    /* column indexes start at 1, thus need to substract 1 */
    if (info->column_order[0] >0)
      j = Table_Index(sTable, i, info->column_order[0]-1);
    if (info->column_order[1] >0)
      d = Table_Index(sTable, i, info->column_order[1]-1);
    if (info->column_order[2] >0)
      F2 = Table_Index(sTable, i, info->column_order[2]-1);
    if (info->column_order[3] >0)
      DWfactor = Table_Index(sTable, i, info->column_order[3]-1);
    if (info->column_order[4] >0)
      w = Table_Index(sTable, i, info->column_order[4]-1);
    if (info->column_order[5] >0 && !(info->column_order[1] >0)) // Only use if d not read already
    { d = Table_Index(sTable, i, info->column_order[5]-1);
      d = (d > 0? 1/d/2 : 0); }
    if (info->column_order[6] >0 && !(info->column_order[1] >0)) // Only use if d not read already
    { q = Table_Index(sTable, i, info->column_order[6]-1);
      d = (q > 0 ? 2*PI/q : 0); }
    if (info->column_order[7] >0  && !F2)
    { F2 = Table_Index(sTable, i, info->column_order[7]-1); F2 *= F2; }
    if (info->column_order[8] >0  && !Epsilon)
    { Epsilon = Table_Index(sTable, i, info->column_order[8]-1)*1e-6; }

    /* assign and check values */
    j        = (j > 0 ? j : 0);
    q        = (d > 0 ? 2*PI/d : 0); /* this is q */
    if (Epsilon && fabs(Epsilon) < 1e6) {
      q     -= Epsilon*q; /* dq/q = -delta_d_d/d = -Epsilon */
    }
    DWfactor = (DWfactor > 0 ? DWfactor : 1);
    w        = (w>0 ? w : 0); /* this is q and d relative spreading */
    F2       = (F2 >= 0 ? F2 : 0);
    if (j == 0 || q == 0) {
      printf("PowderN: %s: line %i has invalid definition\n"
          "         (mult=0 or q=0 or d=0)\n", info->compname, i);
      continue;
    }
    list[list_count].j = j;
    list[list_count].q = q;
    list[list_count].DWfactor = DWfactor;
    list[list_count].w = w;
    list[list_count].F2= F2;
    list[list_count].Epsilon = Epsilon;
    sum_F2 += F2;

    /* adjust multiplicity if j-column + multiple d-spacing lines */
    /* if  d = previous d, increase line duplication index */
    if (!q_count)      q_count  = q;
    if (!j_count)      j_count  = j;
    if (!F2_count)     F2_count = F2;
    if (fabs(q_count-q) < 0.0001*fabs(q)
        && fabs(F2_count-F2) < 0.0001*fabs(F2) && j_count == j) {
      mult_count++; flag=0; }
    else flag=1;
    if (i == size-1) flag=1;
    /* else if d != previous d : just passed equivalent lines */
    if (flag) {
      if (i == size-1) list_count++;
      /*   if duplication index == previous multiplicity */
      /*      set back multiplicity of previous lines to 1 */
      if ((mult_count && list_count>0)
          && (mult_count == list[list_count-1].j
            || ((list_count < size) && (i == size - 1)
              && (mult_count == list[list_count].j))) ) {
        printf("PowderN: %s: Set multiplicity to 1 for lines [%i:%i]\n"
            "         (d-spacing %g is duplicated %i times)\n",
            info->compname, list_count-mult_count, list_count-1, list[list_count-1].q, mult_count);
        for (index=list_count-mult_count; index<list_count; list[index++].j = 1);
        mult_count = 1;
        q_count   = q;
        j_count   = j;
        F2_count  = F2;
      }
      if (i == size-1) list_count--;
      flag=0;
    }
    list_count++;
  } /* end for */

  Table_Free(&sTable);

  if (!sum_F2) {
    MPI_MASTER(
        fprintf(stderr, "PowderN: %s: ERROR: all %i structure factors in file '%s' are null. Check the reflection list.\n",
          info->compname, list_count, SC_file);
        );
    return(0);
  }

  /* sort the list with increasing q */
  qsort(list, list_count, sizeof(struct line_data_union),  PN_list_compare_union);

  printf("PowderN: %s: Read %i reflections from file '%s'\n",
      info->compname, list_count, SC_file);
  // remove temporary F2(hkl) file when giving CFL/CIF/ShelX file
  if (filename && filename != SC_file)
    unlink(filename);

  info->list  = list;
  info->count = list_count;

  return(list_count);
} /* read_line_data_union */


/* computes the number of possible reflections (return value), and the total xsection 'sum' */
/* this routine looks for a pre-computed value in the Nq and sum cache tables               */
/* when found, the earch starts from the corresponding lower element in the table           */
int fluo_calc_xsect(double k, double *q, double *my_s_k2, int count, double *sum,
    struct line_info_struct_union *line_info) {
  int    Nq = 0, line=0, line0=0;
  /*sinth for tthmax=180 is 1.*/
  double sinth=1.0;//sin(DEG2RAD*tth_max*0.5);
  *sum=0;

  /* check if a line_info element has been recorded already - not on OpenACC */
  if (k >= line_info->k_min && k <= line_info->k_max && line_info->photon_passed >= CHAR_BUF_LENGTH) {
    line = (int)floor(k - line_info->k_min)*CHAR_BUF_LENGTH/(line_info->k_max - line_info->k_min);
    Nq    = line_info->xs_Nq[line];
    *sum  = line_info->xs_sum[line];
    if (!Nq && *sum == 0) {
      /* not yet set: we compute the sum up to the corresponding wavevector in the table cache */
      double line_k = line_info->k_min + line*(line_info->k_max - line_info->k_min)/CHAR_BUF_LENGTH;
      for(line0=0; line0<count; line0++) {
        if (q[line0] <= 2*line_k*sinth) { /* q < 2*kf: restrict structural range */
          *sum += my_s_k2[line0];
          if (Nq < line0+1) Nq=line0+1; /* determine maximum line index which can scatter */
        } else break;
      }
      line_info->xs_Nq[line] = Nq;
      line_info->xs_sum[line]= *sum;
      line_info->xs_compute++;
    } else line_info->xs_reuse++;
    line0 = Nq;
  }

  line_info->xs_calls++;

  for(line=line0; line<count; line++) {
    if (q[line] <= 2*k*sinth) { /* q < 2*kf: restrict structural range */
      *sum += my_s_k2[line];
      if (Nq < line+1) Nq=line+1; /* determine maximum line index which can scatter */
    } else break;
  }

  return(Nq);
} /* fluo_calc_xsect */

