/*
 * MATLAB Compiler: 4.10 (R2009a)
 * Date: Fri Jun 18 03:17:45 2010
 * Arguments: "-B" "macro_default" "-m" "-W" "main" "-T" "link:exe" "-v" "-R"
 * "nojit" "-R" "nojvm" "-d"
 * "./add_diagonal_covariance_component_in_full_cov_form_buildfiles"
 * "./add_diagonal_covariance_component_in_full_cov_form.m" 
 */

#include "mclmcrrt.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_add_diagonal_covariance_component_in_full_cov_form_session_key[] = {
    'B', '8', '4', '1', '4', 'F', '7', '9', 'F', '2', '5', 'F', 'D', '6', '2',
    '8', 'B', '0', '9', '7', '2', '3', '4', 'D', '3', '3', '4', '1', '1', '3',
    '2', 'B', '6', '4', '4', 'B', '9', 'B', '1', 'E', 'A', '8', '9', 'F', '5',
    'B', '7', 'D', '9', '9', '9', '7', 'A', '7', '1', 'C', '0', 'E', 'E', '0',
    'B', '3', '0', '6', 'A', '5', 'C', '0', 'B', 'D', '9', 'B', '0', 'E', '1',
    'C', 'E', 'C', 'D', '8', '4', 'B', '0', '7', '0', '9', '3', '2', '4', '3',
    '2', 'B', 'C', '7', '4', '8', 'B', 'D', '5', '2', '9', 'A', '3', '2', '4',
    'B', '7', '0', 'E', '2', 'C', 'F', '1', '7', 'D', 'A', '2', '9', 'B', '0',
    'C', '1', '3', '2', 'F', '4', '1', '5', '1', 'B', '3', '4', 'A', 'C', '3',
    '6', '5', 'A', '0', '8', '2', '3', '9', '1', 'F', '5', '7', 'F', '5', 'D',
    'E', 'F', 'F', 'D', '9', '4', '8', '6', 'E', 'D', '4', 'E', '9', '7', 'B',
    '5', '7', '7', '9', '5', '1', '6', 'A', '2', '7', '7', '2', '5', '6', '1',
    '8', '0', 'C', '2', 'F', '1', '4', '3', 'E', 'A', '2', 'A', 'A', '6', '8',
    '1', '9', 'C', 'A', 'C', 'B', 'B', '2', '6', 'E', '3', 'F', 'B', '8', 'D',
    'B', 'D', 'E', '5', '2', 'F', 'D', '4', 'B', 'B', 'E', '7', '5', 'F', '9',
    '5', '4', 'F', 'A', 'F', '2', '3', 'B', '8', 'F', '6', '2', 'C', '7', 'A',
    'B', '2', '2', 'E', '2', '1', 'B', '9', '3', 'E', '4', 'B', '2', '6', 'A',
    'C', '\0'};

const unsigned char __MCC_add_diagonal_covariance_component_in_full_cov_form_public_key[] = {
    '3', '0', '8', '1', '9', 'D', '3', '0', '0', 'D', '0', '6', '0', '9', '2',
    'A', '8', '6', '4', '8', '8', '6', 'F', '7', '0', 'D', '0', '1', '0', '1',
    '0', '1', '0', '5', '0', '0', '0', '3', '8', '1', '8', 'B', '0', '0', '3',
    '0', '8', '1', '8', '7', '0', '2', '8', '1', '8', '1', '0', '0', 'C', '4',
    '9', 'C', 'A', 'C', '3', '4', 'E', 'D', '1', '3', 'A', '5', '2', '0', '6',
    '5', '8', 'F', '6', 'F', '8', 'E', '0', '1', '3', '8', 'C', '4', '3', '1',
    '5', 'B', '4', '3', '1', '5', '2', '7', '7', 'E', 'D', '3', 'F', '7', 'D',
    'A', 'E', '5', '3', '0', '9', '9', 'D', 'B', '0', '8', 'E', 'E', '5', '8',
    '9', 'F', '8', '0', '4', 'D', '4', 'B', '9', '8', '1', '3', '2', '6', 'A',
    '5', '2', 'C', 'C', 'E', '4', '3', '8', '2', 'E', '9', 'F', '2', 'B', '4',
    'D', '0', '8', '5', 'E', 'B', '9', '5', '0', 'C', '7', 'A', 'B', '1', '2',
    'E', 'D', 'E', '2', 'D', '4', '1', '2', '9', '7', '8', '2', '0', 'E', '6',
    '3', '7', '7', 'A', '5', 'F', 'E', 'B', '5', '6', '8', '9', 'D', '4', 'E',
    '6', '0', '3', '2', 'F', '6', '0', 'C', '4', '3', '0', '7', '4', 'A', '0',
    '4', 'C', '2', '6', 'A', 'B', '7', '2', 'F', '5', '4', 'B', '5', '1', 'B',
    'B', '4', '6', '0', '5', '7', '8', '7', '8', '5', 'B', '1', '9', '9', '0',
    '1', '4', '3', '1', '4', 'A', '6', '5', 'F', '0', '9', '0', 'B', '6', '1',
    'F', 'C', '2', '0', '1', '6', '9', '4', '5', '3', 'B', '5', '8', 'F', 'C',
    '8', 'B', 'A', '4', '3', 'E', '6', '7', '7', '6', 'E', 'B', '7', 'E', 'C',
    'D', '3', '1', '7', '8', 'B', '5', '6', 'A', 'B', '0', 'F', 'A', '0', '6',
    'D', 'D', '6', '4', '9', '6', '7', 'C', 'B', '1', '4', '9', 'E', '5', '0',
    '2', '0', '1', '1', '1', '\0'};

static const char * MCC_add_diagonal_covariance_component_in_full_cov_form_matlabpath_data[] = 
  { "add_diagonal/", "$TOOLBOXDEPLOYDIR/",
    "$TOOLBOXMATLABDIR/general/", "$TOOLBOXMATLABDIR/ops/",
    "$TOOLBOXMATLABDIR/lang/", "$TOOLBOXMATLABDIR/elmat/",
    "$TOOLBOXMATLABDIR/randfun/", "$TOOLBOXMATLABDIR/elfun/",
    "$TOOLBOXMATLABDIR/specfun/", "$TOOLBOXMATLABDIR/matfun/",
    "$TOOLBOXMATLABDIR/datafun/", "$TOOLBOXMATLABDIR/polyfun/",
    "$TOOLBOXMATLABDIR/funfun/", "$TOOLBOXMATLABDIR/sparfun/",
    "$TOOLBOXMATLABDIR/scribe/", "$TOOLBOXMATLABDIR/graph2d/",
    "$TOOLBOXMATLABDIR/graph3d/", "$TOOLBOXMATLABDIR/specgraph/",
    "$TOOLBOXMATLABDIR/graphics/", "$TOOLBOXMATLABDIR/uitools/",
    "$TOOLBOXMATLABDIR/strfun/", "$TOOLBOXMATLABDIR/imagesci/",
    "$TOOLBOXMATLABDIR/iofun/", "$TOOLBOXMATLABDIR/audiovideo/",
    "$TOOLBOXMATLABDIR/timefun/", "$TOOLBOXMATLABDIR/datatypes/",
    "$TOOLBOXMATLABDIR/verctrl/", "$TOOLBOXMATLABDIR/codetools/",
    "$TOOLBOXMATLABDIR/helptools/", "$TOOLBOXMATLABDIR/demos/",
    "$TOOLBOXMATLABDIR/timeseries/", "$TOOLBOXMATLABDIR/hds/",
    "$TOOLBOXMATLABDIR/guide/", "$TOOLBOXMATLABDIR/plottools/",
    "toolbox/local/", "toolbox/shared/dastudio/",
    "$TOOLBOXMATLABDIR/datamanager/", "toolbox/compiler/" };

static const char * MCC_add_diagonal_covariance_component_in_full_cov_form_classpath_data[] = 
  { "" };

static const char * MCC_add_diagonal_covariance_component_in_full_cov_form_libpath_data[] = 
  { "" };

static const char * MCC_add_diagonal_covariance_component_in_full_cov_form_app_opts_data[] = 
  { "" };

static const char * MCC_add_diagonal_covariance_component_in_full_cov_form_run_opts_data[] = 
  { "nojit", "nojvm" };

static const char * MCC_add_diagonal_covariance_component_in_full_cov_form_warning_state_data[] = 
  { "off:MATLAB:dispatcher:nameConflict" };


mclComponentData __MCC_add_diagonal_covariance_component_in_full_cov_form_component_data = { 

  /* Public key data */
  __MCC_add_diagonal_covariance_component_in_full_cov_form_public_key,

  /* Component name */
  "add_diagonal_covariance_component_in_full_cov_form",

  /* Component Root */
  "",

  /* Application key data */
  __MCC_add_diagonal_covariance_component_in_full_cov_form_session_key,

  /* Component's MATLAB Path */
  MCC_add_diagonal_covariance_component_in_full_cov_form_matlabpath_data,

  /* Number of directories in the MATLAB Path */
  38,

  /* Component's Java class path */
  MCC_add_diagonal_covariance_component_in_full_cov_form_classpath_data,
  /* Number of directories in the Java class path */
  0,

  /* Component's load library path (for extra shared libraries) */
  MCC_add_diagonal_covariance_component_in_full_cov_form_libpath_data,
  /* Number of directories in the load library path */
  0,

  /* MCR instance-specific runtime options */
  MCC_add_diagonal_covariance_component_in_full_cov_form_app_opts_data,
  /* Number of MCR instance-specific runtime options */
  0,

  /* MCR global runtime options */
  MCC_add_diagonal_covariance_component_in_full_cov_form_run_opts_data,
  /* Number of MCR global runtime options */
  2,
  
  /* Component preferences directory */
  "add_diagonal_6423968B99592805B52788146334593E",

  /* MCR warning status data */
  MCC_add_diagonal_covariance_component_in_full_cov_form_warning_state_data,
  /* Number of MCR warning status modifiers */
  1,

  /* Path to component - evaluated at runtime */
  NULL

};

#ifdef __cplusplus
}
#endif


