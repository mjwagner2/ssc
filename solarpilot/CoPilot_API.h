/**
BSD 3-Clause License

Copyright (c) Alliance for Sustainable Energy, LLC. See also https://github.com/NREL/ssc/blob/develop/LICENSE
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#ifndef _COPILOT_API_
#define _COPILOT_API_



#include <string>
#include <ctime>
#include "SolarField.h"
#include "interop.h"

#if defined(_WINDLL)
#define SPEXPORT __declspec(dllexport)
#else
#define SPEXPORT
#endif
#ifdef __cplusplus
extern "C" {
#endif

    /** An opaque reference to a structure that holds a collection of variables.  This structure can contain any number of variables referenced by name, and can hold strings, numbers, arrays, and matrices.  Matrices are stored in row-major order, where the array size is nrows*ncols, and the array index is calculated by r*ncols+c. An ssc_data_t object holds all input and output variables for a simulation. It does not distinguish between input, output, and input variables - that is handled at the model context level. */
    typedef void* sp_data_t;

    /** The numeric type used in the SolarPILOT API. */ 
    typedef double sp_number_t;

    SPEXPORT const char* sp_version(sp_data_t p_data);

    SPEXPORT void sp_set_callback(sp_data_t p_data, int (*)(sp_number_t, const char*));

    SPEXPORT void sp_disable_callback(sp_data_t p_data);

    SPEXPORT void sp_cancel_simulation(sp_data_t p_data);

    /** Creates a new data object in memory.  A data object stores a table of named values, where each value can be of any SolarPILOT datatype. */
    SPEXPORT sp_data_t sp_data_create();

    /** Frees the memory associated with a data object, where p_data is the data container to free. */
    SPEXPORT bool sp_data_free(sp_data_t p_data);

    SPEXPORT void var_free_memory(sp_number_t* varptr);

    SPEXPORT bool sp_set_number(sp_data_t p_data, const char* name, sp_number_t v);

    SPEXPORT bool sp_set_string(sp_data_t p_data, const char *name, const char *value);

    SPEXPORT bool sp_set_array(sp_data_t p_data, const char *name, sp_number_t *pvalues, int length);

    /** Assigns value of type @a SSC_MATRIX . Matrices are specified as a continuous array, in row-major order.  Example: the matrix [[5,2,3],[9,1,4]] is stored as [5,2,3,9,1,4]. */
    SPEXPORT bool sp_set_matrix(sp_data_t p_data, const char *name, sp_number_t *pvalues, int nrows, int ncols);

    SPEXPORT sp_number_t sp_get_number(sp_data_t p_data, const char* name);

    /** Returns the value of a @a SSC_STRING variable with the given name. */
    SPEXPORT const char *sp_get_string(sp_data_t p_data, const char *name);

    /** Returns the value of a @a SSC_ARRAY variable with the given name. */
    SPEXPORT sp_number_t* sp_get_array(sp_data_t p_data, const char* name, int* length);
    //SPEXPORT void sp_get_array(sp_data_t p_data, const char *name, sp_number_t* values, int *length);

    SPEXPORT int sp_get_array_size(sp_data_t p_data, const char* name);

    SPEXPORT double sp_get_array_value_by_index(sp_data_t p_data, const char* name, int index);

    SPEXPORT sp_number_t* sp_get_matrix(sp_data_t p_data, const char* name, int* ncols, int* nrows);

    SPEXPORT int sp_get_matrix_row_size(void* p_data, const char* name);

    SPEXPORT int sp_get_matrix_col_size(void* p_data, const char* name);

    SPEXPORT double sp_get_matrix_value_by_index(void* p_data, const char* name, int row_index, int col_index);


    SPEXPORT void sp_reset_geometry(sp_data_t p_data);

    SPEXPORT int sp_add_receiver(sp_data_t p_data, const char* receiver_name);

    SPEXPORT bool sp_drop_receiver(sp_data_t p_data, const char* receiver_name);

    SPEXPORT int sp_add_heliostat_template(sp_data_t p_data, const char* heliostat_name);

    SPEXPORT bool sp_drop_heliostat_template(sp_data_t p_data, const char* heliostat_name);

    SPEXPORT sp_number_t* sp_generate_simulation_days(sp_data_t p_data, int* nrecord, int* ncol);

    SPEXPORT bool sp_update_geometry(sp_data_t p_data);

    SPEXPORT bool sp_generate_layout(sp_data_t p_data, int nthreads);

    SPEXPORT bool sp_assign_layout(sp_data_t p_data, sp_number_t* pvalues, int nrows, int ncols, int nthreads);

    //SPEXPORT sp_number_t* sp_get_layout_info(sp_data_t p_data, int* nhelio, int* ncol);
    SPEXPORT sp_number_t* sp_get_layout_info(sp_data_t p_data, int* nhelio, int* ncol, bool get_corners, bool get_optical_details);

    SPEXPORT const char* sp_get_layout_header(sp_data_t p_data, bool get_corners, bool get_optical_details);

    SPEXPORT bool sp_simulate(sp_data_t p_data, int nthreads, bool update_aimpoints); //bool save_detail,

    SPEXPORT double sp_get_receiver_power(sp_data_t p_data);

    SPEXPORT const char* sp_summary_results(sp_data_t p_data);

    SPEXPORT sp_number_t* sp_detail_results(sp_data_t p_data, int* nrows, int* ncols, sp_number_t* selhel, int nselhel, bool get_corners);
    //SPEXPORT sp_number_t* sp_detail_results(sp_data_t p_data, int* nrows, int* ncols, const char* header, sp_number_t* selhel, int nselhel);

    SPEXPORT const char* sp_detail_results_header(sp_data_t p_data, bool get_corners);

    SPEXPORT sp_number_t* sp_get_fluxmap(sp_data_t p_data, int* nrows, int* ncols, int rec_id, int flux_surf);

    SPEXPORT void sp_optimize(sp_data_t p_data, sp_number_t* pvalues, int nvar);

    SPEXPORT void sp_clear_land(sp_data_t p_data, const char* type);

    SPEXPORT bool sp_add_land(sp_data_t p_data, const char* type, sp_number_t* polygon_points, int* npts, int* ndim, bool is_append);

    SPEXPORT sp_number_t* sp_heliostats_by_region(sp_data_t p_data, const char* coor_sys, int* lenret,
                                    sp_number_t* arguments, int* len_arg, const char* svgfname_data, sp_number_t* svg_opt_tab);

    SPEXPORT bool sp_modify_heliostats(sp_data_t p_data, sp_number_t* helio_data, int* nhel, int* ncols, const char* table_hdr);

    SPEXPORT bool sp_calculate_optical_efficiency_table(sp_data_t p_data, int ud_n_az, int ud_n_zen);

    SPEXPORT sp_number_t* sp_get_optical_efficiency_table(sp_data_t p_data, int* nrows, int* ncols);

    SPEXPORT bool sp_calculate_get_optical_efficiency_table(sp_data_t p_data, const size_t ud_n_az, const size_t ud_n_zen,
        double* azimuths, double* zenith, double* eff_matrix);

    SPEXPORT bool sp_save_optical_efficiency_table(sp_data_t p_data, const char* sp_fname, const char* table_name);

    SPEXPORT bool sp_save_from_script(sp_data_t p_data, const char* sp_fname);

    SPEXPORT bool sp_load_from_script(sp_data_t p_data, const char* sp_fname);

    SPEXPORT bool sp_dump_varmap(sp_data_t p_data, const char* sp_fname);

    SPEXPORT bool sp_export_soltrace(sp_data_t p_data, const char* sp_fname);

    SPEXPORT bool sp_load_soltrace_context(sp_data_t p_data, st_context_t* solt_cxt);

    SPEXPORT void _sp_free_var(sp_number_t* m);

#ifdef __cplusplus
}
#endif


#endif  // _COPILOT_API_

extern int ST_APICallback(st_uint_t ntracedtotal, st_uint_t ntraced, st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages, void* data);

extern int MessageHandler(const char* message, void* data);

extern int ProgressHandler(double progress, const char* message, void* data);