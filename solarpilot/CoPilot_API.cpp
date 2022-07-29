#include <fstream>
#include <cstring>
#include "rapidxml.hpp"

#include "CoPilot_API.h"
#include "IOUtil.h"
#include "shared/lib_weatherfile.h"
#include "API_structures.h"



struct CopilotObject
{
    SolarField solarfield;
    var_map variables;
    sim_results results;
    LayoutSimThread* simthread;
    SimControl sim_control;
    sp_optical_table opttab;

    std::vector<std::string> message_log;

    int (*external_callback)(sp_number_t fraction_complete, const char* notices);
    bool use_api_callback;

    std::string __str_data;  //used to pass string data back to API language     

    CopilotObject()
    {
        variables.reset();
        solarfield.Create(variables);
        results.clear();
        simthread = 0;
        //solarfield.getSimInfoObject()->setCallbackFunction()
        sim_control.soltrace_callback = ST_APICallback;
        sim_control.soltrace_callback_data = (void*)this;
        sim_control.message_callback = MessageHandler;
        sim_control.message_callback_data = (void*)this;
        sim_control.layout_log_callback = ProgressHandler;
        sim_control.layout_log_callback_data = (void*)this;

        use_api_callback = false;

    };
};

SPEXPORT const char* sp_version(sp_data_t p_data)
{
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;
    
    return V->sf.version.val.c_str();
}

SPEXPORT void sp_set_callback(sp_data_t p_data, int(*fcallback)(sp_number_t, const char*))
{
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    mc->external_callback = fcallback;
    mc->use_api_callback = true;
}

SPEXPORT void sp_disable_callback(sp_data_t p_data)
{
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    mc->use_api_callback = false;
}

SPEXPORT void sp_cancel_simulation(sp_data_t p_data)
{
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    mc->sim_control._cancel_simulation = true;
}

SPEXPORT sp_data_t sp_data_create()
{
    return static_cast<sp_data_t> (new CopilotObject);
}

SPEXPORT bool sp_data_free(sp_data_t p_data)
{
    CopilotObject *mc = static_cast<CopilotObject*>(p_data);

    if (mc)
    {
        delete mc;
        return true;
    }
    return false;
}

SPEXPORT void var_free_memory(sp_number_t* varptr)
{
    delete[] varptr;
};

SPEXPORT bool sp_set_number(sp_data_t p_data, const char* name, sp_number_t v)
{
    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;
    
    //make sure the specified variable exists
    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name) + "\nWARNING: Value was not set!";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    //create a string copy of the variable name
    std::string sname = (std::string)name;
    bool is_combo = (mc->variables._varptrs.at(sname)->ctype == "combo");   //is variable a combo?

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (!(dattype == SP_DATTYPE::SP_DOUBLE || dattype == SP_DATTYPE::SP_INT || dattype == SP_DATTYPE::SP_BOOL || (dattype == SP_DATTYPE::SP_STRING && is_combo)))
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_set_number. \nWARNING: Value was not set!";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    //if it's a combo, set value based on map value
    if (is_combo)
    {
        int cbchoice = (int) v;
        if (mc->variables._varptrs.at(sname)->combo_select_by_mapval(cbchoice))
        {
            // successfully set 
            return true;
        }
        else
        {
            std::string msg = "Invalid variable choice for \"" + sname + "\": \"" + my_to_string(v) + "\" is not a valid option. \nWARNING: Value was not set!";
            SC->message_callback(msg.c_str(), SC->message_callback_data);
            return false;
        }
    }
    else
    {
        // set Boolean values
        std::string svalue = "";
        if (dattype == SP_DATTYPE::SP_BOOL) {
            svalue = (v ? "TRUE" : "FALSE");
        }
        else {
            //no problems, just set the variable
            svalue = std::to_string(v);
        }
        mc->variables._varptrs.at(sname)->set_from_string(svalue.c_str());
        return true;
    }
}

SPEXPORT bool sp_set_string(sp_data_t p_data, const char *name, const char *value)
{
    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    //create a string copy of the variable name
    std::string sname = (std::string)name;
    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name) + "\nWARNING: Value was not set!";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (!(dattype == SP_DATTYPE::SP_STRING))
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_set_string. \nWARNING: Value was not set!";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    //if it's a combo, make sure the specified combo choice exists
    if (mc->variables._varptrs.at(sname)->ctype == "combo")
    {
        std::string svalue = my_to_string(value);

        std::vector< std::string > cbchoices = mc->variables._varptrs.at(sname)->combo_get_choices();
        if (std::find(cbchoices.begin(), cbchoices.end(), svalue) != cbchoices.end())
        {
            //valid variable and selection
            mc->variables._varptrs.at(sname)->set_from_string(svalue.c_str());
            return true;
        }
        else
        {
            std::string msg = "Invalid variable choice for \"" + sname + "\": \"" + my_to_string(value) + "\" is not a valid option. \nWARNING: Value was not set!";
            SC->message_callback(msg.c_str(), SC->message_callback_data);
            return false;
        }
    }
    else
    {
        // no problems, just set the variable
        mc->variables._varptrs.at(sname)->set_from_string(value);
        return true;
    }
}

/** Assigns value of type SP_VEC_DOUBLE and SP_VEC_INTEGER */
SPEXPORT bool sp_set_array(sp_data_t p_data, const char *name, sp_number_t *pvalues, int length)
{
    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    //make sure the specified variable exists
    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name) + "\nWARNING: Value was not set!";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (!(dattype == SP_DATTYPE::SP_VEC_DOUBLE || dattype == SP_DATTYPE::SP_VEC_INTEGER))
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_set_array.  \nWARNING: Value was not set!";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    //create a string copy of the variable name
    std::string sname = (std::string) name;

    //collect the array data into a string with correct format
    std::stringstream array_string;
    
    for (size_t i = 0; i < length; i++)
    {
        if (dattype == SP_DATTYPE::SP_VEC_INTEGER)
            array_string << (int)pvalues[i] << ",";
        else
            array_string << pvalues[i] << ",";
    }

    //assign the array
    mc->variables._varptrs.at(name)->set_from_string(array_string.str().c_str());
    return true;
}

/** Assigns value of type @a SSC_MATRIX . Matrices are specified as a continuous array, in row-major order.  Example: the matrix [[5,2,3],[9,1,4]] is stored as [5,2,3,9,1,4]. */
SPEXPORT bool sp_set_matrix(sp_data_t p_data, const char *name, sp_number_t *pvalues, int nrows, int ncols)
{
    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    //make sure the specified variable exists
    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name) + "\nWARNING: Value was not set!";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    //make sure the data type of the variable provided matches the internal data type
    if (mc->variables._varptrs.at(name)->dattype != SP_DATTYPE::SP_MATRIX_T)
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_set_matrix. \nWARNING: Value was not set!";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    //create a string copy of the variable name
    std::string sname = (std::string)name;

    //collect the array data into a string with correct format
    /*
    rows separated by ';'
    cols separated by ','
    */
    std::stringstream matrix_string;
    for (size_t i = 0; i < nrows; i++)
    {
        for (size_t j = 0; j < ncols; j++)
        {
            matrix_string << pvalues[i*ncols + j] << ",";
        }
        matrix_string << ';';
    }
        
    //assign the matrix
    mc->variables._varptrs.at(name)->set_from_string(matrix_string.str().c_str());
    return true;
}

SPEXPORT sp_number_t sp_get_number(sp_data_t p_data, const char* name)
{
    /*
    Return value of type double, int, or bool. Bools are returned as 0 or 1.
    */
    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name);
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<sp_number_t>::quiet_NaN();
    }

    //create a string copy of the variable name
    std::string sname = (std::string)name;
    bool is_combo = (mc->variables._varptrs.at(sname)->ctype == "combo");   //is variable a combo?

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (!(dattype == SP_DATTYPE::SP_DOUBLE || dattype == SP_DATTYPE::SP_INT || dattype == SP_DATTYPE::SP_BOOL || (dattype == SP_DATTYPE::SP_STRING && is_combo)))
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_get_number.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<sp_number_t>::quiet_NaN();
    }

    spbase *var = mc->variables._varptrs[name];
    switch (var->dattype)
    {
    case SP_INT:
    {
        if (var->is_output)
        {
            spout<int>* v = static_cast<spout<int>*>(var);
            return (sp_number_t)v->Val();
        }
        else
        {
            spvar<int>* v = static_cast<spvar<int>*>(var);
            return (sp_number_t)v->val;
        }
    }
    case SP_DOUBLE:
    {
        if (var->is_output)
        {
            spout<double>* v = static_cast<spout<double>*>(var);
            return (sp_number_t)v->Val();
        }
        else
        {
            spvar<double>* v = static_cast<spvar<double>*>(var);
            return (sp_number_t)v->val;
        }
    }
    case SP_BOOL:
    {
        if (var->is_output)
        {
            spout<bool>* v = static_cast<spout<bool>*>(var);
            return (sp_number_t)(v->Val() ? 1. : 0.);
        }
        else
        {
            spvar<bool>* v = static_cast<spvar<bool>*>(var);
            return (sp_number_t)(v->val ? 1. : 0.);
        }
    }
    case SP_STRING: //only if combo
    {
        //spvar<int>* v = static_cast<spvar<int>*>();
        return (sp_number_t) var->mapval();
    }

    default:
        return std::numeric_limits<sp_number_t>::quiet_NaN();
    }
}

/** Returns the value of a @a SP_STRING variable with the given name. */
SPEXPORT const char *sp_get_string(sp_data_t p_data, const char *name)
{
    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name);
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<const char*>::quiet_NaN();
    }

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (dattype != SP_DATTYPE::SP_STRING)
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_get_string.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<const char*>::quiet_NaN();
    }
    spvar<std::string>* ret = static_cast<spvar<std::string>*>(mc->variables._varptrs.at(name));
    
    return ret->val.c_str();
}

/** Returns the value of a @a SSC_ARRAY variable with the given name. */
SPEXPORT sp_number_t *sp_get_array(sp_data_t p_data, const char *name, int *length)
{
    /*
    Populates 'value_array' with 'length' entries
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name);
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<sp_number_t*>::quiet_NaN();
    }

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (!(dattype == SP_DATTYPE::SP_VEC_DOUBLE || dattype == SP_DATTYPE::SP_VEC_INTEGER))
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_get_array.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<sp_number_t*>::quiet_NaN();
    }

    spbase *var = mc->variables._varptrs[name];
    std::string varstr = var->as_string();

    //convert the string formatted vector to a vector<double>
    std::vector<double> Vd;
    spbase::_setv(varstr, Vd);

    //allocate space at the value_array pointer
    sp_number_t *values = new sp_number_t[(int)Vd.size()];
    //set length for return
    *length = (int)Vd.size();

    //convert to to return format
    for (size_t i = 0; i < *length; i++)
        values[i] = Vd.at(i);

    return values;
}

SPEXPORT int sp_get_array_size(sp_data_t p_data, const char* name)
{
    /*
    Returns the length of array 'name'
    */

    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name);
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<int>::quiet_NaN();
    }

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (!(dattype == SP_DATTYPE::SP_VEC_DOUBLE || dattype == SP_DATTYPE::SP_VEC_INTEGER))
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_get_array_size.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<int>::quiet_NaN();
    }

    spbase* var = mc->variables._varptrs[name];
    std::string varstr = var->as_string();

    //convert the string formatted vector to a vector<double>
    std::vector<double> Vd;
    spbase::_setv(varstr, Vd);

    //set length for return
    return (int)Vd.size();
}

SPEXPORT double sp_get_array_value_by_index(sp_data_t p_data, const char* name, int index)
{
    /*
    Returns array value at 'index' position
    */

    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name);
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<double>::quiet_NaN();
    }

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (!(dattype == SP_DATTYPE::SP_VEC_DOUBLE || dattype == SP_DATTYPE::SP_VEC_INTEGER))
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_get_double_array_value_by_index.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<double>::quiet_NaN();
    }

    spbase* var = mc->variables._varptrs[name];
    std::string varstr = var->as_string();

    //convert the string formatted vector to a vector<double>
    std::vector<double> Vd;
    spbase::_setv(varstr, Vd);

    //test index to be within valid range
    if (index<0 || index > (int) Vd.size()-1)
    {
        std::string msg = "Index " + std::to_string(index) + " is outside valid range of array " + std::string(name) + " use sp_get_array_size to get the valid size.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (dattype == SP_DATTYPE::SP_VEC_INTEGER)
    {
        return (double) Vd.at(index);
    }
    else
    {
        return Vd.at(index);
    }
}

SPEXPORT sp_number_t *sp_get_matrix(sp_data_t p_data, const char *name, int *nrows, int* ncols)
{
    /*
    Populates 'value_array' with 'length' entries
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name);
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<sp_number_t*>::quiet_NaN();
    }

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (dattype != SP_DATTYPE::SP_MATRIX_T)
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_get_matrix.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<sp_number_t*>::quiet_NaN();
    }

    spbase *var = mc->variables._varptrs[name];
    std::string varstr = var->as_string();

    //convert the string formatted vector 
    matrix_t<double> Md;
    spbase::_setv(varstr, Md);

    //allocate space at the value_array pointer
    sp_number_t* values = new sp_number_t[(int)(Md.nrows()*Md.ncols())];
    //set lengths for return
    *ncols = (int) Md.ncols();
    *nrows = (int) Md.nrows();

    //convert to to return format
    for (size_t i = 0; i < *nrows; i++)
        for (size_t j = 0; j < *ncols; j++)
            values[j + (*ncols)*i] = Md.at(i,j);

    return values;
}

SPEXPORT int sp_get_matrix_row_size(void* p_data, const char* name)
{
    /*
    Returns size of matrix. nrows is returned, ncols is returned by reference
    */

    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name);
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<int>::quiet_NaN();
    }

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (dattype != SP_DATTYPE::SP_MATRIX_T)
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_get_matrix.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<int>::quiet_NaN();
    }

    spbase* var = mc->variables._varptrs[name];
    std::string varstr = var->as_string();

    //convert the string formatted vector 
    matrix_t<double> Md;
    spbase::_setv(varstr, Md);

    //set lengths for return
    return (int) Md.nrows();
}

SPEXPORT int sp_get_matrix_col_size(void* p_data, const char* name)
{
    /*
    Returns size of matrix. nrows is returned, ncols is returned by reference
    */

    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name);
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<int>::quiet_NaN();
    }

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (dattype != SP_DATTYPE::SP_MATRIX_T)
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_get_matrix.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<int>::quiet_NaN();
    }

    spbase* var = mc->variables._varptrs[name];
    std::string varstr = var->as_string();

    //convert the string formatted vector 
    matrix_t<double> Md;
    spbase::_setv(varstr, Md);

    //set lengths for return
    return (int)Md.ncols();
}

SPEXPORT double sp_get_matrix_value_by_index(void* p_data, const char* name, int row_index, int col_index)
{
    /*
    Populates 'value_array' with 'length' entries
    */

    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        std::string msg = "No such variable: " + std::string(name);
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<double>::quiet_NaN();
    }

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (dattype != SP_DATTYPE::SP_MATRIX_T)
    {
        std::string msg = "Data type of " + std::string(name) + " is not compatible with sp_get_matrix.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<double>::quiet_NaN();
    }

    spbase* var = mc->variables._varptrs[name];
    std::string varstr = var->as_string();

    //convert the string formatted vector 
    matrix_t<double> Md;
    spbase::_setv(varstr, Md);

    //allocate space at the value_array pointer
    sp_number_t* values = new sp_number_t[(int)(Md.nrows() * Md.ncols())];
    //set lengths for return
    int ncols = (int)Md.ncols();
    int nrows = (int)Md.nrows();

    //test index to be within valid range
    if (row_index<0 || row_index > nrows - 1)
    {
        std::string msg = "Row index " + std::to_string(row_index) + " is outside valid range of matrix " + std::string(name) + " use sp_get_matrix_size to get the valid size.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<double>::quiet_NaN();
    }
    else if (col_index<0 || col_index > ncols - 1)
    {
        std::string msg = "Column index " + std::to_string(col_index) + " is outside valid range of matrix " + std::string(name) + " use sp_get_matrix_size to get the valid size.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return std::numeric_limits<double>::quiet_NaN();
    }

    //return value from matrix
    return Md.at(row_index, col_index);
}

SPEXPORT void sp_reset_geometry(sp_data_t p_data)
{
    /*
	Reset the system geometry to the SolarPILOT default values, clearing any changes or variable settings.
	Returns: (void):null
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    mc->variables.reset();
}

SPEXPORT int sp_add_receiver(sp_data_t p_data, const char* receiver_name)
{
    /*
	Add a new receiver, returning the unique index
	Returns: (string:name[, boolean:make selection]):integer
    */

    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    SimControl* SC = &mc->sim_control;

    std::string tname = std::string(receiver_name);
    var_map* V = &mc->variables;

    // check to make sure this isn't a duplicate. Each item needs a unique name
    bool dupe = false;
    for (unsigned int i = 0; i < V->recs.size(); i++)
    {
        if (tname == V->recs.at(i).rec_name.val)
            dupe = true;
    }
    if (dupe)
    {
        std::string msg = "Receiver name '" + tname + "' is not unique.  Please enter a unique name for this geometry.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return -1;
    }

    //Add a receiver
    int ind = (int) V->recs.size();

    V->add_receiver(ind);
    V->recs[ind].rec_name.val = tname;
    
    //Re-create the solar field object
    mc->solarfield.Create(*V);
    
    //update the calculated parameters
    mc->solarfield.updateAllCalculatedParameters(*V);
    mc->solarfield.updateCalculatedReceiverPower(*V);  // (line 1783) unsure if this is necessary

    return V->recs.back().id.val; 
}

SPEXPORT bool sp_drop_receiver(sp_data_t p_data, const char* receiver_name)
{
    /*
	Drop a receiver from the current solar field
	Returns: (integer:index)
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;
    SimControl* SC = &mc->sim_control;

    std::string tname = lower_case(std::string(receiver_name));

    for (size_t i = 0; i < V->recs.size(); i++)
    {
        if (tname == lower_case(V->recs.at(i).rec_name.val))
        {
            //delete the item
            V->drop_receiver(V->recs.at(i).id.val);
            mc->solarfield.Create(*V);
            return true;
        }
    }

    std::string msg = "Receiver name '" + tname + "' was not found.";
    SC->message_callback(msg.c_str(), SC->message_callback_data);
    return false;
}

SPEXPORT int sp_add_heliostat_template(sp_data_t p_data, const char* heliostat_name)
{
    /*
	Add a new heliostat template that can be used in the layout.
	Returns: (string:template name):integer
    */

    //Add a heliostat
    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;
    SimControl* SC = &mc->sim_control;

    //string tname = cxt.arg(0).as_string();
    std::string tname = std::string(heliostat_name);
    bool dupe = false;
    for (unsigned int i = 0; i < V->hels.size(); i++)
    {
        if (tname == V->hels.at(i).helio_name.val)
            dupe = true;
    }
    if (dupe)
    {
        std::string msg = "Heliostat name '" + tname + "' is not unique.  Please enter a unique name for this heliostat template.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return -1;
    }

    int ind = (int) V->hels.size();
    V->add_heliostat(ind);
    V->hels.back().helio_name.val = tname;
    //Re-create the solar field object
    mc->solarfield.Create(*V);
    // F.GetSolarFieldObject()->Create(*V);

    //cxt.result().assign((double)SPFrame::Instance().GetVariablesObject()->hels.back().id.val);
    return V->hels.back().id.val;
}

SPEXPORT bool sp_drop_heliostat_template(sp_data_t p_data, const char* heliostat_name)
{
    /*
	Delete (drop) the specified heliostat template from the current setup. Returns true if successful.
	Returns: (string:template name):bool
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;
    SimControl* SC = &mc->sim_control;


    std::string tname = lower_case(std::string(heliostat_name));

    for (size_t i = 0; i < V->hels.size(); i++)
    {
        if (tname == lower_case(V->hels.at(i).helio_name.val))
        {
            //delete the item
            V->drop_heliostat(V->hels.at(i).id.val);
            mc->solarfield.Create(*V);
            return true;
        }
    }
    std::string msg = "Heliostat template name '" + tname + "' was not found.";
    SC->message_callback(msg.c_str(), SC->message_callback_data);
    return false;
}


bool _load_weather_file(var_map* V, ArrayString &wf_output)
{
    std::string weatherfile_str = std::string(V->amb.weather_file.val);

    Ambient::readWeatherFile(*V);

    //Saving local verison of weather data
    weatherfile wf;
    if (!wf.open(weatherfile_str))
    {
        //std::string msg = "'update_geometry' function cannot find weather file at " + weatherfile_str + "\n Please adjust desired file path or location to be consistent.";
        //SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false; //error
    }

    //Update the weather data
    std::string linef = "%d,%d,%d,%.2f,%.1f,%.1f,%.1f";
    char cline[300];

    int nrec = (int)wf.nrecords();

    //ArrayString local_wfdat;
    wf_output.resize(nrec);

    weather_record wrec;
    for (int i = 0; i < nrec; i++)
    {
        //int year, month, day, hour;
        wf.read(&wrec);
        sprintf(cline, linef.c_str(), wrec.day, wrec.hour, wrec.month, wrec.dn, wrec.tdry, wrec.pres / 1000., wrec.wspd);
        std::string line(cline);
        wf_output.at(i) = line;
    }

    return true;
}


SPEXPORT sp_number_t* sp_generate_simulation_days(sp_data_t p_data, int *nrecord, int *ncol)
{
    /*
    Generate the simulation days and hours needed to evaluate performance of the field over time. This 
    function is provided as a convenience for generating this information, but is called internally
    by the layout algorithm when needed.
    */
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;
    SolarField* SF = &mc->solarfield;
    SimControl* SC = &mc->sim_control;

    if (SF->getHeliostats()->size() == 0)
    {
        //no layout exists, so we should be calling the 'run_layout' method instead
        std::string msg = "No layout exists, so the 'update_geometry' function cannot be executed. Please first create or import a layout using 'run_layout'.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return 0;
    }

    ArrayString local_wfdat;
    if (!_load_weather_file(V, local_wfdat))
    {
        std::string msg = "'update_geometry' function cannot find weather file.\n Please adjust desired file path or location to be consistent.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false; //error
    }

    interop::GenerateSimulationWeatherData(*V, V->sf.des_sim_detail.mapval(), local_wfdat);

    WeatherData* wd = &V->sf.sim_step_data.Val();
    *nrecord = wd->_N_items;

    sp_number_t* simsteps_out = new sp_number_t[6 * (*nrecord)];

    {
        *ncol = 7; 
        for (int i = 0; i < *nrecord; i++)
        {
            int j = 0;
            simsteps_out[i * *ncol + j++] = wd->Month.at(i);
            simsteps_out[i * *ncol + j++] = wd->Day.at(i);
            simsteps_out[i * *ncol + j++] = wd->Hour.at(i);
            simsteps_out[i * *ncol + j++] = wd->DNI.at(i);
            simsteps_out[i * *ncol + j++] = wd->T_db.at(i);
            simsteps_out[i * *ncol + j++] = wd->V_wind.at(i);
            simsteps_out[i * *ncol + j++] = wd->Step_weight.at(i);
        }
    }

    return simsteps_out;
}


SPEXPORT bool sp_update_geometry(sp_data_t p_data)
{
    /*
	Refresh the solar field, receiver, or ambient condition settings based on the current parameter settings.
	Returns: (void):boolean
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;
    SolarField* SF = &mc->solarfield;
    SimControl* SC = &mc->sim_control;
    
    if (SF->getHeliostats()->size() == 0)
    {
        //no layout exists, so we should be calling the 'run_layout' method instead
        std::string msg = "No layout exists, so the 'update_geometry' function cannot be executed. Please first create or import a layout using 'run_layout'.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    //std::string weatherfile_str = std::string(V->amb.weather_file.val);
    ArrayString local_wfdat;

    //Ambient::readWeatherFile(*V);
    //
    ////Saving local verison of weather data
    //weatherfile wf;
    //if (!wf.open(weatherfile_str))
    if(! _load_weather_file(V, local_wfdat))
    {
        std::string msg = "'update_geometry' function cannot find weather file.\n Please adjust desired file path or location to be consistent.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false; //error
    }

    ////Update the weather data
    //std::string linef = "%d,%d,%d,%.2f,%.1f,%.1f,%.1f";
    //char cline[300];

    //int nrec = (int)wf.nrecords();

    //local_wfdat.resize(nrec);

    //weather_record wrec;
    //for (int i = 0; i < nrec; i++)
    //{
    //    //int year, month, day, hour;
    //    wf.read(&wrec);
    //    sprintf(cline, linef.c_str(), wrec.day, wrec.hour, wrec.month, wrec.dn, wrec.tdry, wrec.pres / 1000., wrec.wspd);
    //    std::string line(cline);
    //    local_wfdat.at(i) = line;
    //}

    //F.UpdateDesignSelect(V->sf.des_sim_detail.mapval(), *V);
        // Function seems to only update var_map with simulation data through GenearateSimulationWeatherData()
    interop::GenerateSimulationWeatherData(*V, V->sf.des_sim_detail.mapval(), local_wfdat);

    //Set up the solar field
    SF->Clean();
    SF->Create(*V);

    try
    {
        SolarField::PrepareFieldLayout(mc->solarfield, 0, true);

        if (mc->solarfield.ErrCheck())
        {
            std::string msg = "An error occurred when preparing the updated field geometry in the call 'update_geometry'.";
            SC->message_callback(msg.c_str(), SC->message_callback_data);
            return false;
        }

        SF->calcHeliostatArea();
        SF->updateAllCalculatedParameters(*V);

        double azzen[2];
        mc->solarfield.CalcDesignPtSunPosition(V->sf.sun_loc_des.mapval(), azzen[0], azzen[1]);
        Vect sun = Ambient::calcSunVectorFromAzZen(azzen[0] * D2R, azzen[1] * D2R);

        SF->updateAllTrackVectors(sun);

        if (SF->ErrCheck())
        {
            std::string msg = "An error occurred when preparing the updated field geometry in the call 'update_geometry'.";
            SC->message_callback(msg.c_str(), SC->message_callback_data);
            return false;
        }
    }
    catch (std::exception &e)
    {
        std::string msg = "An error occurred when preparing the updated field geometry in the call 'update_geometry'. \n ERROR: "+ (std::string)e.what();
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }
    catch (...)
    {
        std::string msg = "Unknown error when executing 'update_geometry'.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    return true;
}

SPEXPORT bool sp_generate_layout(sp_data_t p_data, int nthreads = 0) //, bool save_detail = true)
{
    /*
    Create a solar field layout. Options include 'nthreads':integer (default All),'save_detail':boolean (default True)",
    run layout without specified positions SolarPILOT generates layout positions ([table:options])
    Returns: boolean
        
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SolarField* SF = &mc->solarfield;
    var_map* V = &mc->variables;

    SimControl* SC = &mc->sim_control;
    LayoutSimThread* SThread = mc->simthread;

    if (nthreads!=0)
        SC->SetThreadCount(nthreads);    

    // Set heliostat template if only one exists - 
    if ((int)SF->getHeliostatTemplates()->size() == 1) 
    {
        if ((int)V->sf.temp_which.combo_get_count() == 0)
        {
            // Add  
            std::string js = my_to_string(0);
            V->sf.temp_which.combo_add_choice(V->hels.at(0).helio_name.val, js);
        }
        V->sf.temp_which.combo_select_by_choice_index(0);
    } // if single template is selected but there is more than one template available
    else if ( ((int)V->sf.template_rule.combo_get_current_index() == var_solarfield::TEMPLATE_RULE::USE_SINGLE_TEMPLATE)
                && (V->sf.temp_which.as_string() == "") )
    {
        std::string msg = "ERROR: Please set solarfield 'temp_which' value with the heliostat template name to use for field generation.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    std::string weatherfile_str = std::string(V->amb.weather_file.val);
    Ambient::readWeatherFile(*V);

    //Saving local verison of weather data
    weatherfile wf;
    if (!wf.open(weatherfile_str))
    {
        std::string msg = "ERROR: Failed to open weather file.  Please check provided file path.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
    	return false; //error
    }

    //Update the weather data
    std::string linef = "%d,%d,%d,%.2f,%.1f,%.1f,%.1f";
    char cline[300];

    int nrec = (int)wf.nrecords();

    ArrayString local_wfdat;
    local_wfdat.resize(nrec);

    weather_record wrec;
    for (int i = 0; i < nrec; i++)
    {
        //int year, month, day, hour;
        wf.read(&wrec);
        sprintf(cline, linef.c_str(), wrec.day, wrec.hour, wrec.month, wrec.dn, wrec.tdry, wrec.pres / 1000., wrec.wspd);
        std::string line(cline);
        local_wfdat.at(i) = line;
    }

    interop::GenerateSimulationWeatherData(*V, V->sf.des_sim_detail.mapval(), local_wfdat);


    SF->Clean();
    SF->Create(*V);
    bool simok = interop::DoManagedLayout(*SC, *SF, *V, SThread);        //Returns TRUE if successful

    return simok;
}

SPEXPORT bool sp_assign_layout(sp_data_t p_data, sp_number_t* pvalues, int nrows, int ncols, int nthreads = 0) //, bool save_detail = true)
{
    /*
    Run layout with specified positions. User specifies layout positions in the following format, where first 4 columns are required:
        "<template (int)> <location X> <location Y> <location Z> <x focal length> <y focal length> <cant i> <cant j> <cant k> <aim X> <aim Y> <aim Z>" (array:positions)
    Returns: boolean
    */
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;

    //assign layout positions
    V->sf.layout_data.val.clear();

    std::stringstream heliodata;

    size_t npos = nrows;
    for (size_t i = 0; i < npos; i++)
    {
        for (size_t j = 0; j < 12; j++)
        {
            if (j > ncols - 1)
                heliodata << "NULL";
            else
                heliodata << pvalues[j + ncols * i];
            heliodata << (j < 11 ? "," : ";");
        }

    }
    V->sf.layout_data.val = heliodata.str();

    //user specified layout
    V->sf.layout_method.combo_select_by_mapval(var_solarfield::LAYOUT_METHOD::USERDEFINED);

    bool simok = sp_generate_layout(p_data, nthreads);

    return simok;
}

SPEXPORT sp_number_t* sp_get_layout_info(sp_data_t p_data, int* nhelio, int* ncol, bool get_corners = false)
{
    /*
    Get information regarding the heliostat field layout. Returns matrix with each row corresponding to a heliostat.
        "Information includes: [index, position-x, position-y, position-z, template_id, ranking metric value]...[corner positions x,y,z]
    Returns: (void):table
    */
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    //var_map* V = &mc->variables;
    SimControl* SC = &mc->sim_control;

    SolarField* SF = &mc->solarfield;
    Hvector* hels = SF->getHeliostats();

    htemp_map htemps = *SF->getHeliostatTemplates();

    *nhelio = (int) hels->size();
    *ncol = 6;
    if (get_corners)
        *ncol += (int) hels->at(0)->getCornerCoords()->size() * 3;  // assumes all heliostats have a equal number of corners 

    sp_number_t* layoutinfo = new sp_number_t[(*nhelio) * (*ncol)];

    int c;
    for (size_t i = 0; i < (int)hels->size(); i++)
    {
        Heliostat* hel = hels->at(i);
        sp_point* loc = hel->getLocation();

        c = 0;  layoutinfo[i * (*ncol) + c] = hel->getId();
        c++;  layoutinfo[i * (*ncol) + c] = loc->x;
        c++;  layoutinfo[i * (*ncol) + c] = loc->y;
        c++;  layoutinfo[i * (*ncol) + c] = loc->z;
        c++;  layoutinfo[i * (*ncol) + c] = hel->getMasterTemplate()->getVarMap()->id.val;
        c++;  layoutinfo[i * (*ncol) + c] = hel->getRankingMetricValue();
        if (get_corners)
        {
            std::vector<sp_point>* corners = hel->getCornerCoords();
            for (size_t j = 0; j < (int)corners->size(); j++)
            {
                sp_point corner = corners->at(j);
                c++; layoutinfo[i * (*ncol) + c] = corner.x;
                c++; layoutinfo[i * (*ncol) + c] = corner.y;
                c++; layoutinfo[i * (*ncol) + c] = corner.z;
            }

        }

        if ((i == 0) && (c != *ncol - 1))
        {
            std::string msg = "Information was lost check sp_get_layout_info output formating.";
            SC->message_callback(msg.c_str(), SC->message_callback_data);
            delete[] layoutinfo;
            return nullptr;
        }
    }

    return layoutinfo;
}

SPEXPORT const char* sp_get_layout_header(sp_data_t p_data, bool get_corners = false)
{
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    std::string tab_header;

    tab_header.append("id,");
    tab_header.append("x_location,");
    tab_header.append("y_location,");
    tab_header.append("z_location,");
    tab_header.append("helio_template_id,");
    tab_header.append("layout_metric");

    if (get_corners)
    {
        SolarField* SF = &mc->solarfield;
        Heliostat* hel = SF->getHeliostats()->at(0);
        std::vector<std::string > coords{ "x", "y", "z" };
        for (size_t j = 0; j < (int)hel->getCornerCoords()->size(); j++)
        {
            for (size_t i = 0; i < coords.size(); i++)
            {
                tab_header.append(",corner" + std::to_string(j) + "_" + coords.at(i));
            }
        }
    }

    mc->__str_data.clear();
    mc->__str_data = tab_header;

    return mc->__str_data.c_str();
}

SPEXPORT bool sp_simulate(sp_data_t p_data, int nthreads = 1, bool update_aimpoints = true) //bool save_detail = true,
//SPEXPORT void sp_simulate(sp_data_t p_data)
{
    /*
    Calculate heliostat field performance. Options include 'nthreads':integer (default All),
    'save_detail':boolean (default True), 'update_aimpoints':boolean (default True)
    Returns: [table:options]):boolean
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SolarField* SF = &mc->solarfield;
    var_map* V = &mc->variables;
    sim_results* res = &mc->results;
    SimControl* SC = &mc->sim_control;
    
    if (nthreads != 1)
        SC->SetThreadCount(nthreads);
    
    if (update_aimpoints == false)
        V->flux.aim_method.combo_select("Keep existing");

    //Which type of simulation is this?
    int simtype = V->flux.flux_model.mapval();    //0=Delsol, 1=Soltrace

    //Set up field, update aimpoints, and simulate at the performance sun position
    Hvector* helios = SF->getHeliostats();

    if (!interop::PerformanceSimulationPrep(*SF, *helios, simtype)) return false;

    Vect sun = Ambient::calcSunVectorFromAzZen(V->flux.flux_solar_az.Val() * D2R, (90. - V->flux.flux_solar_el.Val())*D2R);

    SF->calcHeliostatShadows(sun);
    if (SF->ErrCheck()) return false;

    res->clear();
    res->resize(1);

    // start timer
    std::clock_t start;
    double duration;

    start = std::clock();

    //Which type of simulation?
    bool simok;
    switch (simtype)
    {
    case var_fluxsim::FLUX_MODEL::HERMITE_ANALYTICAL:
        simok = interop::HermiteFluxSimulationHandler(*res, *SF, *helios);
        break;
    case var_fluxsim::FLUX_MODEL::SOLTRACE:
        simok = interop::SolTraceFluxSimulation(*SC, *res, *SF, *V, *helios);
        break;
    default:
        return false;
    }

    //End timer
    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    std::string msg = "Simulation total time:" + std::to_string(duration) + " seconds";
    SC->message_callback(msg.c_str(), SC->message_callback_data);

    SF->getSimInfoObject()->Reset();

    //F.GetFluxPlotObject()->SetPlotData(*SF, *helios, 0);
    //F.GetFieldPlotObject()->SetPlotData(*SF, FIELD_PLOT::EFF_TOT);

    return simok;
}

SPEXPORT const char *sp_summary_results(sp_data_t p_data)
{
    /*
	Return an array of tables with summary results from each simulation. The array length is greater than 1 for multiple-receiver simulations.
	Returns: (void):array
    */
    
    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;

    grid_emulator_base table;
    sim_results* results = &mc->results;

    std::string ret;    //return string

    if (results->size() < 1)
    {
        ret = "No simulation summary results exist. Please simulate first.";
        mc->sim_control.message_callback(ret.c_str(), mc->sim_control.message_callback_data);
        ret = "Failure";
        return ret.c_str();
    }
    
    // for multiple receivers
    for (size_t i = 0; i < V->recs.size(); i++)
    {
        interop::CreateResultsTable(results->at(i), table);
        
        unordered_map<std::string, double> res_map;

        for (int j = 0; j < table.GetNumberRows(); j++)
        {
            ret.append(table.GetRowLabelValue(j) + ", " + table.GetCellValue(j, 1) + "\n");
            res_map.insert({ table.GetRowLabelValue(j), std::stod(table.GetCellValue(j, 1)) });
        }

        //add a few more summary results
        bool is_soltrace = res_map.find("Shadowing and Cosine efficiency") != res_map.end();

        double Qwf;
        double Qin = Qwf = res_map.at("Power incident on field");

        if (is_soltrace)
        {
            /*
            soltrace
            for this option, the "Shadowing and Cosine efficiency" is already calculated by the
            program. Just make sure the Shading and Cosine efficiencies aren't double counted.
            */
            ret.append("Shading efficiency, 100.\n");
            ret.append("Cosine efficiency, 100.\n");
            ret.append("Shading loss, 0.\n");
            ret.append("Cosine loss, 0.\n");

            double eta_sc = res_map.at("Shadowing and Cosine efficiency") / 100.;
            Qwf *= eta_sc;
        }
        else
        {
            //hermite
            double eta_sc = res_map.at("Shading efficiency") * res_map.at("Cosine efficiency") / 100.;
            ret.append("Shadowing and Cosine efficiency, " + std::to_string(eta_sc) + "\n");

            double eta_s = res_map.at("Shading efficiency") / 100.;
            Qwf *= eta_s;
            ret.append("Shading loss, " + std::to_string(Qin*(1. - eta_s)) + "\n");
            double eta_c = res_map.at("Cosine efficiency") / 100.;
            ret.append("Cosine loss, " + std::to_string(Qwf*(1 - eta_c)) + "\n");
            Qwf *= eta_c;
        }
        ret.append("Shadowing and Cosine loss, " + std::to_string(Qin - Qwf) + "\n");

        double eta_r = res_map.at("Reflection efficiency") / 100.;
        ret.append("Reflection loss, " + std::to_string(Qwf * (1. - eta_r)) + "\n");
        Qwf *= eta_r;
        double eta_b = res_map.at("Blocking efficiency") / 100.;
        ret.append("Blocking loss, " + std::to_string(Qwf*(1. - eta_b)) + "\n");
        Qwf *= eta_b;
        double eta_i = res_map.at("Image intercept efficiency") / 100.;
        ret.append("Image intercept loss, " + std::to_string(Qwf*(1. - eta_i)) + "\n");
        Qwf *= eta_i;
        double eta_a = res_map.at("Absorption efficiency") / 100.;
        ret.append("Absorption loss, " + std::to_string(Qwf*(1. - eta_a)) + "\n");

        ret.append("Receiver name, " + (i == 0 ? "All receivers" : results->at(i).receiver_names.front() ) );
    }
    mc->__str_data.clear();
    mc->__str_data = ret;

    return mc->__str_data.c_str();
}

SPEXPORT double sp_get_receiver_power(sp_data_t p_data)
{
    /*
    Return an array of tables with summary results from each simulation. The array length is greater than 1 for multiple-receiver simulations.
    Returns: (void):array
    */

    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;

    grid_emulator_base table;
    sim_results* results = &mc->results;

    if (results->size() < 1)
    {
        return 0.0;
    }

    // for multiple receivers
    for (size_t i = 0; i < V->recs.size(); i++)
    {
        interop::CreateResultsTable(results->at(i), table);

        unordered_map<std::string, double> res_map;

        for (int j = 0; j < table.GetNumberRows(); j++)
        {
            res_map.insert({ table.GetRowLabelValue(j), std::stod(table.GetCellValue(j, 1)) });
        }

        return res_map.at("Power absorbed by the receiver");
    }
}

SPEXPORT sp_number_t* sp_detail_results(sp_data_t p_data, int* nrows, int* ncols, sp_number_t* selhel = NULL, int nselhel = 0, bool get_corners = false)
{
    /*
    returns a vector with hash entries for each heliostat

    id(integer), location (array), aimpoint (array), tracking_vector (array), layout_metric (double), power_to_receiver (double),
        efficiency (double), cosine (double), intercept (double), reflectance (double), attenuation (double),
        blocking (double), shading (double), clouds (double)


        "[Only valid for Hermite (analytical) simulation engine.]\n"
        "Return an array with detailed heliostat-by-heliostat results from a simulation. "
        "Each entry in the array is a table with entries as follows:\n"
        "{ id(integer), location (array), aimpoint (array), tracking_vector (array), "
        "layout_metric (double), "
        "power_to_receiver (double), "
        "power_reflected (double), "
        "energy (double),"
        "annual efficiency (double),"
        "total efficiency (double), "
        "cosine (double), "
        "intercept (double), "
        "reflectance (double), "
        "attenuation (double), "
        "blocking (double), "
        "shading (double), "
        "clouds (double) }",
        "([array:selected heliostat indices]):array");
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SolarField* SF = &mc->solarfield;
    SimControl* SC = &mc->sim_control;

    if (SF->getHeliostats()->size() > 0)
    {
        Hvector helio_select;

        if (nselhel != 0) //use selected heliostats
        {
            std::vector<int> ids;
            ids.reserve(nselhel);

            for (size_t i = 0; i < nselhel; i++)
                ids.push_back((int)selhel[i]);

            helio_select.reserve(nselhel);

            unordered_map<int, Heliostat*> *heliosid = SF->getHeliostatsByID();

            for (size_t i = 0; i < ids.size(); i++)
                helio_select.push_back(heliosid->at(ids.at(i)));

        }
        else    //use all heliostats
        {
            helio_select.assign(SF->getHeliostats()->begin(), SF->getHeliostats()->end());
        }

        *nrows = (int)helio_select.size();
        *ncols = 23;  //number of results in table
        if (get_corners)
            *ncols += (int)helio_select.at(0)->getCornerCoords()->size() * 3;  // assumes all heliostats have a equal number of corners
        if (SF->getReceivers()->size() > 1)
            *ncols += 1;

        //loop through selected heliostats, gathering information
        sp_number_t* ret = new sp_number_t[(*nrows) * (*ncols)];

        int c;
        for (size_t i = 0; i < helio_select.size(); i++)
        {
            //select specific heliostat object
            Heliostat* H = helio_select.at(i);

            // NOTE: If the order of this is change or if data is added or removed, then
            //              update ncols (above) and sp_detail_results_header() (below)
            c = 0; ret[i*(*ncols) + c] = H->getId();
            c++; ret[i * (*ncols) + c] = H->getLocation()->x;
            c++; ret[i * (*ncols) + c] = H->getLocation()->y;
            c++; ret[i * (*ncols) + c] = H->getLocation()->z;
            
            c++; ret[i * (*ncols) + c] = H->getAimPoint()->x;
            c++; ret[i * (*ncols) + c] = H->getAimPoint()->y;
            c++; ret[i * (*ncols) + c] = H->getAimPoint()->z;
            
            c++; ret[i * (*ncols) + c] = H->getTrackVector()->i;
            c++; ret[i * (*ncols) + c] = H->getTrackVector()->j;
            c++; ret[i * (*ncols) + c] = H->getTrackVector()->k;

            if (get_corners) // adding heliostat corner coordinates
            {
                std::vector<sp_point>* corners = H->getCornerCoords();
                for (size_t j = 0; j < (int)corners->size(); j++)
                {
                    sp_point corner = corners->at(j);
                    c++; ret[i * (*ncols) + c] = corner.x;
                    c++; ret[i * (*ncols) + c] = corner.y;
                    c++; ret[i * (*ncols) + c] = corner.z;
                }

            }
            // finding receiver the heliostat is pointed at (for multi-receiver fields)
            if (SF->getReceivers()->size() > 1)
            {
                int r = 0;
                for (r = 0; r < (int)SF->getReceivers()->size(); r++)
                {
                    if (H->getWhichReceiver() == SF->getReceivers()->at(r))
                    {
                        break;
                    }
                }
                c++; ret[i * (*ncols) + c] = r;
            }

            c++; ret[i * (*ncols) + c] = H->getRankingMetricValue();
            c++; ret[i * (*ncols) + c] = H->getPowerToReceiver() / 1000.;  //kW
            c++; ret[i * (*ncols) + c] = H->getArea()
                *H->getEfficiencyCosine()
                *H->getTotalReflectivity()
                *H->getEfficiencyBlock()
                *H->getEfficiencyShading()
                *H->getEfficiencyCloudiness()
                *SF->getVarMap()->flux.flux_dni.val / 1000.;  //kW
            c++; ret[i * (*ncols) + c] = H->getEnergyValue(); //kWh -- energy delivered over the simulation time period
            c++; ret[i * (*ncols) + c] = H->getAnnualEfficiency();
            c++; ret[i * (*ncols) + c] = H->getEfficiencyTotal();
            c++; ret[i * (*ncols) + c] = H->getEfficiencyCosine();
            c++; ret[i * (*ncols) + c] = H->getEfficiencyIntercept();
            c++; ret[i * (*ncols) + c] = H->getTotalReflectivity();
            c++; ret[i * (*ncols) + c] = H->getEfficiencyAtten();
            c++; ret[i * (*ncols) + c] = H->getEfficiencyBlock();
            c++; ret[i * (*ncols) + c] = H->getEfficiencyShading();
            c++; ret[i * (*ncols) + c] = H->getEfficiencyCloudiness();

            if ((i == 0) && (c != *ncols - 1))
            {
                SC->message_callback("Information was lost check sp_detail_results output formating.", SC->message_callback_data);
                delete[] ret;
                return nullptr;
            }
           
        }

        return ret;
    }
    else
    {
        SC->message_callback("Solarfield object empty... be sure to successfully generate and simulate field before calling sp_detailed_results.", SC->message_callback_data);
        return nullptr;
    }

    return nullptr;
}

SPEXPORT const char* sp_detail_results_header(sp_data_t p_data, bool get_corners = false)
{
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    SolarField* SF = &mc->solarfield;

    if (SF->getHeliostats()->size() > 0)
    {

        std::string tab_header;

        // UPDATE: If table results change
        tab_header.append("id,");
        tab_header.append("x_location,");
        tab_header.append("y_location,");
        tab_header.append("z_location,");
        tab_header.append("x_aimpoint,");
        tab_header.append("y_aimpoint,");
        tab_header.append("z_aimpoint,");
        tab_header.append("i_tracking_vector,");
        tab_header.append("j_tracking_vector,");
        tab_header.append("k_tracking_vector,");

        if (get_corners)
        {
            Heliostat* hel = SF->getHeliostats()->at(0);
            std::vector<std::string > coords{ "x", "y", "z" };
            for (size_t j = 0; j < (int)hel->getCornerCoords()->size(); j++)
            {
                for (size_t i = 0; i < coords.size(); i++)
                {
                    tab_header.append("corner" + std::to_string(j) + "_" + coords.at(i) + ",");
                }
            }
        }

        if (SF->getReceivers()->size() > 1)
            tab_header.append("receiver_map,");

        tab_header.append("layout_metric,");
        tab_header.append("power_to_receiver,");
        tab_header.append("power_reflected,");
        tab_header.append("energy,");
        tab_header.append("efficiency_annual,");
        tab_header.append("efficiency,");
        tab_header.append("cosine,");
        tab_header.append("intercept,");
        tab_header.append("reflectance,");
        tab_header.append("attenuation,");
        tab_header.append("blocking,");
        tab_header.append("shading,");
        tab_header.append("clouds");

        mc->__str_data.clear();
        mc->__str_data = tab_header;

        return mc->__str_data.c_str();
    }

    return nullptr;
}

SPEXPORT sp_number_t* sp_get_fluxmap(sp_data_t p_data, int* nrows, int* ncols, int rec_id = 0)
{
    /*
	Retrieve the receiver fluxmap, optionally specifying the receiver ID to retrieve.
	Returns: ([integer:receiver id]):array
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SolarField* SF = &mc->solarfield;

    Receiver *rec;

    if (rec_id != 0)
    {
        if (rec_id > SF->getReceivers()->size() - 1)
            return nullptr;

        rec = SF->getReceivers()->at(rec_id);
    }
    else
    {
        rec = SF->getReceivers()->front();
    }

    FluxGrid *fg = rec->getFluxSurfaces()->front().getFluxMap();

    *nrows = (int)fg->front().size();
    *ncols = (int)fg->size();
    sp_number_t* fluxmap = new sp_number_t[(*nrows) * (*ncols)];

    //rows
    for (size_t i = 0; i < fg->front().size(); i++)
    {
        //cols
        for (size_t j = 0; j < fg->size(); j++)
        {
            fluxmap[i * (*ncols) + j] = fg->at(fg->size()-1-j).at(i).flux;
        }
    }
    return fluxmap;
}

//TODO: Skipped this function initially 
SPEXPORT void sp_optimize(sp_data_t p_data, sp_number_t* pvalues, int nvar)
{
    /*
    Execute an optimization run, returning the optimized result and iteration information. 
    Variables to be optimized are passed in a vector, with each row containing a table specifying 
    {variable, step, upbound, lowbound, inital}. The table must include the variable key, others are optional. 
    
    The return table includes the following: 'result':table of variable names and associated optimized values, 
    'objective':number, 'flux':number, 'iterations':array of evaluation point, objective, flux. 
    Optional arguments include maxiterations/tolerance/defaultstep/powerpenalty/nthreads.,
    (vector:variable tables[, table:options])
    Returns: table
    */

    //api_helper *mc = static_cast<api_helper*>(p_data);
    //var_map* V = &mc->variables;
    //
    ////get the variable table
    //if (cxt.arg_count() < 1 || cxt.arg(0).type() != lk::vardata_t::VECTOR)
    //    return;
    //
    //lk::vardata_t &vartab = cxt.arg(0);
    //
    //std::stringstream heliodata;
    //
    //size_t nvars = nvar;
    //for (size_t i = 0; i < nvars; i++)
    //{
    //    for (size_t j = 0; j < 12; j++)
    //    {
    //        if (j > ncols - 1)
    //            heliodata << "NULL";
    //        else
    //            heliodata << pvalues[j + ncols * i];
    //        heliodata << (j < 11 ? "," : ";");
    //    }
    //
    //}
    //V->sf.layout_data.val = heliodata.str();
    //
    //
    ////set up options, if provided
    //if (cxt.arg_count() == 2)
    //{
    //    //maxiterations/tolerance/defaultstep/powerpenalty/nthreads
    //
    //    lk::varhash_t *opthash = cxt.arg(1).hash();
    //
    //    if (opthash->find("nthreads") != opthash->end())
    //        F.SetThreadCount(opthash->at("nthreads")->as_integer());
    //
    //    if (opthash->find("maxiterations") != opthash->end())
    //        V->opt.max_iter.val = opthash->at("maxiterations")->as_integer();
    //
    //    if (opthash->find("tolerance") != opthash->end())
    //        V->opt.converge_tol.val = opthash->at("tolerance")->as_number();
    //
    //    if (opthash->find("defaultstep") != opthash->end())
    //        V->opt.max_step.val = opthash->at("defaultstep")->as_number();
    //
    //    if (opthash->find("powerpenalty") != opthash->end())
    //        V->opt.power_penalty.val = opthash->at("powerpenalty")->as_number();
    //}
    //
    ////set up variables
    //int nv = vartab.vec()->size();
    //vector<double*> optvars(nv);
    //vector<double> upper(nv);
    //vector<double> lower(nv);
    //vector<double> stepsize(nv);
    //vector<string> names(nv);
    //
    //for (size_t i = 0; i < nv; i++)
    //{
    //    //check that the specified variable names exist
    //    std::string varname = vartab.vec()->at(i).hash()->at("variable")->as_string();
    //
    //    if (V->_varptrs.find(varname) == V->_varptrs.end())
    //        throw lk::error_t("Specified variable does not exist: " + varname);
    //
    //    //handle the variable
    //    if (V->_varptrs.at(varname)->dattype != SP_DATTYPE::SP_DOUBLE)
    //        throw lk::error_t("Optimized variable must be of type 'double'; discrete or boolean variables are not supported. Variable: " + varname);
    //
    //    spvar<double> *varptr = static_cast<spvar<double>*>(V->_varptrs.at(varname));
    //    optvars.at(i) = &varptr->val;
    //    vector<string> namedat = split(varname, ".");
    //    names.at(i) = namedat.back();
    //
    //    lk::varhash_t *varhash = vartab.vec()->at(i).hash();
    //
    //    //bounds
    //    if (varhash->find("lowbound") == varhash->end())
    //        lower.at(i) = -HUGE_VAL;
    //    else
    //        lower.at(i) = varhash->at("lowbound")->as_number();
    //
    //    if (varhash->find("upbound") == varhash->end())
    //        upper.at(i) = HUGE_VAL;
    //    else
    //        upper.at(i) = varhash->at("upbound")->as_number();
    //
    //    if (varhash->find("initial") != varhash->end())
    //        varptr->val = varhash->at("initial")->as_number();
    //
    //    if (varhash->find("step") == varhash->end())
    //        stepsize.at(i) = V->opt.max_step.val * varptr->val;
    //    else
    //        stepsize.at(i) = varhash->at("step")->as_number();
    //
    //}
    //
    //int n_threads = F.GetThreadCount();
    //ArrayString *local_wfdat = F.GetLocalWeatherDataObject();
    //lk::vardata_t iter_vec;
    //std::vector< double > obj_vals;
    //std::vector< std::vector<double> > flux_vals;
    //std::vector< std::vector< double > > eval_points;
    //
    //if (n_threads > 1)
    //{
    //    AutoPilot_MT *SFopt_MT = new AutoPilot_MT();
    //
    //    SFopt_MT->SetSummaryCallback(LKInfoCallback, SF->getSimInfoObject()->getCallbackData());
    //
    //    //set up the weather data for simulation
    //    vector<string> wdata;
    //    for (int i = 0; i < local_wfdat->size(); i++)
    //        wdata.push_back(local_wfdat->at(i));
    //    SFopt_MT->GenerateDesignPointSimulations(*V, wdata);
    //
    //    //Do the expert setup
    //    SFopt_MT->Setup(*V, true);
    //
    //    //run the optimization
    //    SFopt_MT->Optimize(optvars, upper, lower, stepsize, &names);
    //
    //    //get resulting info
    //    SFopt_MT->GetOptimizationObject()->getOptimizationSimulationHistory(eval_points, obj_vals, flux_vals);
    //
    //    try
    //    {
    //        delete SFopt_MT;
    //    }
    //    catch (...)
    //    {
    //    }
    //}
    //else
    //{
    //
    //    AutoPilot_S *SFopt_S = new AutoPilot_S();
    //    SFopt_S->SetSummaryCallback(LKInfoCallback, SF->getSimInfoObject()->getCallbackData());
    //
    //    //set up the weather data for simulation
    //    vector<string> wdata;
    //    for (int i = 0; i < local_wfdat->size(); i++)
    //        wdata.push_back(local_wfdat->at(i));
    //    SFopt_S->GenerateDesignPointSimulations(*V, wdata);
    //
    //    //Do the expert setup
    //    SFopt_S->Setup(*V, true);
    //
    //    //run the optimization
    //    SFopt_S->Optimize(optvars, upper, lower, stepsize, &names);
    //
    //    //get resulting info
    //    SFopt_S->GetOptimizationObject()->getOptimizationSimulationHistory(eval_points, obj_vals, flux_vals);
    //
    //
    //    try
    //    {
    //        delete SFopt_S;
    //    }
    //    catch (...)
    //    {
    //    }
    //}
    //
    //
    ////set up return structure
    ////result/objective/flux/iterations
    //cxt.result().empty_hash();
    //
    //lk::vardata_t res_hash;
    //res_hash.empty_hash();
    //
    //for (size_t i = 0; i < optvars.size(); i++)
    //{
    //    std::string varname = vartab.vec()->at(i).hash()->at("variable")->as_string();
    //    spvar<double> *varptr = static_cast<spvar<double>*>(V->_varptrs.at(varname));
    //    res_hash.hash_item(varname, varptr->val);
    //}
    //
    //cxt.result().hash_item("result", res_hash);
    //
    //iter_vec.empty_vector();
    //
    //for (size_t i = 0; i < flux_vals.size(); i++)
    //{
    //    iter_vec.vec()->push_back(lk::vardata_t());
    //    iter_vec.vec()->at(i).empty_vector();
    //    for (size_t j = 0; j < eval_points.front().size(); j++)
    //    {
    //        iter_vec.vec()->at(i).vec_append(eval_points.at(i).at(j));
    //    }
    //    iter_vec.vec()->at(i).vec_append(obj_vals.at(i));
    //    for (size_t j = 0; j < flux_vals.front().size(); j++)
    //        iter_vec.vec()->at(i).vec_append(flux_vals.at(i).at(j));
    //}
    //
    //lk::vardata_t fluxresult;
    //fluxresult.empty_vector();
    //for (size_t j = 0; j < flux_vals.back().size(); j++)
    //    fluxresult.vec_append(flux_vals.back().at(j));
    //
    //cxt.result().hash_item("objective", obj_vals.back());
    //cxt.result().hash_item("flux", fluxresult);
    //cxt.result().hash_item("iterations", iter_vec);
    //
}

SPEXPORT void sp_clear_land(sp_data_t p_data, const char* type = NULL)
{
    /*
	Reset the land boundary polygons, clearing any data. Optionally specify 'type' as 'inclusion' or 'exclusion'.
	Returns: ([string:type]):void
    */
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;

    bool clear_inclusions = true;
    bool clear_exclusions = true;

    if (type != NULL)
    {
        std::string arg = lower_case(type);

        if (arg.find("inclusion") != std::string::npos)
            clear_exclusions = false;
        else if (arg.find("exclusion") != std::string::npos)
            clear_inclusions = false;
    }

    if (clear_inclusions)
        V->land.inclusions.val.clear();
    if (clear_exclusions)
        V->land.exclusions.val.clear();

    if (clear_inclusions && clear_exclusions)
        V->land.is_bounds_array.val = false;
}

SPEXPORT bool sp_add_land(sp_data_t p_data, const char* type, sp_number_t* polygon_points, int* npts , int* ndim, bool is_append = true)
{
    /*
	Add land inclusion or a land exclusion region within a specified polygon. Specify the type as 'inclusion' or 'exclusion', and optionally append (true) or overwrite (false) the existing regions.
	Returns: (array:polygon, string:type[, boolean:append=true]):boolean
    */
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;
    var_land& L = V->land;

    std::string type_str = lower_case(type);

    if (type_str.find("incl") != std::string::npos)
        type_str = "inclusion";
    else if (type_str.find("excl") != std::string::npos)
        type_str = "exclusion";
    else
    {
        //invalid argument
        return false;
    }

    //convert the polygon into the required string format
    std::vector< std::string > pt, poly;
    for (size_t i = 0; i < *npts; i++)
    {
        pt.clear();
        for (size_t j = 0; j < *ndim; j++)
            pt.push_back(std::to_string(polygon_points[i*(*ndim) + j]));

        poly.push_back(join(pt, ","));
    }

    std::string spoly = "[POLY]" + join(poly, "[P]");

    if (type_str == "inclusion")
    {
        if (!is_append)
            L.inclusions.val.clear();

        L.inclusions.set_from_string(spoly.c_str());
    }
    else
    {
        if (!is_append)
            L.exclusions.val.clear();

        L.exclusions.set_from_string(spoly.c_str());
    }

    L.is_bounds_array.val = true;

    return true;

}

SPEXPORT sp_number_t* sp_heliostats_by_region(sp_data_t p_data, const char* coor_sys, int* lenret, 
                                                sp_number_t* arguments = NULL, int* len_arg = NULL, 
                                                const char* svgfname_data = NULL, sp_number_t* svg_opt_tab = NULL)
{
    /*
    Returns heliostats that fall within a region. Options are:
    >> all (no additional arguments),
    >> cylindrical (provide [rmin,rmax,azmin,azmax radians]),
    >> cartesian (provide [xmin, xmax, ymin, ymax[, zmin, zmax]]),
    >> polygon (provide [[x1,y1],[x2,y2],...]),
    >> svg (provide string with 'scale-x scale-y;offset-x offset-y;<svg path 1>;<svg path 2>;...',
    >> svgfile (provide string filename, optional table {'offset'=array, 'scale'=array}).
    (string:system, variant:region info[, string:return info - id/location])

    Returns an array of included heliostat ID's and locations. : array
    */

    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    SolarField* SF = &mc->solarfield;
    Hvector *helios = SF->getHeliostats();
    SimControl* SC = &mc->sim_control;

    //which coordinate system?
    std::string system = (std::string) coor_sys;

    //return vector -> length is unknown a priori
    std::vector<double> ret;

    if (helios->size() < 1)
    {
        std::string msg = "ERROR: sp_heliostats_by_region requires a pre-existing field to be modified.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return nullptr;
    }

    double delta = 0.0001;

    if (lower_case(system) == "all")
    {
        for (size_t i = 0; i < helios->size(); i++)
        {
            Heliostat* helio = helios->at(i);

            ret.push_back((double)helio->getId());
            ret.push_back(helio->getLocation()->x);
            ret.push_back(helio->getLocation()->y);
            ret.push_back(helio->getLocation()->z);
        }
    }
    else if (lower_case(system) == "cylindrical")
    {
        double rmin, rmax, azmin, azmax;
        if (*len_arg != 4)
            return nullptr;
     
        rmin = arguments[0] - delta;
        rmax = arguments[1] + delta;
        azmin = arguments[2] - delta;
        azmax = arguments[3] + delta;

        for (size_t i = 0; i < helios->size(); i++)
        {
            double rpos = helios->at(i)->getRadialPos();
            double apos = helios->at(i)->getAzimuthalPos();

            if (rpos > rmin)
                if (rpos < rmax)
                    if (apos > azmin)
                        if (apos < azmax)
                        {
                            ret.push_back((double)helios->at(i)->getId());
                            ret.push_back(helios->at(i)->getLocation()->x);
                            ret.push_back(helios->at(i)->getLocation()->y);
                            ret.push_back(helios->at(i)->getLocation()->z);
                        }
        }

    }
    else if (lower_case(system) == "cartesian")
    {
        if ((*len_arg != 4) && (*len_arg != 6))
            return nullptr;

        double xmin, xmax, ymin, ymax, zmin, zmax;
        xmin = arguments[0] - delta;
        xmax = arguments[1] + delta;
        ymin = arguments[2] - delta;
        ymax = arguments[3] + delta;
        if (*len_arg == 6)
        {
            zmin = arguments[4] - delta;
            zmax = arguments[5] + delta;
        }
        else
        {
            zmin = -9e9;
            zmax = 9e9;
        }


        for (size_t i = 0; i < helios->size(); i++)
        {
            sp_point *loc = helios->at(i)->getLocation();

            if (loc->x > xmin)
                if (loc->x < xmax)
                    if (loc->y > ymin)
                        if (loc->y < ymax)
                            if (loc->z > zmin)
                                if (loc->z < zmax)
                                {
                                    ret.push_back((double)helios->at(i)->getId());
                                    ret.push_back(helios->at(i)->getLocation()->x);
                                    ret.push_back(helios->at(i)->getLocation()->y);
                                    ret.push_back(helios->at(i)->getLocation()->z);
                                }
        }
    }
    else if (lower_case(system) == "polygon")
    {
        //construct a polygon from the listed points
        std::vector< sp_point > polygon;
        for (size_t i = 0; i < *len_arg; i+=2)
        {
            polygon.push_back(sp_point(arguments[i], arguments[i+1], 0.));
        }

        for (size_t i = 0; i < helios->size(); i++)
        {
            if (Toolbox::pointInPolygon(polygon, *helios->at(i)->getLocation()))
            {
                ret.push_back((double)helios->at(i)->getId());
                ret.push_back(helios->at(i)->getLocation()->x);
                ret.push_back(helios->at(i)->getLocation()->y);
                ret.push_back(helios->at(i)->getLocation()->z);
            }
        }

    }
    else if ((lower_case(system) == "svg") || (lower_case(system) == "svgfile"))
    {
        /*
        undocumented feature

        Provide string. Coordinates are space separated, points are comma separated, shapes are semicolon separated.
        First entry is x-scale y-scale, Second entry is x-offset y-offset.

        */

        std::vector< std::string > entries;
        std::vector< std::string > scale_s;
        std::vector< std::string > offset_s;


        if (lower_case(system) == "svgfile")
        {

            if (!ioutil::file_exists(svgfname_data))
            {
                std::string msg = "Invalid SVG file - not found.";
                SC->message_callback(msg.c_str(), SC->message_callback_data);
                return nullptr;
            }

            if (svg_opt_tab != NULL)
            {
                scale_s.push_back(std::to_string(svg_opt_tab[0]));
                scale_s.push_back(std::to_string(svg_opt_tab[1]));

                offset_s.push_back(std::to_string(svg_opt_tab[2]));
                offset_s.push_back(std::to_string(svg_opt_tab[3]));
            }
            else
            {
                scale_s.push_back("1.");
                scale_s.push_back("1.");

                offset_s.push_back("0.");
                offset_s.push_back("0.");
            }

            //load the svg file and parse it as an xml document
            using namespace rapidxml;
            //Read in the file to a string
            std::string file;        //contents of the file
            std::string eol;
            ioutil::read_file((std::string) svgfname_data, file, eol);

            char *fstr = new char[file.size() + 1];
            strncpy(fstr, (const char*)file.c_str(), file.size());
            fstr[file.size()] = 0;    //Null terminator

            xml_document<> doc;
            doc.parse<0>(fstr);
            xml_node<> *top_node = doc.first_node();    //<data>

            xml_node<> *node = top_node->first_node("g");
            xml_node<> *tnode = node->first_node("g");
            if (tnode)
                node = tnode;
            node = node->first_node("path");

            //assume that this is consistent with SVG files created by inkscape.. I don't know whether other SVG creators have a consistent XML structure. 
            entries.clear();
            while (node)
            {
                entries.push_back(node->first_attribute("d")->value());
                node = node->next_sibling("path");
            }
        }
        else
        {
            //get the string data and break it up into units
            if (!svgfname_data)
            {
                std::string msg = "svg data must be provided for the svg option.";
                SC->message_callback(msg.c_str(), SC->message_callback_data);
                return nullptr;
            }
            std::string data = (std::string) svgfname_data;
            entries = split(data, ";");
            scale_s = split(entries.front(), " ");
            offset_s = split(entries.at(1), " ");

            entries.erase(entries.begin(), entries.begin() + 2);
        }

        //get the scale and offset vectors
        double scale_x, scale_y, offset_x, offset_y;

        to_double(scale_s.at(0), &scale_x);
        to_double(scale_s.at(1), &scale_y);
        to_double(offset_s.at(0), &offset_x);
        to_double(offset_s.at(1), &offset_y);

        //allocate the main polygons structure
        std::vector< std::vector< sp_point > > polygons;

        for (size_t i = 0; i < entries.size(); i++)
        {
            polygons.push_back(std::vector< sp_point >());
            std::vector< sp_point > *P = &polygons.back();

            Toolbox::poly_from_svg(entries.at(i), *P, true);

            for (size_t j = 0; j < P->size(); j++)
            {
                P->at(j).x = P->at(j).x * scale_x + offset_x;
                P->at(j).y = P->at(j).y * scale_y + offset_y;
            }
        }

        //check each heliostat to see if it's in any polygon
        for (size_t i = 0; i < helios->size(); i++)
        {
            for (size_t j = 0; j < polygons.size(); j++)
            {
                std::vector< sp_point > *polygon = &polygons.at(j);
                sp_point *loc = helios->at(i)->getLocation();

                if (Toolbox::pointInPolygon(*polygon, *loc))
                {
                    ret.push_back((double)helios->at(i)->getId());
                    ret.push_back(loc->x);
                    ret.push_back(loc->y);
                    ret.push_back(loc->z);
                    //if included, don't need to check other polygons
                    break;
                }
            }
        }
    }
    else
    {
        std::string msg = "Invalid region type specified. Expecting one of [cylindrical, cartesian, polygon]";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return nullptr;
    }
    // pass back size of return vector and vector values
    *lenret = (int)ret.size();
    sp_number_t* retvec = new sp_number_t[*lenret];
    for (int i = 0; i < ret.size(); i++)
        retvec[i] = ret.at(i);

    return retvec;
}

SPEXPORT bool sp_modify_heliostats(sp_data_t p_data, sp_number_t* helio_data, int* nhel, int* ncols, const char* table_hdr)
{

    /*
    Modify attributes of a subset of heliostats in the current layout. Modifiable attributes include 
    location/aimpoint/soiling/reflectivity/enabled, 
    and are specified in the input table by variable name and array pairs. 
    The length of each variable array must match the length of the ID array. For example, 
    modify_heliostats([1,2,3], { 'location'=[[1.,10.],[2.,11.],[3.,12.]] } );, (array:heliostat IDs, table:variables)
    
    Returns: boolean
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    SolarField* SF = &mc->solarfield;
    SimControl* SC = &mc->sim_control;
    
    if (SF->getHeliostats()->size() < 1)
    {
        std::string msg = "ERROR: sp_modify_heliostats requires a pre-existing field to be modified.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    //validity check of provided headers
        //these are the supported options
    std::vector< std::string > attrs = {
        "id",
        "location-x",
        "location-y",
        "location-z",
        "aimpoint-x",
        "aimpoint-y",
        "aimpoint-z",
        "soiling",
        "reflectivity",
        "enabled"
    };
    
    //get the variable table header
    std::string hdr_str(table_hdr);
    std::vector<std::string> vars;
    vars = split(hdr_str, ",");

    if (std::find(vars.begin(), vars.end(), "id") == vars.end())
    {
        std::string msg = "ERROR: sp_modify_heliostats requires heliostat 'id' attribute.";
        SC->message_callback(msg.c_str(), SC->message_callback_data);
        return false;
    }

    //processing user input data
    std::unordered_map<std::string, std::vector<double>> datamap;
    for (size_t i = 0; i < vars.size(); i++)
    {
        if (std::find(attrs.begin(), attrs.end(), vars.at(i)) == attrs.end())
        {
            std::string msg = "Invalid attribute specified: " + vars.at(i) + "\n";
            msg += "Valid attributes names are as follows:\n";
            for (size_t j = 0; j < attrs.size(); j++)
            {
                msg += attrs.at(j) + "\n";
            }
            SC->message_callback(msg.c_str(), SC->message_callback_data);
            return false;
        }

        datamap[vars.at(i)].clear();
        for (size_t h = 0; h < (*nhel); h++)
        {
            if (vars.at(i) == "id")
                datamap[vars.at(i)].push_back((int)helio_data[h * (*ncols) + i]);
            else
                datamap[vars.at(i)].push_back((double)helio_data[h * (*ncols) + i]);       
        }
    }

    //consolidate all heliostat's into a vector by ID
    unordered_map< int, Heliostat* > *hmap = SF->getHeliostatsByID();
    Hvector helios;

    for (size_t i = 0; i < datamap.at("id").size(); i++)
    {
        try
        {
            helios.push_back(hmap->at((int)datamap.at("id").at(i)));
        }
        catch (...)
        {
            std::string msg = "Invalid id specified: " + std::to_string(datamap.at("id").at(i));
            SC->message_callback(msg.c_str(), SC->message_callback_data);
            return false;
        }
    }

    //Push back the rest of the field
    for (unordered_map<int, Heliostat*>::iterator col = hmap->begin(); col != hmap->end(); col++)
    {
        int id = col->first;
        Heliostat* helio = col->second;
        // if id not within datamap.at("id") then add helio to helios
        if (std::find(datamap.at("id").begin(), datamap.at("id").end(), id) == datamap.at("id").end())
        {
            helios.push_back(helio);
        }
    }

    //locations need to be modified through the layout shell object
    layout_shell* layout = SF->getLayoutShellObject();

    //assign location(s)
    layout->clear();

    // Creating layout object with previous set values
    for (size_t i = 0; i < helios.size(); i++)
    {
        //update the layout object
        layout->push_back(layout_obj());
        layout_obj& lobj = layout->back();

        Heliostat* helio = helios.at(i);

        lobj.aim = *helio->getAimPoint();
        lobj.cant = *helio->getCantVector();
        lobj.focal_x = helio->getFocalX();
        lobj.focal_y = helio->getFocalY();
        lobj.helio_type = helio->getMasterTemplate()->getId();
        lobj.location = *helio->getLocation();

        //update enabled/in layout statuses
        lobj.is_enabled = helio->IsEnabled();
        lobj.is_in_layout = helio->IsInLayout();
    }

    for (std::unordered_map<std::string, std::vector<double>>::iterator col = datamap.begin(); col != datamap.end(); col++)
    {
        std::string varname = col->first;
        std::vector<double>& vardata = col->second;

        if (varname == "location-x")
        {
            for (size_t j = 0; j < vardata.size(); j++)
                layout->at(j).location.x = vardata.at(j);
        }
        else if (varname == "location-y")
        {
            for (size_t j = 0; j < vardata.size(); j++)
                layout->at(j).location.y = vardata.at(j);
        }
        else if (varname == "location-z")
        {
            for (size_t j = 0; j < vardata.size(); j++)
                layout->at(j).location.z = vardata.at(j);
        }
        else if (varname == "aimpoint-x")
        {
            for (size_t j = 0; j < vardata.size(); j++)
            {
                layout->at(j).aim.x = vardata.at(j);
                layout->at(j).is_user_aim = true;
            }
        }
        else if (varname == "aimpoint-y")
        {
            for (size_t j = 0; j < vardata.size(); j++)
            {
                layout->at(j).aim.y = vardata.at(j);
                layout->at(j).is_user_aim = true;
            }
        }
        else if (varname == "aimpoint-z")
        {
            for (size_t j = 0; j < vardata.size(); j++)
            {
                layout->at(j).aim.z = vardata.at(j);
                layout->at(j).is_user_aim = true;
            }
        }
        
        else if (varname == "enabled")
        {
            for (size_t j = 0; j < vardata.size(); j++)
            {
                if ((int)vardata.at(j) == 1)
                    layout->at(j).is_enabled = true;
                else
                    layout->at(j).is_enabled = false;
            }
        }
    }

    SF->PrepareFieldLayout(*SF, 0, true);

    //not sure if this is needed...
    Hvector* updated_helios = SF->getHeliostats();

    std::vector<std::string> post_layout_cols = { "soiling", "reflectivity" };

    for (size_t s = 0; s < post_layout_cols.size(); s++)
    {
        std::string varname = post_layout_cols.at(s);

        if (datamap.find(varname) == datamap.end())
            continue;

        std::vector<double>& vardata = datamap.at(varname);

        if (varname == "soiling")
        {
            for (size_t j = 0; j < vardata.size(); j++)
                updated_helios->at(j)->getEfficiencyObject()->soiling = vardata.at(j);
        }
        else if (varname == "reflectivity")
        {
            for (size_t j = 0; j < vardata.size(); j++)
                updated_helios->at(j)->getEfficiencyObject()->reflectivity = vardata.at(j);
        }
    }

    return true;
}


SPEXPORT bool sp_calculate_optical_efficiency_table(sp_data_t p_data, int ud_n_az = NULL, int ud_n_zen = NULL)
{
    /*
    Generates a solar field optical efficiency table based an even spacing of azimuths and zeniths angles.
    
    Param: ud_n_az: User-defined number of azimuth sampling points
    Param: ud_n_zen: User-defined number of zenith sampling points

    Returns: boolean
    */

    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;
    SolarField* SF = &mc->solarfield;
    SimControl* SC = &mc->sim_control;
    sp_optical_table* opttab = &mc->opttab;

    SC->_cancel_simulation = false;

    //set the solar positions to calculate
    if (ud_n_az == NULL && ud_n_zen == NULL) {
        //set the solar positions for calculation to the default values
        double def_az[] = { 0.,  30.,  60.,  90., 120., 150., 180., 210., 240., 270., 300., 330. };
        double def_zen[] = { 0.50,   7.,  15.,  30.,  45.,  60.,  75.,  85. };

        //set the solar positions based on defaults
        opttab->azimuths.assign(def_az, def_az + 12);
        opttab->zeniths.assign(def_zen, def_zen + 8);
    }
    else {
        //set the solar positions based on even spacing and user-defined size
        opttab->azimuths.clear();
        double az_delta = 360. / double(ud_n_az);
        for (int i = 0; i < ud_n_az; i++)
            opttab->azimuths.push_back(0 + i * az_delta);

        opttab->zeniths.clear();
        double zen_delta = 85. / double(ud_n_zen - 1); // removing low angles solar positions
        for (int i = 0; i < ud_n_zen; i++)
            opttab->zeniths.push_back(0 + i * zen_delta);
    }

    double dni = SF->getVarMap()->sf.dni_des.val;
    sim_params P;
    P.dni = dni;
    P.Tamb = 25.;
    P.Patm = 1.;

    //Which type of simulation is this?
    int simtype = V->flux.flux_model.mapval();    //0=Delsol, 1=Soltrace
    Hvector* helios = SF->getHeliostats();

    sim_results* results = &mc->results;
    results->clear();
    results->resize(1);

    int neff_az = opttab->azimuths.size();
    int neff_zen = opttab->zeniths.size();
    int k = 0;
    bool allsimok = true;
    for (int j = 0; j < neff_zen; j++) {
        for (int i = 0; i < neff_az; i++) {
            //Update the solar position
            V->flux.flux_solar_az_in.set_from_string(my_to_string(opttab->azimuths.at(i)).c_str());
            V->flux.flux_solar_el_in.set_from_string(my_to_string(90. - opttab->zeniths.at(j)).c_str());
            V->flux.flux_time_type.combo_select_by_mapval(var_fluxsim::FLUX_TIME_TYPE::SUN_POSITION);

            //Set-up simulation
            if (!interop::PerformanceSimulationPrep(*SF, *helios, simtype)) return false;
            Vect sun = Ambient::calcSunVectorFromAzZen(V->flux.flux_solar_az.Val() * D2R, (90. - V->flux.flux_solar_el.Val()) * D2R);
            SF->calcHeliostatShadows(sun);
            if (SF->ErrCheck()) return false;
            
            //Which type of simulation?
            bool simok;
            switch (simtype)
            {
            case var_fluxsim::FLUX_MODEL::HERMITE_ANALYTICAL:
                simok = interop::HermiteFluxSimulationHandler(*results, *SF, *helios);
                break;
            case var_fluxsim::FLUX_MODEL::SOLTRACE:
                simok = interop::SolTraceFluxSimulation(*SC, *results, *SF, *V, *helios);
                break;
            default:
                return false;
            }

            allsimok = allsimok && simok;
            results->push_back(sim_result());

            if (SC->_cancel_simulation)
                return false;
        }
    }

    //collect all of the results and process into the efficiency table data structure
    opttab->eff_data.clear();
    k = 0;
    for (int j = 0; j < neff_zen; j++) {
        std::vector<double> row;
        for (int i = 0; i < neff_az; i++) {
            row.push_back(results->at(k).eff_total_sf.wtmean / results->at(k).eff_absorption.wtmean); // reporting field optical efficiency only
            k++;
        }
        opttab->eff_data.push_back(row);
    }

    return allsimok;
}


SPEXPORT sp_number_t* sp_get_optical_efficiency_table(sp_data_t p_data, int* nrows, int* ncols)
{
    /*
    Gets optical efficiency table in the following format:
        First row: Elevation angles (with leading zero)
        First column: Azimuth angles
        Rest of columns: Optical efficiency corresponding to elevation angle

    Param: nrows: number of rows
    Param: ncols: number of cols

    Returns: boolean
    */
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    sp_optical_table* opttab = &mc->opttab;

    // add one to account for zenith (row) and azimuth (col) 
    *nrows = (int)opttab->azimuths.size() + 1;
    *ncols = (int)opttab->zeniths.size() + 1;

    sp_number_t* efftab = new sp_number_t[(*nrows) * (*ncols)];

    //rows
    for (size_t i = 0; i < *nrows; i++)
    {
        //cols
        for (size_t j = 0; j < *ncols; j++)
        {
            if (i == 0 && j == 0) {
                efftab[i * (*ncols) + j] = 0; // leading zero
            }
            else if (i == 0) {
                efftab[i * (*ncols) + j] = (90. - opttab->zeniths.at(j - 1)); // Elevation angle
            }
            else if (i != 0 && j == 0) {
                efftab[i * (*ncols) + j] = opttab->azimuths.at(i - 1);
            }
            else {
                efftab[i * (*ncols) + j] = opttab->eff_data.at(j - 1).at(i - 1);
            }
        }
    }
    return efftab;
}


SPEXPORT bool sp_save_optical_efficiency_table(sp_data_t p_data, const char* sp_fname, const char* table_name)
{
    /*
    Saves optical efficiency table as a CSV file in the following format:
        First row: Elevation angles (with leading zero)
        First column: Azimuth angles
        Rest of columns: Optical efficiency corresponding to elevation angle        

    Param: sp_fname: location to save efficiency table file
    Param: table_name: specific table name for Modelica table output file. 
        If equal to "none", then table format follows sp_get_optical_efficiency_table. 
        Otherwise, table format is consistent with Modelica with extra header lines and elevation angle to be in ascending order.

    Returns: boolean
    */

    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    sp_optical_table* opttab = &mc->opttab;

    std::string fname(sp_fname);
    if (!ioutil::dir_exists(ioutil::path_only(fname).c_str()))
    {
        return false;
    }

    try
    {
        std::ofstream of(fname);
        std::string sep = ", ";

        int neff_zen = opttab->zeniths.size();
        int neff_az = opttab->azimuths.size();

        if (of.is_open())      //check whether the file is accessible
        {
            // Specific Modelica table file formatting
            std::string tablename(table_name);
            if (tablename != "none") {
                of << "#1\n";
                of << "double " << tablename << "(" << to_string(neff_az + 1) << "," << to_string(neff_zen + 1) << ")\n";
            }

            // First row is elevation angles with leading 0
            of << "0";
            for (int i = 0; i < neff_zen; i++) {
                int ind = i;
                if (tablename != "none") {
                    ind = (neff_zen - (int)1) - i; // ascending order
                }
                std::string el_val = to_string(90. - opttab->zeniths.at(ind));  
                of << sep << el_val;
            }
            of << "\n";

            // Rest of rows lead with azimuth then efficiency at azimuth and corresponding zenith
            for (int i = 0; i < neff_az; i++) {
                std::string az_val = to_string(opttab->azimuths.at(i));
                of << az_val;
                for (int j = 0; j < neff_zen; j++) {
                    int ind = j;
                    if (tablename != "none") {
                        ind = (neff_zen - (int)1) - j; // ascending order
                    }
                    std::string eff_val = to_string(opttab->eff_data.at(ind).at(i)); // ascending order
                    of << sep << eff_val;
                }
                of << "\n";
            }
        }
        of.close();

        return true;
    }
    catch (...)
    {
        return false;
    }
}

SPEXPORT bool sp_save_from_script(sp_data_t p_data, const char* sp_fname)
{
    /*
	Save the current case as a SolarPILOT .spt file. Returns true if successful.
	Returns: (string:path):boolean
    */

    CopilotObject *mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;

    std::string fname(sp_fname);
    if (!ioutil::dir_exists(ioutil::path_only(fname).c_str()))
    {
        return false;
    }

    try
    {
        parametric p;
        optimization o;

        ioutil::saveXMLInputFile(fname, *V, p, o,"C++_API");
        return true;
    }
    catch (...)
    {

    }
    return false;
}

SPEXPORT bool sp_load_from_script(sp_data_t p_data, const char* sp_fname)
{
    /*
    Loads a SolarPILOT .spt file. Returns true if successful.
    Returns: boolean
    */

    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;

    if (!ioutil::file_exists(sp_fname))
    {
        return false;
    }

    try
    {
        parametric p;
        optimization o;

        ioutil::parseXMLInputFile(sp_fname, *V, p, o);
        return true;
    }
    catch (...)
    {

    }
    return false;
}

SPEXPORT bool sp_dump_varmap(sp_data_t p_data, const char* sp_fname)
{
    /*
	Dump the variable structure to a text csv file. Returns true if successful.
	Returns: (string:path):boolean
    */
    CopilotObject* mc = static_cast<CopilotObject*>(p_data);
    var_map* V = &mc->variables;
    SimControl* SC = &mc->sim_control;

    //valid path?
    std::string fname(sp_fname);
    if (!ioutil::dir_exists(ioutil::path_only(fname).c_str()))
    {
        return false;
    }


    try
    {
        //create a list of all variable keys
        std::vector< std::string > names;
        names.reserve(V->_varptrs.size());

        for (unordered_map< std::string, spbase* >::iterator it = V->_varptrs.begin(); it != V->_varptrs.end(); it++)
            names.push_back(it->first);

        sort(names.begin(), names.end());   //output the variables alphabetically

        std::ofstream of( fname );
        std::string sep = " ,";

        if (of.is_open())      //check whether the file is accessible
        {
            for (size_t i = 0; i < names.size(); i++)
            {
                spbase *var = V->_varptrs.at(names.at(i));

                std::string val = var->as_string();

                if (val.size() > 30)      //tuncate very long values
                {
                    val.erase(val.begin() + 20, val.end());
                    val += "... (truncated)";
                }

                //replace all commas
                std::string::size_type n = 0;
                while ((n = val.find(",", n)) != std::string::npos)
                {
                    val.replace(n, 1, "'");
                    n += 1;
                }

                std::string units = var->units;
                //replace all commas
                n = 0;
                while ((n = units.find(",", n)) != std::string::npos)
                {
                    units.replace(n, 1, "'");
                    n += 1;
                }

                of << names.at(i) << sep << val << sep << units << sep << var->short_desc << "\n";
            }

            of.close();
            return true;  //success
        }

    }
    catch (std::exception &e)
    {
        SC->message_callback(e.what(), SC->message_callback_data);
        //std::runtime_error(e.what());
    }
    return false;
}

SPEXPORT void _sp_free_var(sp_number_t* m)
{
    try
    {
        if(m!=0)
            delete[] m;

    }
    catch (...)
    {
        return;
    }
}


int ST_APICallback(st_uint_t ntracedtotal, st_uint_t ntraced, st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages, void* data)
{
    /* 
    
    return 1 if ok, return 0 to kill soltrace simulation
    */
    CopilotObject* api= static_cast<CopilotObject*>(data);

    if (api->use_api_callback)
    {
        std::string messages = "";
        if (api->message_log.size() != 0)
        {
            messages = join(api->message_log, "\n");
        }
        api->message_log.clear();

        (*api->external_callback)((sp_number_t)((double)ntraced / (double)(std::max((int)ntotrace, 1))), messages.c_str());
    }
    return 1;
};

int MessageHandler(const char* message, void* data)
{
    /*
    type        |   name        | Description
    -------------------------------------------
    const char* | message       | Messages produced by the simulation
    void*       | data          | Pointer to data (to be type-cast) accompanying the message
    */
    
    
    CopilotObject* api = static_cast<CopilotObject*>(data);

    if (api->use_api_callback)
    {
        (*api->external_callback)((sp_number_t)0., message);
    }
    return 1;
}

int ProgressHandler(double progress, const char* message, void* data)
{
    /*
    type        |   name        | Description
    -------------------------------------------
    double      | progress      | Progress of the simulation [0...1]
    const char* | message       | Messages produced by the simulation
    void*       | data          | Pointer to data (to be type-cast) accompanying the message
    */


    CopilotObject* api = static_cast<CopilotObject*>(data);

    if (api->use_api_callback)
    {
        (*api->external_callback)((sp_number_t)progress, message);
    }
    return 1;
}