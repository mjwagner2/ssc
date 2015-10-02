#include <vector>
#include <string>
#include "csp_solver_core.h"
#include <unordered_map>
#include "lp_lib.h" 
//#include "glpk\src\glpk.h"

using namespace std;

#ifndef _CSP_DISPATCH
#define _CSP_DISPATCH

class csp_dispatch_opt
{
    int  m_nstep_opt;              //number of time steps in the optimized array
    bool m_is_weather_setup;  //bool indicating whether the weather has been copied
    
    void clear_output_arrays();

public:
    bool m_last_opt_successful;   //last optimization run was successful?
    int m_current_read_step;        //current step to read from optimization results
    vector<double> price_signal;    //IN [- or $/MWh] Price factor indicating market value of generated energy
    C_csp_weatherreader m_weather;       //Local copy of weather reader object

    struct s_solver_params
    {
        bool is_abort_flag;         //optimization flagged for abort
        int iter_count;             //branch and bound iteration count
        string log_message;
        double obj_relaxed;

        //user settings
        int max_bb_iter;            //Maximum allowable iterations for B&B algorithm
        double mip_gap;             //convergence tolerance - gap between relaxed MIP solution and current best solution
        double solution_timeout;    //[s] Max solve time for each solution
        int presolve_type;
        int bb_type;  
        int scaling_type;

        s_solver_params()
        {
            bb_type = -1;
            presolve_type = -1;
            scaling_type = -1;
        };

        void reset()
        {
            is_abort_flag = false;
            iter_count = 0;
            log_message.clear();
            obj_relaxed = 0.;
        };

    } solver_params;

    struct s_params
    {
        bool is_rec_operating0;     //receiver is operating at the initial time step
        bool is_pb_operating0;      //Power block is operating at the initial time step
        bool is_pb_standby0;        //Power block is in standby at the initial time step
        double dt;                  //Time step (hr)
        double e_tes_init;          //[kWht] current stored energy capacity
        double e_tes_min;           //[kWht] minimum allowable energy capacity in TES
        double e_tes_max;           //[kWht] maximum allowable energy capacity in TES
        double e_pb_startup_cold;   //[kWht] energy requirement to start up the power block from cold state
        double e_pb_startup_hot;    //[kWht] energy requirement to start up the power block from standby
        double e_rec_startup;       //[kWht] energy requirement to start up the reciever
        double dt_pb_startup_cold;  //[hr] time requiremeent to start up the power block from cold state
        double dt_pb_startup_hot;   //[hr] time requiremeent to start up the power block from hot state
        double dt_rec_startup;      //[hr] time requirement to start up the reciever
        double tes_degrade_rate;    //IN [1/hr] Fractional energy loss from tes per hour
        double q_pb_standby;        //[kWt] power requirement to maintain the power block in standby mode
        double q_pb_des;               //[kWe] design cycle thermal power input
        double q_pb_max;            //[kWt] Maximum allowable thermal energy rate to the cycle
        double q_pb_min;            //[kWt] Minimum allowable thermal energy rate to the cycle
        double q_rec_min;           //[kWt] Minimum allowable power delivery by the receiver when operating
        double w_rec_pump;          //[kWe/kWt] Pumping parasitic power per thermal energy produced
        double sf_effadj;           //[-] 0..1 Solar field efficiency adjustment
        double info_time;           //[s] time of the year at sim start. informational only.

        C_csp_solver_sim_info *siminfo;     //Pointer to existing simulation info object
        C_csp_collector_receiver *col_rec;   //Pointer to collector/receiver object
        C_csp_messages *messages;   //Pointer to message structure

        struct s_efftable
        {
        private:
            struct s_effmember
            {
                double x;
                double eta;

                s_effmember(){};
                s_effmember(double _x, double _eta)
                {
                    x = _x;
                    eta = _eta;
                };
            };
            vector<s_effmember> table;

        public:

            void clear()
            {
                table.clear();
            }

            void add_point(double x, double eta)
            {
                table.push_back( s_effmember(x, eta) );
            };

            bool get_point(int index, double &x, double &eta)
            {
                if( index > table.size()-1 || index < 0 ) return false;

                x = table.at(index).x;
                eta = table.at(index).eta;
            }

            double get_point_eff(int index)
            {
                return table.at(index).eta;
            }

            double get_point_x(int index)
            {
                return table.at(index).x;
            }

            size_t get_size()
            {
                return table.size();
            }

            double interpolate(double x)
            {

                double eff = table.front().eta;

                int ind = 0;
                int ni = (int)table.size();
                while( true )
                {
                    if( ind ==  ni-1 )
                    {
                        eff = table.back().eta;
                        break;
                    }

                    if( x < table.at(ind).x )
                    {
                        if(ind == 0)
                        {
                            eff = table.front().eta;
                        }
                        else
                        {
                            eff = table.at(ind-1).eta + (table.at(ind).eta - table.at(ind-1).eta)*(x - table.at(ind-1).x)/(table.at(ind).x - table.at(ind-1).x);
                        }
                        break;
                    }

                    ind ++;
                }

                return eff;
            }

        } eff_table_load, eff_table_Tdb;        //Efficiency of the power cycle
        
    } params;

    struct s_outputs
    {
        double objective;
        double objective_relaxed;
        vector<bool> rec_operation;   //receiver startup ok?
        vector<bool> pb_operation;    //power block startup ok?
        vector<bool> pb_standby;    //power block standby ok?
        vector<double> q_pb_target;       //optimized energy generation (less startup loss)
        vector<double> q_pb_standby;      //standby energy allowed
        vector<double> q_sfavail_expected;       //Expected available solar field energy
        vector<double> q_sf_expected;           //Expected solar field energy generation
        vector<double> eta_pb_expected;     //Expected power cycle conversion efficiency (normalized)
        vector<double> eta_sf_expected;     //Expected solar field thermal efficiency (normalized)
        vector<double> tes_charge_expected;     //Expected thermal energy storage charge state
        vector<double> q_pb_startup;    //thermal power going to startup
        vector<double> q_rec_startup;   //thermal power going to startup

        int solve_iter;             //Number of iterations required to solve
        int solve_state;
    } outputs;
    
    //----- public member functions ----

    csp_dispatch_opt();

    //check parameters and inputs to make sure everything has been set up correctly
    bool check_setup(int nstep);

    //copy the weather data over
    bool copy_weather_data(C_csp_weatherreader &weather_source);

    //Predict performance out nstep values. 
    bool predict_performance(int step_start, int nstep);    

    //declare dispatch function in csp_dispatch.cpp
    static int optimize_demo();

    static int optimize_demo2();

    bool optimize();

    //Get optimized variable states for this timestep
    //s_outputs *get_step_vars(int step);

    
    

};

template <class T>
class Array_base
{
protected:
    T *dat;
    int *starts;
    int dim_size;
    int mem_size;
    void index_exception();

public:
    string var_name;
    Array_base();
    Array_base(int n, string varname="");
    T *data_array();
    string getvarname();
    int getmemsize();
    void deallocate();

    virtual void allocate(int n, bool zeros=true) = 0;
    virtual T &at(int row, int col) = 0;
};

// ----------------------------------------
template <class T>
class Array_T : public Array_base<T>
{
public:
    void allocate(int n, bool zeros=true);
    T &at(int row, int col=0);
    T &at(int row);
};

template <class T>
class Array_2T : public Array_base<T>
{
public:
    void allocate(int n, bool zeros=true);
    T &at(int row, int col);
};

template <class T>
class Array_2T_Tri : public Array_base<T>
{
public:
    void allocate(int n, bool zeros=true);
    T &at(int row, int col);
};

// ----------------------------------------



class optimization_vars
{
    int current_mem_pos;
    int alloc_mem_size; 

    REAL *data;
public:
    struct opt_var
    {
        string name;
        int var_type;
        int var_dim;
        int var_dim_size;
        int var_dim_size2;
        int ind_start;
        int ind_end;
        REAL upper_bound;
        REAL lower_bound;
    };
private: 
    vector<opt_var> var_objects;

    unordered_map<string, opt_var*> var_by_name;

public:
    struct VAR_TYPE { enum A {REAL_T, INT_T, BINARY_T}; };
    struct VAR_DIM { enum A {DIM_T, DIM_NT, DIM_T2, DIM_2T_TRI}; };

    optimization_vars();
    //~optimization_vars();

    void add_var(char *vname, int var_type /* VAR_TYPE enum */, int var_dim /* VAR_DIM enum */, int var_dim_size, REAL lowbo=-DEF_INFINITE, REAL upbo=DEF_INFINITE);
    void add_var(char *vname, int var_type /* VAR_TYPE enum */, int var_dim /* VAR_DIM enum */, int var_dim_size, int var_dim_size2, REAL lowbo=-DEF_INFINITE, REAL upbo=DEF_INFINITE);

    bool construct();

    int get_num_varobjs();
    int get_total_var_count();

    REAL &operator()(char *varname, int ind);    //Access for 1D var
    REAL &operator()(char *varname, int ind1, int ind2);     //Access for 2D var
    REAL &operator()(int varindex, int ind);    
    REAL &operator()(int varindex, int ind1, int ind2);

    int column(char *varname, int ind);
    int column(char *varname, int ind1, int ind2);
    int column(int varindex, int ind);
    int column(int varindex, int ind1, int ind2);

    REAL *get_variable_array(); 

    opt_var *get_var(char *varname);
    opt_var *get_var(int varindex);
};




#endif