/**
BSD-3-Clause
Copyright 2019 Alliance for Sustainable Energy, LLC
Redistribution and use in source and binary forms, with or without modification, are permitted provided 
that the following conditions are met :
1.	Redistributions of source code must retain the above copyright notice, this list of conditions 
and the following disclaimer.
2.	Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.	Neither the name of the copyright holder nor the names of its contributors may be used to endorse 
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER, CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES 
DEPARTMENT OF ENERGY, NOR ANY OF THEIR EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "core.h"

// for adjustment factors
#include "common.h"

// solarpilot header files
#include "AutoPilot_API.h"
#include "SolarField.h"
#include "IOUtil.h"
#include "csp_common.h"

// Can probably delete these headers later...
#include "csp_solver_util.h"
#include "csp_solver_core.h"
#include "csp_solver_pt_sf_perf_interp.h"
//#include "csp_solver_mspt_receiver_222.h"
#include "csp_solver_mspt_receiver.h"
#include "csp_solver_mspt_collector_receiver.h"
#include "csp_solver_tower_collector_receiver.h"
#include "csp_solver_pc_Rankine_indirect_224.h"
#include "csp_solver_pc_indirect_Gen3.h"
#include "csp_solver_pc_sco2.h"
#include "csp_solver_two_tank_tes.h"
#include "csp_solver_tou_block_schedules.h"

#include "csp_system_costs.h"

static var_info _cm_vtab_tcsmolten_salt[] = {
    // VARTYPE       DATATYPE    NAME                                  LABEL                                                                                                                                      UNITS           META                                 GROUP                                       REQUIRED_IF                                                         CONSTRAINTS      UI_HINTS
    { SSC_INPUT,     SSC_STRING, "solar_resource_file",                "Local weather file path",                                                                                                                 "",             "",                                  "Location and Resource",                    "?",                                                                "LOCAL_FILE",    ""},
    { SSC_INPUT,     SSC_TABLE,  "solar_resource_data",                "Weather resource data in memory",                                                                                                         "",             "",                                  "Location and Resource",                    "?",                                                                "",              ""},

    { SSC_INPUT,     SSC_NUMBER, "ppa_multiplier_model",               "PPA multiplier model",                                                                                                                    "0/1",          "0=diurnal,1=timestep",              "Time of Delivery Factors",                 "?=0",                                                              "INTEGER,MIN=0", ""},
    { SSC_INPUT,     SSC_ARRAY,  "dispatch_factors_ts",                "Dispatch payment factor array",                                                                                                           "",             "",                                  "Time of Delivery Factors",                 "ppa_multiplier_model=1",                                           "",              ""},

    //{ SSC_INPUT,     SSC_NUMBER, "field_model_type",                   "0=design field and tower/receiver geometry, 1=design field, 2=user specified field, 3=user performance maps vs solar position",           "",             "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "gross_net_conversion_factor",        "Estimated gross to net conversion factor",                                                                                                "",             "",                                  "System Design",                            "*",                                                                "",              ""},

    { SSC_INPUT,     SSC_NUMBER, "helio_width",                        "Heliostat width",                                                                                                                         "m",            "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "helio_height",                       "Heliostat height",                                                                                                                        "m",            "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    //{ SSC_INPUT,     SSC_NUMBER, "helio_optical_error_mrad",           "Heliostat optical error",                                                                                                                 "mrad",         "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "helio_active_fraction",              "Heliostat active fraction",                                                                                                               "",             "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dens_mirror",                        "Ratio of heliostat reflective area to profile",                                                                                           "",             "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "helio_reflectance",                  "Heliostat reflectance",                                                                                                                   "",             "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "rec_absorptance",                    "Receiver absorptance",                                                                                                                    "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "rec_hl_perm2",                       "Receiver design heatloss",                                                                                                                "kW/m2",        "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "land_max",                           "Land max boundary",                                                                                                                       "-ORm",         "",                                  "Heliostat Field",                          "?=7.5",                                                            "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "land_min",                           "Land min boundary",                                                                                                                       "-ORm",         "",                                  "Heliostat Field",                          "?=0.75",                                                           "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "land_bound_table",                   "Land boundary table",                                                                                                                     "m",            "",                                  "Heliostat Field",                          "?",                                                                "",              ""},
    { SSC_INPUT,     SSC_ARRAY,  "land_bound_list",                    "Land boundary table listing",                                                                                                             "",             "",                                  "Heliostat Field",                          "?",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "p_start",                            "Heliostat startup energy",                                                                                                                "kWe-hr",       "",                                  "Heliostat Field",                          "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "p_track",                            "Heliostat tracking energy",                                                                                                               "kWe",          "",                                  "Heliostat Field",                          "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "hel_stow_deploy",                    "Stow/deploy elevation angle",                                                                                                             "deg",          "",                                  "Heliostat Field",                          "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "v_wind_max",                         "Heliostat max wind velocity",                                                                                                             "m/s",          "",                                  "Heliostat Field",                          "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "interp_nug",                         "Interpolation nugget",                                                                                                                    "-",            "",                                  "Heliostat Field",                          "?=0",                                                              "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "interp_beta",                        "Interpolation beta coef.",                                                                                                                "-",            "",                                  "Heliostat Field",                          "?=1.99",                                                           "",              "" },
    { SSC_INPUT,     SSC_MATRIX, "eta_map",                            "Field efficiency array",                                                                                                                  "",             "",                                  "Heliostat Field",                          "?",                                                                "",              "" },

    { SSC_INPUT,     SSC_NUMBER, "water_usage_per_wash",               "Water usage per wash",                                                                                                                    "L/m2_aper",    "",                                  "Heliostat Field",                          "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "washing_frequency",                  "Mirror washing frequency",                                                                                                                "none",         "",                                  "Heliostat Field",                          "*",                                                                "",              "" },

    { SSC_INPUT,     SSC_NUMBER, "tower_fixed_cost",                   "Tower fixed cost",                                                                                                                        "$",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "tower_exp",                          "Tower cost scaling exponent",                                                                                                             "",             "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "foundation_fixed_cost",              "Tower foundation fixed cost",                                                                                                             "$",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "foundation_cost_scaling_quadratic",  "Tower foundation cost scaling quadratic",                                                                                                 "$/m2",         "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "foundation_cost_scaling_linear",     "Tower foundation cost scaling linear",                                                                                                    "$/m",          "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "particle_lift_cost",                 "Particle lift fixed cost",                                                                                                                "$",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "riser_and_downcomer_cost",           "Combined riser and downcomer fixed cost",                                                                                                 "$",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "rec_ref_cost",                       "Receiver reference cost",                                                                                                                 "$",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "rec_ref_area",                       "Receiver reference area for cost scale",                                                                                                  "",             "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "rec_cost_exp",                       "Receiver cost scaling exponent",                                                                                                          "",             "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "site_spec_cost",                     "Site improvement cost",                                                                                                                   "$/m2",         "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "heliostat_spec_cost",                "Heliostat field cost",                                                                                                                    "$/m2",         "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "plant_spec_cost",                    "Power cycle specific cost",                                                                                                               "$/kWe",        "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "bop_spec_cost",                      "BOS specific cost",                                                                                                                       "$/kWe",        "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "tes_spec_cost",                      "Thermal energy storage cost",                                                                                                             "$/kWht",       "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "land_spec_cost",                     "Total land area cost",                                                                                                                    "$/acre",       "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "contingency_rate",                   "Contingency for cost overrun",                                                                                                            "%",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "sales_tax_rate",                     "Sales tax rate",                                                                                                                          "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "sales_tax_frac",                     "Percent of cost to which sales tax applies",                                                                                              "%",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "cost_sf_fixed",                      "Solar field fixed cost",                                                                                                                  "$",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "fossil_spec_cost",                   "Fossil system specific cost",                                                                                                             "$/kWe",        "",                                  "System Costs",                             "*",                                                                "",              "" },

        //other costs needed for optimization update
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.epc.per_acre",           "EPC cost per acre",                                                                                                                       "$/acre",       "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.epc.percent.smaller",    "EPC cost percent of direct, smaller plants",                                                                                              "%",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.epc.percent.larger",     "EPC cost percent of direct, larger plants",                                                                                               "%",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.epc.per_watt",           "EPC cost per watt",                                                                                                                       "$/W",          "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.epc.fixed.smaller",      "EPC fixed, smaller plants",                                                                                                               "$",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.epc.fixed.larger",       "EPC fixed, larger plants",                                                                                                                "$",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.plm.percent",            "PLM cost percent of direct",                                                                                                              "%",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.plm.per_watt",           "PLM cost per watt",                                                                                                                       "$/W",          "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.cost.plm.fixed",              "PLM fixed",                                                                                                                               "$",            "",                                  "System Costs",                             "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.sf.fixed_land_area",          "Fixed land area",                                                                                                                         "acre",         "",                                  "Heliostat Field",                          "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.sf.land_overhead_factor",     "Land overhead factor",                                                                                                                    "",             "",                                  "Heliostat Field",                          "*",                                                                "",              "" },

        // System Design
    //    { SSC_INPUT,     SSC_NUMBER, "T_htf_cold_des",                     "Cold HTF inlet temperature at design conditions",                                                                                         "C",            "",                                  "System Design",                            "*",                                                                "",              ""},
    //    { SSC_INPUT,     SSC_NUMBER, "T_htf_hot_des",                      "Hot HTF outlet temperature at design conditions",                                                                                         "C",            "",                                  "System Design",                            "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "P_ref",                              "Reference output electric power at design condition",                                                                                     "MW",           "",                                  "System Design",                            "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "design_eff",                         "Power cycle efficiency at design",                                                                                                        "none",         "",                                  "System Design",                            "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "tshours",                            "Equivalent full-load thermal storage hours",                                                                                              "hr",           "",                                  "System Design",                            "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "solarm",                             "Solar multiple",                                                                                                                          "-",            "",                                  "System Design",                            "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "dni_des",                            "Design-point DNI",                                                                                                                        "W/m2",         "",                                  "System Design",                            "*",                                                                "",              ""},

        // Receiver (type 222) parameters
    { SSC_INPUT,     SSC_MATRIX, "rec_efficiency_lookup",              "Receiver efficiency lookup table: eta(load, Tin)",                                                                                        "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_MATRIX, "rec_pressure_lookup",                "Receiver pressure drop lookup table: dP(load, Pin)",                                                                                      "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },

    { SSC_INPUT,     SSC_MATRIX, "field_fl_props",                     "User defined field fluid property data",                                                                                                  "-",            "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "rec_htf",                            "Receiver HTF, 17=Salt (60% NaNO3, 40% KNO3) 10=Salt (46.5% LiF 11.5% NaF 42% KF) 50=Lookup tables",                                       "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "f_rec_min",                          "Minimum receiver mass flow rate turn down fraction",                                                                                      "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "rec_su_delay",                       "Fixed startup delay time for the receiver",                                                                                               "hr",           "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "rec_qf_delay",                       "Energy-based receiver startup delay (fraction of rated thermal power)",                                                                   "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.rec.max_oper_frac",           "Maximum receiver mass flow rate fraction",                                                                                                "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "eta_pump",                           "Receiver HTF pump efficiency",                                                                                                            "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "piping_loss",                        "Thermal loss per meter of piping",                                                                                                        "Wt/m",         "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "piping_length_mult",                 "Piping length multiplier",                                                                                                                "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "piping_length_const",                "Piping constant length",                                                                                                                  "m",            "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "piping_loss_coeff",                  "Wetted loss coefficient for riser or downcomer",                                                                                          "W/m2/K",       "",                                  "Tower and Receiver",                       "?=5.0",                                                            "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "piping_riser_diam",                  "Piping riser inner diameter",                                                                                                             "m",            "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "piping_downcomer_diam",              "Piping downcomer inner diameter",                                                                                                         "m",            "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "T_rec_cold_des",                     "Individual receiver design cold temperature",                                                                                             "C",            "",                                  "Tower and Receiver",                       "",                                                                 "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "T_rec_hot_des",                      "Individual receiver design hot temperature",                                                                                              "C",            "",                                  "Tower and Receiver",                       "",                                                                 "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "L_recHX",                            "Length of charge HX",                                                                                                                     "m",            "",                                  "Tower and Receiver",                       "",                                                                 "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "n_cells_recHX",                      "Number of cells in charge HX",                                                                                                            "-",            "",                                  "Tower and Receiver",                       "",                                                                 "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "is_rec_recirc_available",            "1: Receiver has option to use recirculator, 0: receiver cannot produce heat unless PC is ON",                                             "",             "",                                  "Tower and Receiver",                       "?=0",                                                              "",              ""},

    
    // TES parameters - general
    { SSC_INPUT,     SSC_NUMBER, "store_htf",                          "Storage HTF, 17=Salt (60% NaNO3, 40% KNO3) 10=Salt (46.5% LiF 11.5% NaF 42% KF) 50=Lookup tables",                                        "",             "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "store_fl_props",                     "User defined storage fluid property data",                                                                                                "-",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "T_tes_hot_des",                      "TES tank design hot temperature",                                                                                                         "C",            "",                                  "Thermal Storage",                          "",                                                                 "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "T_tes_warm_des",                     "TES (virtual) tank design warm temperature, between the high-temp and low-temp HXs",                                                      "C",            "",                                  "Thermal Storage",                          "",                                                                 "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "T_tes_cold_des",                     "TES tank design cold temperature",                                                                                                        "C",            "",                                  "Thermal Storage",                          "",                                                                 "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "csp.pt.tes.init_hot_htf_percent",    "Initial fraction of available volume that is hot",                                                                                        "%",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "h_tank",                             "Total height of tank (height of HTF when tank is full)",                                                                                  "m",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "cold_tank_max_heat",                 "Rated heater capacity for cold tank heating",                                                                                             "MW",           "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dt_charging",                        "Charging HX approach temp",                                                                                                               "C",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dt_ht_discharging",                  "High temperature discharging HX approach temp",                                                                                           "C",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dt_lt_discharging",                  "Low temperature discharging HX approach temp",                                                                                            "C",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "u_tank",                             "Loss coefficient from the tank",                                                                                                          "W/m2-K",       "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "tank_pairs",                         "Number of equivalent tank pairs",                                                                                                         "",             "",                                  "Thermal Storage",                          "*",                                                                "INTEGER",       ""},
    { SSC_INPUT,     SSC_NUMBER, "cold_tank_Thtr",                     "Minimum allowable cold tank HTF temperature",                                                                                             "C",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    // TES parameters - 2 tank
    { SSC_INPUT,     SSC_NUMBER, "h_tank_min",                         "Minimum allowable HTF height in storage tank",                                                                                            "m",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "hot_tank_Thtr",                      "Minimum allowable hot tank HTF temperature",                                                                                              "C",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "hot_tank_max_heat",                  "Rated heater capacity for hot tank heating",                                                                                              "MW",           "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "L_LTHX",                             "Length of low-temp HX",                                                                                                                   "m",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "L_HTHX",                             "Length of high-temp HX",                                                                                                                  "m",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "n_cells_LTHX",                       "Number of cells in low-temp HX",                                                                                                          "-",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "n_cells_HTHX",                       "Number of cells in high-temp HX",                                                                                                         "-",            "",                                  "Thermal Storage",                          "*",                                                                "",              ""},

    // Power Cycle Inputs
    { SSC_INPUT,     SSC_NUMBER, "T_pc_cold_des",                      "Power cycle design cold temperature",                                                                                                     "C",            "",                                  "Power Cycle",                              "",                                                                 "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "T_pc_hot_des",                       "Power cycle design hot temperature",                                                                                                      "C",            "",                                  "Power Cycle",                              "",                                                                 "",              ""},
    //{ SSC_INPUT,     SSC_NUMBER, "pc_config",                          "PC configuration 0=Steam Rankine (224), 1=user defined, 2=sCO2 Recompression (424)",                                                      "",             "",                                  "Power Cycle",                              "?=0",                                                              "INTEGER",       ""},
    { SSC_INPUT,     SSC_NUMBER, "is_udpc_co2",                        "Is the user-defined power cycle an integrated sco2 cycle, where the HTF is also the cycle working fluid",                                 "",             "",                                  "Power Cycle",                              "?=0",                                                              "INTEGER",       ""},
    { SSC_INPUT,     SSC_NUMBER, "startup_time",                       "Time needed for power block startup",                                                                                                     "hr",           "",                                  "Power Cycle",                              "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "startup_frac",                       "Fraction of design thermal power needed for startup",                                                                                     "none",         "",                                  "Power Cycle",                              "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "cycle_max_frac",                     "Maximum turbine over design operation fraction",                                                                                          "",             "",                                  "Power Cycle",                              "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "cycle_cutoff_frac",                  "Minimum turbine operation fraction before shutdown",                                                                                      "",             "",                                  "Power Cycle",                              "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "q_sby_frac",                         "Fraction of thermal power required for standby",                                                                                          "",             "",                                  "Power Cycle",                              "*",                                                                "",              ""},

    // User Defined cycle
    { SSC_INPUT,     SSC_NUMBER, "ud_T_amb_des",                       "Ambient temperature at user-defined power cycle design point",                                                                            "C",            "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_f_W_dot_cool_des",                "Percent of user-defined power cycle design gross output consumed by cooling",                                                             "%",            "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_m_dot_water_cool_des",            "Mass flow rate of water required at user-defined power cycle design point",                                                               "kg/s",         "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_T_htf_low",                       "Low level HTF inlet temperature for T_amb parametric",                                                                                    "C",            "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_T_htf_high",                      "High level HTF inlet temperature for T_amb parametric",                                                                                   "C",            "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_T_amb_low",                       "Low level ambient temperature for HTF mass flow rate parametric",                                                                         "C",            "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_T_amb_high",                      "High level ambient temperature for HTF mass flow rate parametric",                                                                        "C",            "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_m_dot_htf_low",                   "Low level normalized HTF mass flow rate for T_HTF parametric",                                                                            "",             "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "ud_m_dot_htf_high",                  "High level normalized HTF mass flow rate for T_HTF parametric",                                                                           "",             "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "ud_T_htf_ind_od",                    "Off design table of user-defined power cycle performance formed from parametric on T_htf_hot [C]",                                        "",             "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "ud_T_amb_ind_od",                    "Off design table of user-defined power cycle performance formed from parametric on T_amb [C]",                                            "",             "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "ud_m_dot_htf_ind_od",                "Off design table of user-defined power cycle performance formed from parametric on m_dot_htf [ND]",                                       "",             "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "ud_ind_od",                          "Off design user-defined power cycle performance as function of T_htf, m_dot_htf [ND], and T_amb",                                         "",             "",                                  "User Defined Power Cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "ud_ind_od_off_sun",                  "Off design for OFF SUN user-defined power cycle performance as function of T_htf, m_dot_htf [ND], and T_amb",                                         "",             "",                      "User Defined Power Cycle",                 "*",                                                      "",              "" },

    // Direct CO2 User Defined
    { SSC_INPUT,     SSC_NUMBER, "P_phx_in_co2_des",                   "CO2 PHX inlet pressure",                                                                                                                  "kPa",          "",                                  "User Defined Power cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "P_turb_in_co2_des",                  "CO2 turbine inlet pressure",                                                                                                              "kPa",          "",                                  "User Defined Power cycle",                 "*",                                                      "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "P_turb_in_co2_off_sun_des",          "CO2 turbine inlet pressure at off sun",                                                                                                   "kPa",          "",                                  "User Defined Power cycle",                 "*",                                                      "",              ""},

    // System Control
    { SSC_INPUT,     SSC_NUMBER, "time_start",                         "Simulation start time",                                                                                                                   "s",            "",                                  "System Control",                           "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "time_stop",                          "Simulation stop time",                                                                                                                    "s",            "",                                  "System Control",                           "?=31536000",                                                       "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "time_steps_per_hour",                "Number of simulation time steps per hour",                                                                                                "",             "",                                  "System Control",                           "?=-1",                                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "vacuum_arrays",                      "Allocate arrays for only the required number of steps",                                                                                   "",             "",                                  "System Control",                           "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "pb_fixed_par",                       "Fixed parasitic load - runs at all times",                                                                                                "MWe/MWcap",    "",                                  "System Control",                           "*",                                                                "",              ""},
    //{ SSC_INPUT,     SSC_NUMBER, "aux_par",                            "Aux heater, boiler parasitic",                                                                                                            "MWe/MWcap",    "",                                  "System Control",                           "*",                                                                "",              ""},
    //{ SSC_INPUT,     SSC_NUMBER, "aux_par_f",                          "Aux heater, boiler parasitic - multiplying fraction",                                                                                     "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    //{ SSC_INPUT,     SSC_NUMBER, "aux_par_0",                          "Aux heater, boiler parasitic - constant coefficient",                                                                                     "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    //{ SSC_INPUT,     SSC_NUMBER, "aux_par_1",                          "Aux heater, boiler parasitic - linear coefficient",                                                                                       "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    //{ SSC_INPUT,     SSC_NUMBER, "aux_par_2",                          "Aux heater, boiler parasitic - quadratic coefficient",                                                                                    "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "bop_par",                            "Balance of plant parasitic power fraction",                                                                                               "MWe/MWcap",    "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "bop_par_f",                          "Balance of plant parasitic power fraction - mult frac",                                                                                   "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "bop_par_0",                          "Balance of plant parasitic power fraction - const coeff",                                                                                 "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "bop_par_1",                          "Balance of plant parasitic power fraction - linear coeff",                                                                                "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "bop_par_2",                          "Balance of plant parasitic power fraction - quadratic coeff",                                                                             "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_ARRAY,  "f_turb_tou_periods",                 "Dispatch logic for turbine load fraction",                                                                                                "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "weekday_schedule",                   "12x24 CSP operation Time-of-Use Weekday schedule",                                                                                        "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "weekend_schedule",                   "12x24 CSP operation Time-of-Use Weekend schedule",                                                                                        "",             "",                                  "System Control",                           "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "is_tod_pc_target_also_pc_max",       "Is the TOD target cycle heat input also the max cycle heat input?",                                                                       "",             "",                                  "System Control",                           "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "is_dispatch",                        "Allow dispatch optimization?",                                                                                                            "",             "",                                  "System Control",                           "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_horizon",                       "Time horizon for dispatch optimization",                                                                                                  "hour",         "",                                  "System Control",                           "is_dispatch=1",                                                    "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_frequency",                     "Frequency for dispatch optimization calculations",                                                                                        "hour",         "",                                  "System Control",                           "is_dispatch=1",                                                    "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_steps_per_hour",                "Time steps per hour for dispatch optimization calculations",                                                                              "",             "",                                  "System Control",                           "?=1",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_max_iter",                      "Max number of dispatch optimization iterations",                                                                                          "",             "",                                  "System Control",                           "is_dispatch=1",                                                    "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_timeout",                       "Max dispatch optimization solve duration",                                                                                                "s",            "",                                  "System Control",                           "is_dispatch=1",                                                    "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_mip_gap",                       "Dispatch optimization solution tolerance",                                                                                                "",             "",                                  "System Control",                           "is_dispatch=1",                                                    "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_spec_bb",                       "Dispatch optimization B&B heuristic",                                                                                                     "",             "",                                  "System Control",                           "?=-1",                                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_reporting",                     "Dispatch optimization reporting level",                                                                                                   "",             "",                                  "System Control",                           "?=-1",                                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_spec_presolve",                 "Dispatch optimization presolve heuristic",                                                                                                "",             "",                                  "System Control",                           "?=-1",                                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_spec_scaling",                  "Dispatch optimization scaling heuristic",                                                                                                 "",             "",                                  "System Control",                           "?=-1",                                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_time_weighting",                "Dispatch optimization future time discounting factor",                                                                                    "",             "",                                  "System Control",                           "?=0.99",                                                           "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "is_write_ampl_dat",                  "Write AMPL data files for dispatch run",                                                                                                  "",             "",                                  "System Control",                           "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_STRING, "ampl_data_dir",                      "AMPL data file directory",                                                                                                                "",             "",                                  "System Control",                           "?=''",                                                             "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "is_ampl_engine",                     "Run dispatch optimization with external AMPL engine",                                                                                     "",             "",                                  "System Control",                           "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_STRING, "ampl_exec_call",                     "System command to run AMPL code",                                                                                                         "",             "",                                  "System Control",                           "?='ampl sdk_solution.run'",                                        "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_rsu_cost",                      "Receiver startup cost",                                                                                                                   "$",            "",                                  "System Control",                           "is_dispatch=1",                                                    "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_csu_cost",                      "Cycle startup cost",                                                                                                                      "$",            "",                                  "System Control",                           "is_dispatch=1",                                                    "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "disp_pen_delta_w",                   "Dispatch cycle production change penalty",                                                                                                "$/kWe-change", "",                                  "System Control",                           "is_dispatch=1",                                                    "",              "" },
    { SSC_INPUT,     SSC_NUMBER, "disp_inventory_incentive",           "Dispatch storage terminal inventory incentive multiplier",                                                                                "",             "",                                  "System Control",                           "?=0.0",                                                            "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "q_rec_standby",                      "Receiver standby energy consumption",                                                                                                     "kWt",          "",                                  "System Control",                           "?=9e99",                                                           "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "q_rec_heattrace",                    "Receiver heat trace energy consumption during startup",                                                                                   "kWe-hr",       "",                                  "System Control",                           "?=0.0",                                                            "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "is_wlim_series",                     "Use time-series net electricity generation limits",                                                                                       "",             "",                                  "System Control",                           "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_ARRAY,  "wlim_series",                        "Time series net electicity generation limits",                                                                                            "kWe",          "",                                  "System Control",                           "is_wlim_series=1",                                                 "",              ""},


    // Financial inputs
    { SSC_INPUT,     SSC_MATRIX, "dispatch_sched_weekday",             "PPA pricing weekday schedule, 12x24",                                                                                                     "",             "",                                  "Time of Delivery Factors",                 "?=[[1]]",                                                          "",              ""},
    { SSC_INPUT,     SSC_MATRIX, "dispatch_sched_weekend",             "PPA pricing weekend schedule, 12x24",                                                                                                     "",             "",                                  "Time of Delivery Factors",                 "?=[[1]]",                                                          "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor1",                   "Dispatch payment factor 1",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "?=1",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor2",                   "Dispatch payment factor 2",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "?=1",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor3",                   "Dispatch payment factor 3",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "?=1",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor4",                   "Dispatch payment factor 4",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "?=1",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor5",                   "Dispatch payment factor 5",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "?=1",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor6",                   "Dispatch payment factor 6",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "?=1",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor7",                   "Dispatch payment factor 7",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "?=1",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor8",                   "Dispatch payment factor 8",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "?=1",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "dispatch_factor9",                   "Dispatch payment factor 9",                                                                                                               "",             "",                                  "Time of Delivery Factors",                 "?=1",                                                              "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "is_dispatch_series",                 "Use time-series dispatch factors",                                                                                                        "",             "",                                  "System Control",                           "?=0",                                                              "",              ""},
    { SSC_INPUT,     SSC_ARRAY,  "dispatch_series",                    "Time series dispatch factors",                                                                                                            "",             "",                                  "System Control",                           "",                                                                 "",              ""},


    // optimized outputs updated depending on run type 
    { SSC_INOUT,     SSC_NUMBER, "rec_height",                         "Receiver height",                                                                                                                         "m",            "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INOUT,     SSC_NUMBER, "D_rec",                              "The overall outer diameter of the receiver",                                                                                              "m",            "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    { SSC_INOUT,     SSC_NUMBER, "h_tower",                            "Tower height",                                                                                                                            "m",            "",                                  "Tower and Receiver",                       "*",                                                                "",              ""},
    //{ SSC_INOUT,     SSC_NUMBER, "N_hel",                              "Number of heliostats",                                                                                                                    "",             "",                                  "Heliostat Field",                          "",                                                                 "",              ""},
    //{ SSC_INOUT,     SSC_MATRIX, "helio_positions",                    "Heliostat position table",                                                                                                                "",             "",                                  "Heliostat Field",                          "*",                                                                "",              "COL_LABEL=XY_POSITION"},
    { SSC_INOUT,     SSC_NUMBER, "land_area_base",                     "Base land area occupied by heliostats",                                                                                                   "acre",         "",                                  "Heliostat Field",                          "*",                                                                "",              ""},
    
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.total_land_area",        "Total land area",                                                                                                                         "acre",         "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.site_improvements",      "Site improvement cost",                                                                                                                   "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.heliostats",             "Heliostat cost",                                                                                                                          "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.tower",                  "Tower cost",                                                                                                                              "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.receiver",               "Receiver cost",                                                                                                                           "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.storage",                "TES cost",                                                                                                                                "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.power_block",            "Power cycle cost",                                                                                                                        "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.rad_field",              "Radiative field cost"                                                                                                                     "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.rad_fluid",              "Radiative fluid cost"                                                                                                                     "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.rad_storage",            "Cold storage cost"                                                                                                                        "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.bop",                    "BOP cost",                                                                                                                                "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.fossil",                 "Fossil backup cost",                                                                                                                      "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "ui_direct_subtotal",                 "Direct capital precontingency cost",                                                                                                      "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.contingency",            "Contingency cost",                                                                                                                        "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "total_direct_cost",                  "Total direct cost",                                                                                                                       "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.epc.total",              "EPC and owner cost",                                                                                                                      "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.plm.total",              "Total land cost",                                                                                                                         "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.sales_tax.total",        "Sales tax cost",                                                                                                                          "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "total_indirect_cost",                "Total indirect cost",                                                                                                                     "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "total_installed_cost",               "Total installed cost",                                                                                                                    "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "csp.pt.cost.installed_per_capacity", "Estimated installed cost per cap",                                                                                                        "$",            "",                                  "System Costs",                             "*",                                                                "",              ""},

        // Construction financing inputs/outputs (SSC variable table from cmod_cb_construction_financing)
    { SSC_INPUT,     SSC_NUMBER, "const_per_interest_rate1",           "Interest rate, loan 1",                                                                                                                   "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_interest_rate2",           "Interest rate, loan 2",                                                                                                                   "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_interest_rate3",           "Interest rate, loan 3",                                                                                                                   "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_interest_rate4",           "Interest rate, loan 4",                                                                                                                   "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_interest_rate5",           "Interest rate, loan 5",                                                                                                                   "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_months1",                  "Months prior to operation, loan 1",                                                                                                       "",             "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_months2",                  "Months prior to operation, loan 2",                                                                                                       "",             "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_months3",                  "Months prior to operation, loan 3",                                                                                                       "",             "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_months4",                  "Months prior to operation, loan 4",                                                                                                       "",             "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_months5",                  "Months prior to operation, loan 5",                                                                                                       "",             "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_percent1",                 "Percent of total installed cost, loan 1",                                                                                                 "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_percent2",                 "Percent of total installed cost, loan 2",                                                                                                 "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_percent3",                 "Percent of total installed cost, loan 3",                                                                                                 "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_percent4",                 "Percent of total installed cost, loan 4",                                                                                                 "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_percent5",                 "Percent of total installed cost, loan 5",                                                                                                 "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_upfront_rate1",            "Upfront fee on principal, loan 1",                                                                                                        "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_upfront_rate2",            "Upfront fee on principal, loan 2",                                                                                                        "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_upfront_rate3",            "Upfront fee on principal, loan 3",                                                                                                        "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_upfront_rate4",            "Upfront fee on principal, loan 4",                                                                                                        "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_INPUT,     SSC_NUMBER, "const_per_upfront_rate5",            "Upfront fee on principal, loan 5",                                                                                                        "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal1",               "Principal, loan 1",                                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal2",               "Principal, loan 2",                                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal3",               "Principal, loan 3",                                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal4",               "Principal, loan 4",                                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal5",               "Principal, loan 5",                                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest1",                "Interest cost, loan 1",                                                                                                                   "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest2",                "Interest cost, loan 2",                                                                                                                   "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest3",                "Interest cost, loan 3",                                                                                                                   "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest4",                "Interest cost, loan 4",                                                                                                                   "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest5",                "Interest cost, loan 5",                                                                                                                   "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_total1",                   "Total financing cost, loan 1",                                                                                                            "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_total2",                   "Total financing cost, loan 2",                                                                                                            "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_total3",                   "Total financing cost, loan 3",                                                                                                            "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_total4",                   "Total financing cost, loan 4",                                                                                                            "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_total5",                   "Total financing cost, loan 5",                                                                                                            "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_percent_total",            "Total percent of installed costs, all loans",                                                                                             "%",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_principal_total",          "Total principal, all loans",                                                                                                              "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "const_per_interest_total",           "Total interest costs, all loans",                                                                                                         "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "construction_financing_cost",        "Total construction financing cost",                                                                                                       "$",            "",                                  "Financial Parameters",                     "*",                                                                "",              ""},



    // ****************************************************************************************************************************************
    // Outputs here:
    // ****************************************************************************************************************************************
        // Simulation outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "time_hr",                            "Time at end of timestep",                                                                                                                 "hr",           "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "solzen",                             "Resource solar zenith",                                                                                                                   "deg",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "solaz",                              "Resource solar azimuth",                                                                                                                  "deg",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "beam",                               "Resource beam normal irradiance",                                                                                                         "W/m2",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "tdry",                               "Resource dry Bulb temperature",                                                                                                           "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "twet",                               "Resource wet Bulb temperature",                                                                                                           "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "rh",                                 "Resource relative humidity",                                                                                                              "%",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "wspd",                               "Resource wind velocity",                                                                                                                  "m/s",          "",                                  "",                                         "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_ARRAY,  "A_sf1",                              "Solar field area assigned receiver 1",                                                                                                    "m^2",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "A_sf2",                              "Solar field area assigned receiver 2",                                                                                                    "m^2",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "A_sf3",                              "Solar field area assigned receiver 3",                                                                                                    "m^2",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "q_sf_inc",                           "Field incident thermal power",                                                                                                            "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "eta_field_tot",                      "Field optical efficiency total",                                                                                                          "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "eta_field1",                         "Field optical efficiency 1",                                                                                                              "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "eta_field2",                         "Field optical efficiency 2",                                                                                                              "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "eta_field3",                         "Field optical efficiency 3",                                                                                                              "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "defocus",                            "Field optical focus fraction",                                                                                                            "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "sf_adjust_out",                      "Field availability adjustment factor",                                                                                                    "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_rec_inc",                      "Tower incident thermal power",                                                                                                            "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "eta_therm",                          "Tower thermal efficiency",                                                                                                                "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "Q_thermal",                          "Tower thermal power to HTF less piping loss",                                                                                             "MWt",          "",                                  "",                                         "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_rec",                          "Tower mass flow rate",                                                                                                                    "kg/s",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_startup",                          "Tower startup thermal energy consumed",                                                                                                   "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_in",                           "Tower HTF inlet temperature",                                                                                                             "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_out",                          "Tower HTF outlet temperature",                                                                                                            "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_piping_losses",                    "Tower header/tower piping losses",                                                                                                        "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_thermal_loss",                     "Tower convection and emission losses",                                                                                                    "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_out_end",                      "Tower HTF outlet temperature at end of timestep",                                                                                         "C",            "",                                "CR",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_out_max",                      "Tower maximum HTF outlet temperature during timestep",                                                                                    "C",            "",                                "CR",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_panel_out_max",                    "Tower panel maximum HTF outlet temperature during timestep",                                                                              "C",            "",                                "CR",                                         "*",                                                                "",              ""},    
    { SSC_OUTPUT,    SSC_ARRAY,  "T_wall_rec_inlet",                   "Tower inlet panel wall temperature at end of timestep",                                                                                   "C",            "",                                "CR",                                         "*",                                                                "",              ""},    
    { SSC_OUTPUT,    SSC_ARRAY,  "T_wall_rec_outlet",                  "Tower outlet panel wall temperature at end of timestep",                                                                                  "C",            "",                                "CR",                                         "*",                                                                "",              ""},    
    { SSC_OUTPUT,    SSC_ARRAY,  "T_wall_riser",                       "Tower riser wall temperature at end of timestep",                                                                                         "C",            "",                                "CR",                                         "*",                                                                "",              ""},    
    { SSC_OUTPUT,    SSC_ARRAY,  "T_wall_downcomer",                   "Tower downcomer wall temperature at end of timestep",                                                                                     "C",            "",                                "CR",                                         "*",                                                                "",              ""},    

    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_rec_inc1",                     "Receiver 1 incident thermal power",                                                                                                       "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_rec_inc2",                     "Receiver 2 incident thermal power",                                                                                                       "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_rec_inc3",                     "Receiver 3 incident thermal power",                                                                                                       "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_rec3",                         "Receiver 3 mass flow rate",                                                                                                               "kg/hr",        "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_rec1",                         "Receiver 1 mass flow rate",                                                                                                               "kg/hr",        "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_rec2",                         "Receiver 2 mass flow rate",                                                                                                               "kg/hr",        "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_in1",                          "Receiver 1 HTF inlet temperature",                                                                                                        "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_in2",                          "Receiver 2 HTF inlet temperature",                                                                                                        "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_in3",                          "Receiver 3 HTF inlet temperature",                                                                                                        "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_out1",                         "Receiver 1 HTF outlet temperature",                                                                                                       "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_out2",                         "Receiver 2 HTF outlet temperature",                                                                                                       "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rec_out3",                         "Receiver 3 HTF outlet temperature",                                                                                                       "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_HX1",                          "HX 1 transferred power",                                                                                                                "MWt",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_HX2",                          "HX 2 transferred power",                                                                                                                "MWt",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_HX3",                          "HX 3 transferred power",                                                                                                                "MWt",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_HX_tes1",                      "HX 1 TES mass flow rate",                                                                                                             "kg/hr",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_HX_tes2",                      "HX 2 TES mass flow rate",                                                                                                             "kg/hr",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_HX_tes3",                      "HX 3 TES mass flow rate",                                                                                                             "kg/hr",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_HX_tes_out1",                      "HX 1 TES outlet temperature",                                                                                                             "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_HX_tes_out2",                      "HX 2 TES outlet temperature",                                                                                                             "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_HX_tes_out3",                      "HX 3 TES outlet temperature",                                                                                                             "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "eta_rec_therm1",                     "Receiver 1 thermal efficiency",                                                                                                           "-",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "eta_rec_therm2",                     "Receiver 2 thermal efficiency",                                                                                                           "-",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "eta_rec_therm3",                     "Receiver 3 thermal efficiency",                                                                                                           "-",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "dp_rec1",                            "Receiver 1 pressure drop",                                                                                                                "kPa",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "dp_rec2",                            "Receiver 2 pressure drop",                                                                                                                "kPa",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "dp_rec3",                            "Receiver 3 pressure drop",                                                                                                                "kPa",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "dp_riser",                           "Riser pressure drop",                                                                                                                     "kPa",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "dp_downcomer",                       "Downcomer pressure drop",                                                                                                                 "kPa",          "",                                  "",                                         "*",                                                                "",              "" },

    { SSC_OUTPUT,    SSC_ARRAY,  "W_dot_recirc",                       "Receiver CO2 recirculator power",                                                                                                         "MWe",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "Q_dot_to_particles",                 "Total CHX heat transfer to particles",                                                                                                    "MWt",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "Q_dot_rec_therm_1",                  "Rec 1 CO2 heat absorbed",                                                                                                                 "MWt",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "Q_dot_rec_therm_2",                  "Rec 2 CO2 heat absorbed",                                                                                                                 "MWt",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "Q_dot_rec_therm_3",                  "Rec 3 CO2 heat absorbed",                                                                                                                 "MWt",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "dp_CO2_HX_1",                        "HX 1 CO2 pressure drop",                                                                                                                  "kPa",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "dp_CO2_HX_2",                        "HX 2 CO2 pressure drop",                                                                                                                  "kPa",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "dp_CO2_HX_3",                        "HX 3 CO2 pressure drop",                                                                                                                  "kPa",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_downcomer",                    "Downcomer thermal losses",                                                                                                                "MWt",          "",                                  "",                                         "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_riser",                        "Riser thermal losses",                                                                                                                    "MWt",          "",                                  "",                                         "*",                                                                "",              "" },


        // Power cycle outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "eta",                                "PC efficiency, gross",                                                                                                                    "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_pb",                               "PC input energy",                                                                                                                         "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_pc",                           "PC HTF mass flow rate",                                                                                                                   "kg/s",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_pc_startup",                       "PC startup thermal energy",                                                                                                               "MWht",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_pc_startup",                   "PC startup thermal power",                                                                                                                "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_cycle",                            "PC electrical power output, gross",                                                                                                       "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_pc_in",                            "PC HTF inlet temperature",                                                                                                                "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_pc_out",                           "PC HTF outlet temperature",                                                                                                               "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_water_pc",                     "PC water consumption, makeup + cooling",                                                                                                  "kg/s",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_cond_out",                         "PC condenser water outlet temperature",                                                                                                   "C",            "",                                  "PC",                                       "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_cold",                             "Cold storage cold temperature",                                                                                                           "C",            "",                                  "PC",                                       "?",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_cold",                             "Cold storage cold tank mass",                                                                                                             "kg",           "",                                  "PC",                                       "?",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_warm",                             "Cold storage warm tank mass",                                                                                                             "kg",           "",                                  "PC",                                       "?",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_warm",                             "Cold storage warm tank temperature",                                                                                                      "C",            "",                                  "PC",                                       "?",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_rad_out",                          "Radiator outlet temperature",                                                                                                             "C",            "",                                  "PC",                                       "?",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "A_radfield",                         "Radiator field surface area",                                                                                                             "m^2",          "",                                  "PC",                                       "?",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_cond",                             "PC condensing presssure",                                                                                                                 "Pa",           "",                                  "PC",                                       "?",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "radcool_control",                    "Radiative cooling status code",                                                                                                           "-",            "",                                  "PC",                                       "?",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "P_phx_in",                           "Primary HX/receiver inlet pressure",                                                                                                      "kPa",          "",                                  "PC",                                       "?",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_ARRAY,  "P_phx_out",                          "Primary HX return/cycle inlet pressure",                                                                                                  "kPa",          "",                                  "PC",                                       "?",                                                                "",              "" },

        // Thermal energy storage outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "tank_losses",                        "TES thermal losses",                                                                                                                      "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_heater",                           "TES freeze protection power",                                                                                                             "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_tes_hot",                          "TES hot temperature",                                                                                                                     "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "T_tes_cold",                         "TES cold temperature",                                                                                                                    "C",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dc_tes",                           "TES discharge thermal power",                                                                                                             "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_ch_tes",                           "TES charge thermal power",                                                                                                                "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "e_ch_tes",                           "TES charge state",                                                                                                                        "MWht",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_tes_dc",                       "TES discharge mass flow rate",                                                                                                            "kg/s",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_tes_ch",                       "TES charge mass flow rate",                                                                                                               "kg/s",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "V_hot_tank_initial",                 "TES initial hot tank fill fraction",                                                                                                      "-",            "",                                  "",                                         "*",                                                                "",              ""},

        // Parasitics outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "pparasi",                            "Parasitic power heliostat drives",                                                                                                        "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_tower_pump",                       "Parasitic power receiver/tower HTF pump",                                                                                                 "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "htf_pump_power",                     "Parasitic power TES and cycle HTF pump",                                                                                                  "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_cooling_tower_tot",                "Parasitic power condenser operation",                                                                                                     "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_fixed",                            "Parasitic power fixed load",                                                                                                              "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_plant_balance_tot",                "Parasitic power generation-dependent load",                                                                                               "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "P_rec_heattrace",                    "Receiver heat trace parasitic load",                                                                                                      "MWe",          "",                                  "System",                                   "*",                                                                "",              ""},

        // System outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "P_out_net",                          "Total electric power to grid",                                                                                                            "MWe",          "",                                  "",                                         "*",                                                                "",              ""},

        // Controller outputs
    { SSC_OUTPUT,    SSC_ARRAY,  "tou_value",                          "CSP operating time-of-use value",                                                                                                         "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "pricing_mult",                       "PPA price multiplier",                                                                                                                    "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "n_op_modes",                         "Operating modes in reporting timestep",                                                                                                   "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "op_mode_1",                          "1st operating mode",                                                                                                                      "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "op_mode_2",                          "2nd operating mode, if applicable",                                                                                                       "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "op_mode_3",                          "3rd operating mode, if applicable",                                                                                                       "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "m_dot_balance",                      "Relative mass flow balance error",                                                                                                        "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_balance",                          "Relative energy balance error",                                                                                                           "",             "",                                  "",                                         "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_ARRAY,  "disp_solve_state",                   "Dispatch solver state",                                                                                                                   "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_solve_iter",                    "Dispatch iterations count",                                                                                                               "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_objective",                     "Dispatch objective function value",                                                                                                       "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_obj_relax",                     "Dispatch objective function - relaxed max",                                                                                               "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_qsf_expected",                  "Dispatch expected solar field available energy",                                                                                          "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_qsfprod_expected",              "Dispatch expected solar field generation",                                                                                                "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_qsfsu_expected",                "Dispatch expected solar field startup enegy",                                                                                             "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_tes_expected",                  "Dispatch expected TES charge level",                                                                                                      "MWht",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_pceff_expected",                "Dispatch expected power cycle efficiency adj.",                                                                                           "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_thermeff_expected",             "Dispatch expected SF thermal efficiency adj.",                                                                                            "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_qpbsu_expected",                "Dispatch expected power cycle startup energy",                                                                                            "MWht",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_wpb_expected",                  "Dispatch expected power generation",                                                                                                      "MWe",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_rev_expected",                  "Dispatch expected revenue factor",                                                                                                        "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_presolve_nconstr",              "Dispatch number of constraints in problem",                                                                                               "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_presolve_nvar",                 "Dispatch number of variables in problem",                                                                                                 "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "disp_solve_time",                    "Dispatch solver time",                                                                                                                    "sec",          "",                                  "",                                         "*",                                                                "",              ""},


        // These outputs correspond to the first csp-solver timestep in the reporting timestep.
        //     Subsequent csp-solver timesteps within the same reporting timestep are not tracked
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_pc_sb",                        "Thermal power for PC standby",                                                                                                            "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_pc_min",                       "Thermal power for PC min operation",                                                                                                      "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_pc_max",                       "Max thermal power to PC",                                                                                                                 "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_pc_target",                    "Target thermal power to PC",                                                                                                              "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "is_rec_su_allowed",                  "Is receiver startup allowed",                                                                                                             "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "is_pc_su_allowed",                   "Is power cycle startup allowed",                                                                                                          "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "is_pc_sb_allowed",                   "Is power cycle standby allowed",                                                                                                          "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_est_cr_su",                    "Estimated receiver startup thermal power",                                                                                                "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_est_cr_on",                    "Estimated receiver thermal power TO HTF",                                                                                                 "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_est_tes_dc",                   "Estimated max TES discharge thermal power",                                                                                               "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "q_dot_est_tes_ch",                   "Estimated max TES charge thermal power",                                                                                                  "MWt",          "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "operating_modes_a",                  "First 3 operating modes tried",                                                                                                           "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "operating_modes_b",                  "Next 3 operating modes tried",                                                                                                            "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_ARRAY,  "operating_modes_c",                  "Final 3 operating modes tried",                                                                                                           "",             "",                                  "",                                         "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_ARRAY,  "gen",                                "Total electric power to grid with available derate",                                                                                      "kWe",          "",                                  "",                                         "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_MATRIX, "ud_T_htf_ind_od_out",                "T_htf_hot cycle off design",                                                                                                              "",             "",                                  "",                                         "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]", "",              "COL_LABEL=UDPC_T_HTF_HOT,ROW_LABEL=NO_ROW_LABEL" },
    { SSC_OUTPUT,    SSC_MATRIX, "ud_T_amb_ind_od_out",                "T_amb cycle off design",                                                                                                                  "",             "",                                  "",                                         "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]", "",              "COL_LABEL=UDPC_T_AMB,ROW_LABEL=NO_ROW_LABEL" },
    { SSC_OUTPUT,    SSC_MATRIX, "ud_m_dot_htf_ind_od_out",            "M_dot_htf cycle off design",                                                                                                              "",             "",                                  "",                                         "?=[[0,1,2,3,4,5,6,7,8,9,10,11,12][0,1,2,3,4,5,6,7,8,9,10,11,12]]", "",              "COL_LABEL=UDPC_M_DOT_HTF,ROW_LABEL=NO_ROW_LABEL" },
    { SSC_OUTPUT,    SSC_MATRIX, "sco2_preprocess_table_out",          "sCO2 cycle preprocessed data in UDPC format",                                                                                             "",             "",                                  "",                                         "?=[[0]]",                                                          "",              "COL_LABEL=UDPC_SCO2_PREPROC,ROW_LABEL=NO_ROW_LABEL"},

    // Annual single-value outputs
    { SSC_OUTPUT,    SSC_NUMBER, "annual_energy",                      "Annual total electric power to grid",                                                                                                     "kWhe",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "annual_W_cycle_gross",               "Electrical source - power cycle gross output",                                                                                            "kWhe",         "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "annual_W_cooling_tower",             "Total of condenser operation parasitics",                                                                                                 "kWhe",         "",                                  "PC",                                       "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "annual_q_rec_inc",                   "Annual receiver incident thermal power after reflective losses",                                                                          "MWt-hr",       "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_NUMBER, "annual_q_rec_loss",                  "Annual receiver convective and radiative losses",                                                                                         "MWt-hr",       "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_NUMBER, "annual_eta_rec_th",                  "Annual receiver thermal efficiency ignoring rec reflective loss",                                                                         "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,    SSC_NUMBER, "annual_eta_rec_th_incl_refl",        "Annual receiver thermal efficiency including reflective loss",                                                                            "",             "",                                  "Tower and Receiver",                       "*",                                                                "",              "" },

    { SSC_OUTPUT,    SSC_NUMBER, "conversion_factor",                  "Gross to net conversion factor",                                                                                                          "%",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "capacity_factor",                    "Capacity factor",                                                                                                                         "%",            "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "kwh_per_kw",                         "First year kWh/kW",                                                                                                                       "kWh/kW",       "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "annual_total_water_use",             "Total annual water usage, cycle + mirror washing",                                                                                        "m3",           "",                                  "",                                         "*",                                                                "",              ""},

    { SSC_OUTPUT,    SSC_NUMBER, "disp_objective_ann",                 "Annual sum of dispatch objective function value",                                                                                         "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "disp_iter_ann",                      "Annual sum of dispatch solver iterations",                                                                                                "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "disp_presolve_nconstr_ann",          "Annual sum of dispatch problem constraint count",                                                                                         "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "disp_presolve_nvar_ann",             "Annual sum of dispatch problem variable count",                                                                                           "",             "",                                  "",                                         "*",                                                                "",              ""},
    { SSC_OUTPUT,    SSC_NUMBER, "disp_solve_time_ann",                "Annual sum of dispatch solver time",                                                                                                      "",             "",                                  "",                                         "*",                                                                "",              ""},


    var_info_invalid };

class cm_tcsmolten_salt : public compute_module
{
public:

    cm_tcsmolten_salt()
    {
        add_var_info(_cm_vtab_tcsmolten_salt);
        add_var_info(vtab_adjustment_factors);
        add_var_info(vtab_sf_adjustment_factors);
    } 

    bool relay_message(string &msg, double percent)
    {
        log(msg);
        return update(msg, (float)percent);
    }

    void exec() override
    {
        //       FILE* fp = fopen("input_log.txt", "w");

        //       std::vector<var_info*> tables = { _cm_vtab_tcsmolten_salt, vtab_adjustment_factors, vtab_sf_adjustment_factors };

        //       for (size_t t = 0; t < 3; t++)
        //       {
        //           var_info* table = tables[t];

        //           for (size_t i = 0; true; i++)
        //           {
        //               if (table[i].data_type == SSC_INVALID)
        //                   break;

        //               if (! (table[i].var_type == SSC_INPUT || table[i].var_type == SSC_INOUT))
        //                   continue;
        //               if (!is_assigned(table[i].name))
        //                   continue;

        //               fprintf(fp, "var(\"%s\", ", table[i].name);
        //           
        //               switch (table[i].data_type)
        //               {
        //               case SSC_NUMBER:
        //               {
        //                   double val = as_double(table[i].name);
        //                   if(val > 9e12)
        //                       fprintf(fp, "1e38");
        //                   else
        //                       fprintf(fp, "%f", val);
        //                   break;
        //               }
        //               case SSC_STRING:
        //                   fprintf(fp, "%s", as_string(table[i].name));
        //                   break;
        //               case SSC_ARRAY:
        //               {
        //                   size_t nval;
        //                   ssc_number_t* pval = as_array(table[i].name, &nval);
        //                   fprintf(fp, "[");
        //                   for (int j = 0; j < nval; j++)
        //                   {
        //                       if (pval[j] > 9e12)
        //                           fprintf(fp, "1e38%s", j < nval - 1 ? "," : "");
        //                       else
        //                           fprintf(fp, "%f%s", pval[j], j < nval - 1 ? "," : "");
        //                   }
        //                   fprintf(fp, "]");
        //                   break;
        //               }
        //               case SSC_MATRIX:
        //               {
        //                   size_t nrow, ncol;
        //                   ssc_number_t* pval = as_matrix(table[i].name, &nrow, &ncol);
        //                   fprintf(fp, "[");
        //                   for (int j = 0; j < nrow; j++)
        //                   {
        //                       fprintf(fp, "[");
        //                       for (int k = 0; k < ncol; k++)
        //                           fprintf(fp, "%f%s", pval[j*ncol + k], k < ncol - 1 ? "," : "");
        //                       fprintf(fp, "]%s", j < nrow-1 ? "," : "");
        //                   }
        //                   fprintf(fp, "]");
        //                   break;
        //               }
        //               default:
        //                   break;
        //               }

        //               fprintf(fp, ");\n");
        //           }
        //       }

        //       fclose(fp);

               //FILE* fp = fopen("cmod_to_lk_script.lk", "w");

               //write_cmod_to_lk_script(fp, m_vartab);



               // Weather reader
        C_csp_weatherreader weather_reader;
        if (is_assigned("solar_resource_file")) {
            weather_reader.m_weather_data_provider = make_shared<weatherfile>(as_string("solar_resource_file"));
            if (weather_reader.m_weather_data_provider->has_message()) log(weather_reader.m_weather_data_provider->message(), SSC_WARNING);
        }
        if (is_assigned("solar_resource_data")) {
            weather_reader.m_weather_data_provider = make_shared<weatherdata>(lookup("solar_resource_data"));
            if (weather_reader.m_weather_data_provider->has_message()) log(weather_reader.m_weather_data_provider->message(), SSC_WARNING);
        }

        weather_reader.m_trackmode = 0;
        weather_reader.m_tilt = 0.0;
        weather_reader.m_azimuth = 0.0;
        // Initialize to get weather file info
        weather_reader.init();
        if (weather_reader.has_error()) throw exec_error("tcsmolten_salt", weather_reader.get_error());

        // Get info from the weather reader initialization
        double site_elevation = weather_reader.ms_solved_params.m_elev;     //[m]


        int tes_type = 1;
        const int N_rec = 3;        // number of receivers and corresponding fields

        assign("rec_aspect", as_number("rec_height") / as_number("D_rec"));

        double q_design = as_number("P_ref") / as_number("design_eff") * as_number("solarm");
        assign("q_design", q_design);

        // Calculate system capacity instead of pass in
        double system_capacity = as_double("P_ref") * as_double("gross_net_conversion_factor") * 1.E3;       //[kWe]

        //forced assignments here

        // 'sf_model_type'
        // 3 = user performance maps vs solar position
        assign("field_model_type", 3);
        assign("is_rec_model_trans", 0);
        assign("is_rec_startup_trans", 0);

        assign("is_optimize", 0);
        bool is_optimize = false;

        //assignments for field_model_type == 3 / The following optional inputs must be set here:
        assign("calc_fluxmaps", 0);

        if (tes_type != 1)
        {
            throw exec_error("MSPT CSP Solver", "Thermocline thermal energy storage is not yet supported by the new CSP Solver and Dispatch Optimization models.\n");
        }


        // Set steps per hour
        C_csp_solver::S_sim_setup sim_setup;
        sim_setup.m_sim_time_start = as_double("time_start");       //[s] time at beginning of first time step
        sim_setup.m_sim_time_end = as_double("time_stop");          //[s] time at end of last time step

        int steps_per_hour = (int)as_double("time_steps_per_hour");     //[-]

        //if the number of steps per hour is not provided (=-1), then assign it based on the weather file step
        if (steps_per_hour < 0)
        {
            double sph_d = 3600. / weather_reader.m_weather_data_provider->step_sec();
            steps_per_hour = (int)(sph_d + 1.e-5);
            if ((double)steps_per_hour != sph_d)
                throw spexception("The time step duration must be evenly divisible within an hour.");
        }

        size_t n_steps_fixed = (size_t)steps_per_hour * 8760;   //[-]
        if (as_boolean("vacuum_arrays"))
        {
            n_steps_fixed = steps_per_hour * (size_t)((sim_setup.m_sim_time_end - sim_setup.m_sim_time_start) / 3600.);
        }
        //int n_steps_fixed = (int)( (sim_setup.m_sim_time_end - sim_setup.m_sim_time_start) * steps_per_hour / 3600. ) ; 
        sim_setup.m_report_step = 3600.0 / (double)steps_per_hour;  //[s]


        // *********************************************
        //    System Configuration
        // *********************************************
        bool are_rec_pc_directly_coupled = true;



        /***********************************************
             Power cycle
         ***********************************************/

         // Steam Rankine and User Defined power cycle classes
        C_pc_Rankine_indirect_224 rankine_pc;
        C_pc_indirect_Gen3 indirect_pc;
        C_csp_power_cycle* p_csp_power_cycle;

        // Select cycle model based on direct or indirect configuration
        if (are_rec_pc_directly_coupled) {
            C_pc_Rankine_indirect_224::S_params* pc = &rankine_pc.ms_params;
            pc->m_P_ref = as_double("P_ref");
            pc->m_eta_ref = as_double("design_eff");
            pc->m_T_htf_hot_ref = as_double("T_pc_hot_des");
            pc->m_T_htf_cold_ref = as_double("T_pc_cold_des");

            pc->m_P_phx_in_co2_des = as_double("P_phx_in_co2_des") / 1000.;  //convert to MPa from kPa
            pc->m_P_turb_in_co2_des = as_double("P_turb_in_co2_des") / 1000.;  //convert to MPa from kPa
            pc->m_P_turb_in_co2_off_sun_des = as_double("P_turb_in_co2_off_sun_des") / 1000.;  //convert to MPa from kPa

            pc->m_cycle_max_frac = as_double("cycle_max_frac");
            pc->m_cycle_cutoff_frac = as_double("cycle_cutoff_frac");
            pc->m_q_sby_frac = as_double("q_sby_frac");
            pc->m_startup_time = as_double("startup_time");
            pc->m_startup_frac = as_double("startup_frac");
            pc->m_htf_pump_coef = 0.0;          // no pc htf pump for co2 system as particle-to-cycle is gravity-fed    as_double("pb_pump_coef");
            pc->m_pc_fl = as_integer("rec_htf");                            // power cycle HTF is same as receiver HTF
            pc->m_pc_fl_props = as_matrix("field_fl_props");

            pc->m_is_user_defined_pc = true;
            pc->m_is_udpc_co2 = as_boolean("is_udpc_co2");

            // User-Defined Cycle Parameters
            pc->m_T_amb_des = as_double("ud_T_amb_des");    //[C]
            pc->m_W_dot_cooling_des = as_double("ud_f_W_dot_cool_des") / 100.0 * as_double("P_ref");  //[MWe]
            pc->m_m_dot_water_des = as_double("ud_m_dot_water_cool_des");       //[kg/s]

            // Also need lower and upper levels for the 3 independent variables...
            pc->m_T_htf_low = as_double("ud_T_htf_low");            //[C]
            pc->m_T_htf_high = as_double("ud_T_htf_high");          //[C]
            pc->m_T_amb_low = as_double("ud_T_amb_low");            //[C]
            pc->m_T_amb_high = as_double("ud_T_amb_high");          //[C]
            pc->m_m_dot_htf_low = as_double("ud_m_dot_htf_low");    //[-]
            pc->m_m_dot_htf_high = as_double("ud_m_dot_htf_high");  //[-]

            // User-Defined Cycle Off-Design Tables 
            pc->mc_T_htf_ind = as_matrix("ud_T_htf_ind_od");
            pc->mc_T_amb_ind = as_matrix("ud_T_amb_ind_od");
            pc->mc_m_dot_htf_ind = as_matrix("ud_m_dot_htf_ind_od");
            pc->mc_combined_ind = as_matrix("ud_ind_od");

            if (is_assigned("ud_ind_od_off_sun"))
            {
                pc->mc_combined_ind_OS = as_matrix("ud_ind_od_off_sun");
            }

            // Set pointer to parent class
            p_csp_power_cycle = &rankine_pc;
        }
        else {
            C_pc_indirect_Gen3::S_params* pc = &indirect_pc.ms_params;
            pc->m_P_ref = as_double("P_ref");
            pc->m_eta_ref = as_double("design_eff");
            pc->m_T_htf_hot_ref = as_double("T_pc_hot_des");
            pc->m_T_htf_cold_ref = as_double("T_pc_cold_des");

            pc->m_cycle_max_frac = as_double("cycle_max_frac");
            pc->m_cycle_cutoff_frac = as_double("cycle_cutoff_frac");
            pc->m_q_sby_frac = as_double("q_sby_frac");
            pc->m_startup_time = as_double("startup_time");
            pc->m_startup_frac = as_double("startup_frac");
            pc->m_htf_pump_coef = 0.0;          // no pc htf pump for co2 system as particle-to-cycle is gravity-fed    as_double("pb_pump_coef");

            // For peaker, use TES htf
            pc->m_pc_fl = as_integer("store_htf");                            // power cycle HTF is same as receiver HTF
            pc->m_pc_fl_props = as_matrix("store_fl_props");

            // Can test with this false so pc code will use Rankine ND performance
            pc->m_is_user_defined_pc = false;

            // User-Defined Cycle Parameters
            pc->m_T_amb_des = as_double("ud_T_amb_des");    //[C]
            pc->m_W_dot_cooling_des = as_double("ud_f_W_dot_cool_des") / 100.0 * as_double("P_ref");  //[MWe]
            pc->m_m_dot_water_des = as_double("ud_m_dot_water_cool_des");       //[kg/s]

            // Also need lower and upper levels for the 3 independent variables...
            pc->m_T_htf_low = as_double("ud_T_htf_low");            //[C]
            pc->m_T_htf_high = as_double("ud_T_htf_high");          //[C]
            pc->m_T_amb_low = as_double("ud_T_amb_low");            //[C]
            pc->m_T_amb_high = as_double("ud_T_amb_high");          //[C]
            pc->m_m_dot_htf_low = as_double("ud_m_dot_htf_low");    //[-]
            pc->m_m_dot_htf_high = as_double("ud_m_dot_htf_high");  //[-]

            // User-Defined Cycle Off-Design Tables 
            pc->mc_T_htf_ind = as_matrix("ud_T_htf_ind_od");
            pc->mc_T_amb_ind = as_matrix("ud_T_amb_ind_od");
            pc->mc_m_dot_htf_ind = as_matrix("ud_m_dot_htf_ind_od");
            pc->mc_combined_ind = as_matrix("ud_ind_od");

            p_csp_power_cycle = &indirect_pc;
        }

        // Set power cycle outputs common to all power cycle technologies
        p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_ETA_THERMAL, allocate("eta", n_steps_fixed), n_steps_fixed);
        p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_Q_DOT_HTF, allocate("q_pb", n_steps_fixed), n_steps_fixed);
        p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_M_DOT_HTF, allocate("m_dot_pc", n_steps_fixed), n_steps_fixed);
        p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_Q_DOT_STARTUP, allocate("q_dot_pc_startup", n_steps_fixed), n_steps_fixed);
        p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_W_DOT, allocate("P_cycle", n_steps_fixed), n_steps_fixed);
        p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_T_HTF_IN, allocate("T_pc_in", n_steps_fixed), n_steps_fixed);
        p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_T_HTF_OUT, allocate("T_pc_out", n_steps_fixed), n_steps_fixed);
        p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_M_DOT_WATER, allocate("m_dot_water_pc", n_steps_fixed), n_steps_fixed);
        p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_T_COND_OUT, allocate("T_cond_out", n_steps_fixed), n_steps_fixed);
        p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_P_PHX_IN, allocate("P_phx_in", n_steps_fixed), n_steps_fixed);
        p_csp_power_cycle->assign(C_pc_Rankine_indirect_224::E_P_PHX_OUT, allocate("P_phx_out", n_steps_fixed), n_steps_fixed);

        //heliostat field class
        C_pt_sf_perf_interp heliostatfield;

        heliostatfield.ms_params.m_p_start = as_double("p_start");      //[kWe-hr] Heliostat startup energy
        heliostatfield.ms_params.m_p_track = as_double("p_track");      //[kWe] Heliostat tracking power
        heliostatfield.ms_params.m_hel_stow_deploy = as_double("hel_stow_deploy");  // N/A
        heliostatfield.ms_params.m_v_wind_max = as_double("v_wind_max");            // N/A
        heliostatfield.ms_params.m_n_flux_x = 3;
        heliostatfield.ms_params.m_n_flux_y = 1;      // sp match

        heliostatfield.ms_params.m_eta_map_aod_format = false; // as_boolean("eta_map_aod_format");

        // Set 'n_flux_x' and 'n_flux_y' here, for now
        assign("n_flux_y", 1);
        assign("n_flux_x", 3);



        //Load the solar field adjustment factors
        sf_adjustment_factors sf_haf(this);
        size_t n_steps_full = weather_reader.m_weather_data_provider->nrecords(); //steps_per_hour * 8760;
        if (!sf_haf.setup((int)n_steps_full))
            throw exec_error("tcsmolten_salt", "failed to setup sf adjustment factors: " + sf_haf.error());
        //allocate array to pass to tcs
        heliostatfield.ms_params.m_sf_adjust.resize(sf_haf.size());
        for (int i = 0; i < sf_haf.size(); i++)
            heliostatfield.ms_params.m_sf_adjust.at(i) = sf_haf(i);

        // Set callback information
        heliostatfield.mf_callback = ssc_cmod_solarpilot_callback;
        heliostatfield.m_cdata = (void*)this;

        C_pt_sf_perf_interp heliostatfield2 = heliostatfield;
        C_pt_sf_perf_interp heliostatfield3 = heliostatfield;

        //set the performance for each heliostat field
        {
            util::matrix_t<double> all_eta = as_matrix("eta_map");
            /*
            format is:
            az      zen     eta1    eta2    eta3    nh1     nh2     nh3
            */

            util::matrix_t<double> eta_sub;
            util::matrix_t<double> nhelio_sub;
            util::matrix_t<double> area_sub;
            eta_sub.resize(all_eta.nrows(), 3);
            nhelio_sub.resize(all_eta.nrows(), 3);
            area_sub.resize(all_eta.nrows(), 3);

            //fill all vectors with the sun positions
            for (size_t j = 0; j < all_eta.nrows(); j++)
            {
                eta_sub.at(j, 0) = all_eta.at(j, 0);
                eta_sub.at(j, 1) = all_eta.at(j, 1);
                nhelio_sub.at(j, 0) = all_eta.at(j, 0);
                nhelio_sub.at(j, 1) = all_eta.at(j, 1);
                area_sub.at(j, 0) = all_eta.at(j, 0);
                area_sub.at(j, 1) = all_eta.at(j, 1);
            }

            //calculate area of single heliostat
            double A_helio = as_number("helio_width") * as_number("helio_height") * as_number("dens_mirror");
            //also keep track of total field area
            double A_sf = 0.;

            //create a vector of all heliostat field pointers
            std::vector<C_pt_sf_perf_interp*> hfs = { &heliostatfield, &heliostatfield2, &heliostatfield3 };

            //populate
            for (size_t i = 0; i < 3; i++)
            {
                for (size_t j = 0; j < all_eta.nrows(); j++)
                {
                    eta_sub.at(j, 2) = all_eta.at(j, 2 + i);
                    nhelio_sub.at(j, 2) = all_eta.at(j, 5 + i);
                    area_sub.at(j, 2) = all_eta.at(j, 5 + i) * A_helio;
                }

                //set matrices for heliostat field 'i'
                hfs.at(i)->ms_params.m_eta_map = eta_sub;
                hfs.at(i)->ms_params.m_nhelio_map = nhelio_sub;
                hfs.at(i)->ms_params.m_area_map = area_sub;
                //dummy flux maps
                hfs.at(i)->ms_params.m_flux_maps.resize_fill(all_eta.nrows(), 4, 0.25);
                A_sf += area_sub.at(0, 2);
            }
            assign("A_sf", A_sf);
        }

        //// *********************************************************
        ////      Now set Type 222 parameters
        //// *********************************************************
        double H_rec = as_double("rec_height");
        double D_rec = as_double("D_rec");
        double A_rec = H_rec * D_rec;


        std::unique_ptr<C_pt_receiver> receiver, receiver2, receiver3;

        std::unique_ptr<C_mspt_receiver> trans_receiver = std::unique_ptr<C_mspt_receiver>(new C_mspt_receiver());    // transient receiver

        trans_receiver->m_w_rec = D_rec;
        trans_receiver->m_h_rec = H_rec;
        trans_receiver->m_field_fl = as_integer("rec_htf");
        trans_receiver->m_field_fl_props = as_matrix("field_fl_props");
        trans_receiver->m_pipe_loss_per_m = as_double("piping_loss");                       //[Wt/m]
        trans_receiver->m_pipe_length_add = as_double("piping_length_const") / N_rec;   //[m]
        trans_receiver->m_pipe_length_mult = as_double("piping_length_mult");       //[-]
        trans_receiver->m_T_salt_hot_target = as_double("T_rec_hot_des");
        trans_receiver->m_hel_stow_deploy = as_double("hel_stow_deploy");

        trans_receiver->m_efficiency_lookup_input = as_matrix("rec_efficiency_lookup");
        trans_receiver->m_pressure_lookup_input = as_matrix("rec_pressure_lookup");

        std::unique_ptr<C_mspt_receiver> trans_receiver2{ new C_mspt_receiver{*trans_receiver} };
        std::unique_ptr<C_mspt_receiver> trans_receiver3{ new C_mspt_receiver{*trans_receiver} };
        receiver = std::move(trans_receiver);
        receiver2 = std::move(trans_receiver2);
        receiver3 = std::move(trans_receiver3);

        receiver->m_h_tower = receiver2->m_h_tower = receiver3->m_h_tower = as_double("h_tower");
        receiver->m_epsilon = receiver2->m_epsilon = receiver3->m_epsilon = std::numeric_limits<double>::quiet_NaN(); // as_double("epsilon");
        receiver->m_T_htf_hot_des = receiver2->m_T_htf_hot_des = receiver3->m_T_htf_hot_des = as_double("T_rec_hot_des");          //[C]
        receiver->m_T_htf_cold_des = receiver2->m_T_htf_cold_des = receiver3->m_T_htf_cold_des = as_double("T_rec_cold_des");      //[C]
        receiver->m_f_rec_min = receiver2->m_f_rec_min = receiver3->m_f_rec_min = as_double("f_rec_min");
        receiver->m_q_rec_des = receiver2->m_q_rec_des = receiver3->m_q_rec_des = (as_double("P_ref") / N_rec) / as_double("design_eff") * as_double("solarm");
        receiver->m_rec_su_delay = receiver2->m_rec_su_delay = receiver3->m_rec_su_delay = as_double("rec_su_delay");
        receiver->m_rec_qf_delay = receiver2->m_rec_qf_delay = receiver3->m_rec_qf_delay = as_double("rec_qf_delay");
        receiver->m_m_dot_htf_max_frac = receiver2->m_m_dot_htf_max_frac = receiver3->m_m_dot_htf_max_frac = as_double("csp.pt.rec.max_oper_frac");
        receiver->m_eta_pump = receiver2->m_eta_pump = receiver3->m_eta_pump = as_double("eta_pump");
        receiver->m_night_recirc = receiver2->m_night_recirc = receiver3->m_night_recirc = 0;
        //receiver->m_P_cold_des = receiver2->m_P_cold_des = receiver3->m_P_cold_des = as_double("P_phx_in_co2_des") - P_drop_low_temp_tes;       // [kPa]

        // Now try to instantiate mspt_collector_receiver
        C_csp_mspt_collector_receiver collector_receiver(heliostatfield, *receiver);
        C_csp_mspt_collector_receiver collector_receiver2(heliostatfield2, *receiver2);
        C_csp_mspt_collector_receiver collector_receiver3(heliostatfield3, *receiver3);

        // Instantiate tower
        std::vector<C_csp_mspt_collector_receiver> collector_receivers{ collector_receiver, collector_receiver2, collector_receiver3 };
        C_csp_tower_collector_receiver tower(collector_receivers);
        tower.m_field_fl = as_integer("rec_htf");
        tower.m_field_fl_props = as_matrix("field_fl_props");
        tower.m_tes_fl = as_integer("store_htf");
        tower.m_tes_fl_props = as_matrix("store_fl_props");
        tower.hx_duty = as_double("P_ref") / as_double("design_eff") * 1.E6;  //[Wt], indiv. HX duty = q_pb_design
        tower.m_dt_hot = as_double("dt_charging");               //[C]
        tower.T_rec_hot_des = as_double("T_rec_hot_des");        //[C]
        tower.T_rec_cold_des = as_double("T_rec_cold_des");      //[C]
        tower.T_hx_cold_des = as_double("T_tes_cold_des");       //[C]
        tower.h_lift = as_double("h_tower") * 1.1;               //[m] Lift height is 10% taller than receiver optical midline
        tower.riser_length = as_double("h_tower") * as_double("piping_length_mult") + as_double("piping_length_const");  //[m]
        tower.riser_diam = as_double("piping_riser_diam");       //[m]
        tower.downcomer_diam = as_double("piping_downcomer_diam");  //[m]
        tower.L_recHX = as_double("L_recHX");                    //[m]
        tower.n_cells_recHX = as_double("n_cells_recHX");        //[-]
        tower.pipe_loss_per_m = as_double("piping_loss");        //[Wt/m]

        tower.m_is_rec_recirc_available = as_boolean("is_rec_recirc_available");        //[-]
        // If indirect system, then must have a recirculator
        if (!are_rec_pc_directly_coupled) {
            tower.m_is_rec_recirc_available = true;
        }

        // Calculate tower inlet pressure using power block outlet pressure and a temporary HX representing the tes low-temp HX
        double T_in_pb = as_double("T_pc_hot_des");                // [C]
        double T_out_pb = as_double("T_pc_cold_des");              // [C]
        double T_avg_pb = 0.5 * (T_in_pb + T_out_pb);              // [C]
        double P_avg_pb = 0.5 * (as_double("P_phx_in_co2_des") + as_double("P_turb_in_co2_des")); // [kPa]
        double Q_pb = as_double("P_ref") / as_double("design_eff") * 1.e3;                      // [kWt]
        sco2Properties sco2_properties;
        double cp_pb = sco2_properties.Cp(T_avg_pb + 273.15, P_avg_pb);                         // [kJ/kg-K]
        double m_dot_pb = Q_pb / (cp_pb * (T_in_pb - T_out_pb));                                // [kg/s]
        double m_dot_hx = m_dot_pb;
        double T_avg_hx = 0.5 * (as_double("T_pc_cold_des") + as_double("T_rec_cold_des"));     // [C]
        double P_in_hx = as_double("P_phx_in_co2_des");                                         // [kPa]
        C_heat_exchanger* tes_lowtemp_hx_for_calcs = new C_heat_exchanger;
        tes_lowtemp_hx_for_calcs->length_ = as_double("L_LTHX");
        tes_lowtemp_hx_for_calcs->n_cells_ = as_double("n_cells_LTHX");
        double P_drop_low_temp_tes = tes_lowtemp_hx_for_calcs->PressureDropFrac(T_avg_hx, m_dot_hx) * P_in_hx; // [kPa]
        tower.P_cold_des = as_double("P_phx_in_co2_des") - P_drop_low_temp_tes;                 // [kPa]
        delete tes_lowtemp_hx_for_calcs;

        // *******************************************************
        // *******************************************************
        // Set receiver outputs
        //float *p_q_thermal_copy = allocate("Q_thermal_123", n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_FIELD_AREA1, allocate("A_sf1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_FIELD_AREA2, allocate("A_sf2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_FIELD_AREA3, allocate("A_sf3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_FIELD_Q_DOT_INC, allocate("q_sf_inc", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_FIELD_ETA_TOT, allocate("eta_field_tot", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_FIELD_ETA_OPT1, allocate("eta_field1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_FIELD_ETA_OPT2, allocate("eta_field2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_FIELD_ETA_OPT3, allocate("eta_field3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_FIELD_ADJUST, allocate("sf_adjust_out", n_steps_fixed), n_steps_fixed);

        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_INC, allocate("q_dot_rec_inc", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_ETA_THERMAL, allocate("eta_therm", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_THERMAL, allocate("Q_thermal", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_M_DOT_HTF, allocate("m_dot_rec", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_STARTUP, allocate("q_startup", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HTF_IN, allocate("T_rec_in", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HTF_OUT, allocate("T_rec_out", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_PIPE_LOSS, allocate("q_piping_losses", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_LOSS, allocate("q_thermal_loss", n_steps_fixed), n_steps_fixed);

        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_P_HEATTRACE, allocate("P_rec_heattrace", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HTF_OUT_END, allocate("T_rec_out_end", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HTF_OUT_MAX, allocate("T_rec_out_max", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HTF_PANEL_OUT_MAX, allocate("T_panel_out_max", n_steps_fixed), n_steps_fixed);

        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_WALL_INLET, allocate("T_wall_rec_inlet", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_WALL_OUTLET, allocate("T_wall_rec_outlet", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_RISER, allocate("T_wall_riser", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_DOWNC, allocate("T_wall_downcomer", n_steps_fixed), n_steps_fixed);

        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_INC1, allocate("q_dot_rec_inc1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_INC2, allocate("q_dot_rec_inc2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_INC3, allocate("q_dot_rec_inc3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_M_DOT_HTF1, allocate("m_dot_rec1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_M_DOT_HTF2, allocate("m_dot_rec2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_M_DOT_HTF3, allocate("m_dot_rec3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HTF_IN1, allocate("T_rec_in1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HTF_IN2, allocate("T_rec_in2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HTF_IN3, allocate("T_rec_in3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HTF_OUT1, allocate("T_rec_out1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HTF_OUT2, allocate("T_rec_out2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HTF_OUT3, allocate("T_rec_out3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_HX1, allocate("q_dot_HX1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_HX2, allocate("q_dot_HX2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_HX3, allocate("q_dot_HX3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_M_DOT_HX1, allocate("m_dot_HX_tes1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_M_DOT_HX2, allocate("m_dot_HX_tes2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_M_DOT_HX3, allocate("m_dot_HX_tes3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HX_OUT1, allocate("T_HX_tes_out1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HX_OUT2, allocate("T_HX_tes_out2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_T_HX_OUT3, allocate("T_HX_tes_out3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_ETA_THERM1, allocate("eta_rec_therm1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_ETA_THERM2, allocate("eta_rec_therm2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_ETA_THERM3, allocate("eta_rec_therm3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_DP_REC1, allocate("dp_rec1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_DP_REC2, allocate("dp_rec2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_DP_REC3, allocate("dp_rec3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_DP_RISER, allocate("dp_riser", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_DP_DOWNCOMER, allocate("dp_downcomer", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_W_DOT_CO2_RECIRC, allocate("W_dot_recirc", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_PARTICLES, allocate("Q_dot_to_particles", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_REC_THERM_1, allocate("Q_dot_rec_therm_1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_REC_THERM_2, allocate("Q_dot_rec_therm_2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_REC_THERM_3, allocate("Q_dot_rec_therm_3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_DP_CO2_HX_1, allocate("dp_CO2_HX_1", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_DP_CO2_HX_2, allocate("dp_CO2_HX_2", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_DP_CO2_HX_3, allocate("dp_CO2_HX_3", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_DOWNCOMER, allocate("q_dot_downcomer", n_steps_fixed), n_steps_fixed);
        tower.mc_reported_outputs.assign(C_csp_tower_collector_receiver::E_Q_DOT_RISER, allocate("q_dot_riser", n_steps_fixed), n_steps_fixed);

        // Thermal energy storage 
        C_csp_two_tank_two_hx_tes storage;
        tower.set_tes(&storage);        // give storage reference to tower
        C_csp_two_tank_two_hx_tes::S_params* tes = &storage.ms_params;

        if (are_rec_pc_directly_coupled) {

            tower.m_is_T_particle_cold_from_tes_ref = true;

            tes->m_is_hx = true;                                   // hardcode = true

            tes->m_field_fl = as_integer("rec_htf");
            tes->m_field_fl_props = as_matrix("field_fl_props");
        }
        else
        {
            tower.m_is_T_particle_cold_from_tes_ref = false;

            tes->m_is_hx = false;                                   // hardcode = true

            tes->m_field_fl = as_integer("store_htf");
            tes->m_field_fl_props = as_matrix("store_fl_props");
        }

        tes->m_tes_fl = as_integer("store_htf");
        tes->m_tes_fl_props = as_matrix("store_fl_props");
        tes->m_W_dot_pc_design = as_double("P_ref");        //[MWe]
        tes->m_eta_pc = as_double("design_eff");                //[-]
        tes->m_solarm = as_double("solarm");
        tes->m_ts_hours = as_double("tshours");
        tes->m_h_tank = as_double("h_tank");
        tes->m_u_tank = as_double("u_tank");
        tes->m_tank_pairs = as_integer("tank_pairs");
        tes->m_hot_tank_Thtr = as_double("hot_tank_Thtr");
        tes->m_hot_tank_max_heat = as_double("hot_tank_max_heat");
        tes->m_cold_tank_Thtr = as_double("cold_tank_Thtr");
        tes->m_cold_tank_max_heat = as_double("cold_tank_max_heat");
        tes->m_dt_hthx = as_double("dt_ht_discharging");        //[C]
        tes->m_dt_lthx = as_double("dt_lt_discharging");        //[C]
        tes->m_T_tes_hot_des = as_double("T_tes_hot_des");
        tes->m_T_tes_warm_des = as_double("T_tes_warm_des");
        tes->m_T_tes_cold_des = as_double("T_tes_cold_des");
        tes->m_T_ht_in_des = as_double("T_rec_cold_des");       //[C]  TES high-temp HX HTF inlet temperature on field side opposite particles
        tes->m_T_lt_in_des = as_double("T_pc_cold_des");        //[C]  TES low-temp HX HTF inlet temperature on power cycle side opposite particles
        tes->m_T_tank_hot_ini = as_double("T_tes_hot_des");
        tes->m_T_tank_cold_ini = as_double("T_tes_cold_des");
        tes->m_h_tank_min = as_double("h_tank_min");
        tes->m_f_V_hot_ini = as_double("csp.pt.tes.init_hot_htf_percent");
        tes->m_htf_pump_coef = 0.0;     // as_double("pb_pump_coef");
        tes->m_tes_pump_coef = 0.0;     // as_double("tes_pump_coef");             //[kWe/kg/s]
        tes->eta_pump = as_double("eta_pump");                  //[-]
        tes->P_avg = 0.5 * (as_double("P_phx_in_co2_des") + as_double("P_turb_in_co2_des"));  //[kPa]
        tes->L_LTHX = as_double("L_LTHX");                 //[m]
        tes->L_HTHX = as_double("L_HTHX");                 //[m]
        tes->n_cells_LTHX = as_double("n_cells_LTHX");     //[-]
        tes->n_cells_HTHX = as_double("n_cells_HTHX");     //[-]
        //storage.L_LTHX = as_double("L_LTHX");                 //[m]
        //storage.L_HTHX = as_double("L_HTHX");                 //[m]
        //storage.n_cells_LTHX = as_double("n_cells_LTHX");     //[-]
        //storage.n_cells_HTHX = as_double("n_cells_HTHX");     //[-]

        // TOU parameters
        C_csp_tou_block_schedules tou;
        C_csp_tou_block_schedules::S_params *tou_params = &tou.ms_params;
        tou_params->mc_csp_ops.mc_weekdays = as_matrix("weekday_schedule");
        tou_params->mc_csp_ops.mc_weekends = as_matrix("weekend_schedule");
        tou_params->mc_pricing.mc_weekdays = as_matrix("dispatch_sched_weekday");
        if (tou_params->mc_pricing.mc_weekdays.ncells() == 1) { tou_params->mc_pricing.mc_weekdays.resize_fill(12, 24, 1.); };
        tou_params->mc_pricing.mc_weekends = as_matrix("dispatch_sched_weekend");
        if (tou_params->mc_pricing.mc_weekends.ncells() == 1) { tou_params->mc_pricing.mc_weekends.resize_fill(12, 24, 1.); };
        tou.mc_dispatch_params.m_dispatch_optimize = as_boolean("is_dispatch");
        tou.mc_dispatch_params.m_is_write_ampl_dat = as_boolean("is_write_ampl_dat");
        tou.mc_dispatch_params.m_is_ampl_engine = as_boolean("is_ampl_engine");
        tou.mc_dispatch_params.m_ampl_data_dir = as_string("ampl_data_dir");
        tou.mc_dispatch_params.m_ampl_exec_call = as_string("ampl_exec_call");
        
        if( tou.mc_dispatch_params.m_dispatch_optimize )
        {
            tou.mc_dispatch_params.m_optimize_frequency = as_integer("disp_frequency");
            tou.mc_dispatch_params.m_disp_steps_per_hour = as_integer("disp_steps_per_hour");
            tou.mc_dispatch_params.m_optimize_horizon = as_integer("disp_horizon");
            tou.mc_dispatch_params.m_max_iterations = as_integer("disp_max_iter");
            tou.mc_dispatch_params.m_solver_timeout = as_double("disp_timeout");
            tou.mc_dispatch_params.m_mip_gap = as_double("disp_mip_gap");
            tou.mc_dispatch_params.m_presolve_type = as_integer("disp_spec_presolve");
            tou.mc_dispatch_params.m_bb_type = as_integer("disp_spec_bb");
            tou.mc_dispatch_params.m_disp_reporting = as_integer("disp_reporting");
            tou.mc_dispatch_params.m_scaling_type = as_integer("disp_spec_scaling");
            tou.mc_dispatch_params.m_disp_time_weighting = as_double("disp_time_weighting");
            tou.mc_dispatch_params.m_rsu_cost = as_double("disp_rsu_cost");
            tou.mc_dispatch_params.m_csu_cost = as_double("disp_csu_cost");
            tou.mc_dispatch_params.m_pen_delta_w = as_double("disp_pen_delta_w");
            tou.mc_dispatch_params.m_disp_inventory_incentive = as_double("disp_inventory_incentive");
            
            tou.mc_dispatch_params.m_q_rec_standby = as_double("q_rec_standby");
            tou.mc_dispatch_params.m_w_rec_ht = as_double("q_rec_heattrace");

            if (as_boolean("is_wlim_series"))
            {
                size_t n_wlim_series = 0;
                ssc_number_t* wlim_series = as_array("wlim_series", &n_wlim_series);
                if (n_wlim_series != n_steps_full)
                    throw exec_error("tcsmolten_salt", "Invalid net electricity generation limit series dimension. Matrix must have "+util::to_string((int)n_steps_full)+" rows.");
                for (size_t i = 0; i < n_steps_full; i++)
                    tou.mc_dispatch_params.m_w_lim_full.at(i) = (double)wlim_series[i];
            }

    
        }
        tou.mc_dispatch_params.m_is_tod_pc_target_also_pc_max = as_boolean("is_tod_pc_target_also_pc_max");
        tou.mc_dispatch_params.m_is_block_dispatch = ! tou.mc_dispatch_params.m_dispatch_optimize;      //mw
        tou.mc_dispatch_params.m_use_rule_1 = true;
        tou.mc_dispatch_params.m_standby_off_buffer = 2.0;
        tou.mc_dispatch_params.m_use_rule_2 = false;
        tou.mc_dispatch_params.m_q_dot_rec_des_mult = -1.23;
        tou.mc_dispatch_params.m_f_q_dot_pc_overwrite = -1.23;

        size_t n_f_turbine = 0;
        ssc_number_t *p_f_turbine = as_array("f_turb_tou_periods", &n_f_turbine);
        tou_params->mc_csp_ops.mvv_tou_arrays[C_block_schedule_csp_ops::TURB_FRAC].resize(n_f_turbine,0.0);
        //tou_params->mv_t_frac.resize(n_f_turbine, 0.0);
        for( size_t i = 0; i < n_f_turbine; i++ )
            tou_params->mc_csp_ops.mvv_tou_arrays[C_block_schedule_csp_ops::TURB_FRAC][i] = (double)p_f_turbine[i];

        bool is_timestep_input = (as_integer("ppa_multiplier_model") == 1);
        tou_params->mc_pricing.mv_is_diurnal = !(is_timestep_input);
        if (is_timestep_input)
        {
            size_t nmultipliers;
            ssc_number_t *multipliers = as_array("dispatch_factors_ts", &nmultipliers);
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE].resize(nmultipliers, 0.0);
            for (size_t ii = 0; ii < nmultipliers; ii++)
                tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][ii] = multipliers[ii];
        }
        else // standard diuranal input
        {
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE].resize(9, 0.0);
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][0] = as_double("dispatch_factor1");
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][1] = as_double("dispatch_factor2");
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][2] = as_double("dispatch_factor3");
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][3] = as_double("dispatch_factor4");
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][4] = as_double("dispatch_factor5");
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][5] = as_double("dispatch_factor6");
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][6] = as_double("dispatch_factor7");
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][7] = as_double("dispatch_factor8");
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][8] = as_double("dispatch_factor9");
        }

        // System parameters
        C_csp_solver::S_csp_system_params system;
        system.m_pb_fixed_par = as_double("pb_fixed_par");
        system.m_bop_par = as_double("bop_par");
        system.m_bop_par_f = as_double("bop_par_f");
        system.m_bop_par_0 = as_double("bop_par_0");
        system.m_bop_par_1 = as_double("bop_par_1");
        system.m_bop_par_2 = as_double("bop_par_2");

        system.are_rec_pc_directly_coupled = are_rec_pc_directly_coupled;

        // Instantiate Solver       
        C_csp_solver csp_solver(weather_reader, 
                        tower, 
                        *p_csp_power_cycle, 
                        storage, 
                        tou, 
                        system,
                        ssc_cmod_update,
                        (void*)(this));


        // Set solver reporting outputs
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TIME_FINAL, allocate("time_hr", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::ERR_M_DOT, allocate("m_dot_balance", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::ERR_Q_DOT, allocate("q_balance", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::N_OP_MODES, allocate("n_op_modes", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::OP_MODE_1, allocate("op_mode_1", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::OP_MODE_2, allocate("op_mode_2", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::OP_MODE_3, allocate("op_mode_3", n_steps_fixed), n_steps_fixed);


        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TOU_PERIOD, allocate("tou_value", n_steps_fixed), n_steps_fixed);            
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::PRICING_MULT, allocate("pricing_mult", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::PC_Q_DOT_SB, allocate("q_dot_pc_sb", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::PC_Q_DOT_MIN, allocate("q_dot_pc_min", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::PC_Q_DOT_TARGET, allocate("q_dot_pc_target", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::PC_Q_DOT_MAX, allocate("q_dot_pc_max", n_steps_fixed), n_steps_fixed);
        
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::CTRL_IS_REC_SU, allocate("is_rec_su_allowed", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::CTRL_IS_PC_SU, allocate("is_pc_su_allowed", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::CTRL_IS_PC_SB, allocate("is_pc_sb_allowed", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::EST_Q_DOT_CR_SU, allocate("q_dot_est_cr_su", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::EST_Q_DOT_CR_ON, allocate("q_dot_est_cr_on", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::EST_Q_DOT_DC, allocate("q_dot_est_tes_dc", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::EST_Q_DOT_CH, allocate("q_dot_est_tes_ch", n_steps_fixed), n_steps_fixed);
        
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::CTRL_OP_MODE_SEQ_A, allocate("operating_modes_a", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::CTRL_OP_MODE_SEQ_B, allocate("operating_modes_b", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::CTRL_OP_MODE_SEQ_C, allocate("operating_modes_c", n_steps_fixed), n_steps_fixed);
        
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_SOLVE_STATE, allocate("disp_solve_state", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_SOLVE_ITER, allocate("disp_solve_iter", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_SOLVE_OBJ, allocate("disp_objective", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_SOLVE_OBJ_RELAX, allocate("disp_obj_relax", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_QSF_EXPECT, allocate("disp_qsf_expected", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_QSFPROD_EXPECT, allocate("disp_qsfprod_expected", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_QSFSU_EXPECT, allocate("disp_qsfsu_expected", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_TES_EXPECT, allocate("disp_tes_expected", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_PCEFF_EXPECT, allocate("disp_pceff_expected", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_SFEFF_EXPECT, allocate("disp_thermeff_expected", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_QPBSU_EXPECT, allocate("disp_qpbsu_expected", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_WPB_EXPECT, allocate("disp_wpb_expected", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_REV_EXPECT, allocate("disp_rev_expected", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_PRES_NCONSTR, allocate("disp_presolve_nconstr", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_PRES_NVAR, allocate("disp_presolve_nvar", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::DISPATCH_SOLVE_TIME, allocate("disp_solve_time", n_steps_fixed), n_steps_fixed);

        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::SOLZEN, allocate("solzen", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::SOLAZ, allocate("solaz", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::BEAM, allocate("beam", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TDRY, allocate("tdry", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TWET, allocate("twet", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::RH, allocate("RH", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::WSPD, allocate("wspd", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::CR_DEFOCUS, allocate("defocus", n_steps_fixed), n_steps_fixed);

        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TES_Q_DOT_LOSS, allocate("tank_losses", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TES_W_DOT_HEATER, allocate("q_heater", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TES_T_HOT, allocate("T_tes_hot", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TES_T_COLD, allocate("T_tes_cold", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TES_Q_DOT_DC, allocate("q_dc_tes", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TES_Q_DOT_CH, allocate("q_ch_tes", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TES_E_CH_STATE, allocate("e_ch_tes", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TES_M_DOT_DC, allocate("m_dot_tes_dc", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TES_M_DOT_CH, allocate("m_dot_tes_ch", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::TES_V_INITIAL, allocate("V_hot_tank_initial", n_steps_fixed), n_steps_fixed);

        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::COL_W_DOT_TRACK, allocate("pparasi", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::CR_W_DOT_PUMP, allocate("P_tower_pump", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::SYS_W_DOT_PUMP, allocate("htf_pump_power", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::PC_W_DOT_COOLING, allocate("P_cooling_tower_tot", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::SYS_W_DOT_FIXED, allocate("P_fixed", n_steps_fixed), n_steps_fixed);
        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::SYS_W_DOT_BOP, allocate("P_plant_balance_tot", n_steps_fixed), n_steps_fixed);

        csp_solver.mc_reported_outputs.assign(C_csp_solver::C_solver_outputs::W_DOT_NET, allocate("P_out_net", n_steps_fixed), n_steps_fixed);



        update("Initialize MSPT model...", 0.0);

        int out_type = -1;
        std::string out_msg = "";
        try
        {
            // Initialize Solver
            csp_solver.init();
        }
        catch( C_csp_exception &csp_exception )
        {
            // Report warning before exiting with error
            while( csp_solver.mc_csp_messages.get_message(&out_type, &out_msg) )
            {
                log(out_msg, out_type);
            }

            throw exec_error("tcsmolten_salt", csp_exception.m_error_message);
        }

        // If no exception, then report messages
        while (csp_solver.mc_csp_messages.get_message(&out_type, &out_msg))
        {
            log(out_msg, out_type);
        }


        //if the pricing schedule is provided as hourly, overwrite the tou schedule
        if( as_boolean("is_dispatch_series") )
        {
            size_t n_dispatch_series;
            ssc_number_t *dispatch_series = as_array("dispatch_series", &n_dispatch_series);

       //     if( n_dispatch_series != n_steps_fixed)
                //throw exec_error("tcsmolten_salt", "Invalid dispatch pricing series dimension. Array length must match number of simulation time steps ("+my_to_string(n_steps_fixed)+").");
                
            //resize the m_hr_tou array
            if( tou_params->mc_pricing.m_hr_tou != 0 )
                delete [] tou_params->mc_pricing.m_hr_tou;
            tou_params->mc_pricing.m_hr_tou = new double[n_steps_fixed];
            //set the tou period as unique for each time step
            for(size_t i=0; i<n_steps_fixed; i++)
                tou_params->mc_pricing.m_hr_tou[i] = i+1;
            //allocate reported arrays
            tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE].resize(n_steps_fixed);
            for( size_t i=0; i<n_steps_fixed; i++)
                tou_params->mc_pricing.mvv_tou_arrays[C_block_schedule_pricing::MULT_PRICE][i] = dispatch_series[i];
        }

        update("Begin timeseries simulation...", 0.0);

        try
        {
            // Simulate !
            csp_solver.Ssimulate(sim_setup);
        }
        catch(C_csp_exception &csp_exception)
        {
            // Report warning before exiting with error
            while( csp_solver.mc_csp_messages.get_message(&out_type, &out_msg) )
            {
                log(out_msg);
            }

            throw exec_error("tcsmolten_salt", csp_exception.m_error_message);
        }

        // If no exception, then report messages
        while (csp_solver.mc_csp_messages.get_message(&out_type, &out_msg))
        {
            log(out_msg, out_type);
        }

        // ******* Re-calculate system costs here ************
        C_mspt_system_costs sys_costs;

        sys_costs.ms_par.A_sf_refl = as_double("A_sf");
        sys_costs.ms_par.site_improv_spec_cost = as_double("site_spec_cost");
        sys_costs.ms_par.heliostat_spec_cost = as_double("heliostat_spec_cost");
        sys_costs.ms_par.heliostat_fixed_cost = as_double("cost_sf_fixed");

        sys_costs.ms_par.h_tower = as_double("h_tower");
        sys_costs.ms_par.h_rec = H_rec;
        sys_costs.ms_par.h_helio = as_double("helio_height");
        sys_costs.ms_par.tower_fixed_cost = as_double("tower_fixed_cost");
        sys_costs.ms_par.tower_cost_scaling_exp = as_double("tower_exp");
        sys_costs.ms_par.foundation_fixed_cost = as_double("foundation_fixed_cost");
        sys_costs.ms_par.foundation_cost_scaling_quadratic = as_double("foundation_cost_scaling_quadratic");
        sys_costs.ms_par.foundation_cost_scaling_linear = as_double("foundation_cost_scaling_linear");
        sys_costs.ms_par.particle_lift_cost = as_double("particle_lift_cost");
        sys_costs.ms_par.riser_and_downcomer_cost = as_double("riser_and_downcomer_cost");

        sys_costs.ms_par.A_rec = A_rec;
        sys_costs.ms_par.rec_ref_cost = as_double("rec_ref_cost");
        sys_costs.ms_par.A_rec_ref = as_double("rec_ref_area");
        sys_costs.ms_par.rec_cost_scaling_exp = as_double("rec_cost_exp");

        sys_costs.ms_par.Q_storage = as_double("P_ref") / as_double("design_eff")*as_double("tshours");
        sys_costs.ms_par.tes_spec_cost = as_double("tes_spec_cost");

        sys_costs.ms_par.W_dot_design = as_double("P_ref");
        sys_costs.ms_par.power_cycle_spec_cost = as_double("plant_spec_cost");

        sys_costs.ms_par.bop_spec_cost = as_double("bop_spec_cost");

        sys_costs.ms_par.fossil_backup_spec_cost = as_double("fossil_spec_cost");

        sys_costs.ms_par.contingency_rate = as_double("contingency_rate");

        //land area
        sys_costs.ms_par.total_land_area = as_double("land_area_base") * as_double("csp.pt.sf.land_overhead_factor") + as_double("csp.pt.sf.fixed_land_area")+ sys_costs.ms_par.radfield_area/4046.86 /*acres/m^2*/ ;
        assign("csp.pt.cost.total_land_area", (ssc_number_t)sys_costs.ms_par.total_land_area);

        sys_costs.ms_par.plant_net_capacity = system_capacity / 1000.0;         //[MWe], convert from kWe
        sys_costs.ms_par.EPC_land_spec_cost = as_double("csp.pt.cost.epc.per_acre");
        sys_costs.ms_par.EPC_land_perc_direct_cost_smaller = as_double("csp.pt.cost.epc.percent.smaller");
        sys_costs.ms_par.EPC_land_perc_direct_cost_larger = as_double("csp.pt.cost.epc.percent.larger");
        sys_costs.ms_par.EPC_land_per_power_cost = as_double("csp.pt.cost.epc.per_watt");
        sys_costs.ms_par.EPC_land_fixed_cost_smaller = as_double("csp.pt.cost.epc.fixed.smaller");
        sys_costs.ms_par.EPC_land_fixed_cost_larger = as_double("csp.pt.cost.epc.fixed.larger");
        sys_costs.ms_par.total_land_spec_cost = as_double("land_spec_cost");
        sys_costs.ms_par.total_land_perc_direct_cost = as_double("csp.pt.cost.plm.percent");
        sys_costs.ms_par.total_land_per_power_cost = as_double("csp.pt.cost.plm.per_watt");
        sys_costs.ms_par.total_land_fixed_cost = as_double("csp.pt.cost.plm.fixed");
        sys_costs.ms_par.sales_tax_basis = as_double("sales_tax_frac");
        sys_costs.ms_par.sales_tax_rate = as_double("sales_tax_rate");

        try
        {
            sys_costs.calculate_costs();
        }
        catch (C_csp_exception &)
        {
            throw exec_error("MSPT system costs", util::format("System cost calculations failed. Check that all inputs are properly defined"));
        }

        // 1.5.2016 twn: financial model needs an updated total_installed_cost, remaining are for reporting only
        assign("total_installed_cost", (ssc_number_t)sys_costs.ms_out.total_installed_cost);

        assign("csp.pt.cost.site_improvements", (ssc_number_t)sys_costs.ms_out.site_improvement_cost);
        assign("csp.pt.cost.heliostats", (ssc_number_t)sys_costs.ms_out.heliostat_cost);
        assign("csp.pt.cost.tower", (ssc_number_t)sys_costs.ms_out.tower_cost);
        assign("csp.pt.cost.receiver", (ssc_number_t)sys_costs.ms_out.receiver_cost);
        assign("csp.pt.cost.storage", (ssc_number_t)sys_costs.ms_out.tes_cost);
        assign("csp.pt.cost.power_block", (ssc_number_t)sys_costs.ms_out.power_cycle_cost);
        
        assign("csp.pt.cost.bop", (ssc_number_t)sys_costs.ms_out.bop_cost);
        assign("csp.pt.cost.fossil", (ssc_number_t)sys_costs.ms_out.fossil_backup_cost);
        assign("ui_direct_subtotal", (ssc_number_t)sys_costs.ms_out.direct_capital_precontingency_cost);
        assign("csp.pt.cost.contingency", (ssc_number_t)sys_costs.ms_out.contingency_cost);
        assign("total_direct_cost", (ssc_number_t)sys_costs.ms_out.total_direct_cost);
        assign("csp.pt.cost.epc.total", (ssc_number_t)sys_costs.ms_out.epc_and_owner_cost);
        assign("csp.pt.cost.plm.total", (ssc_number_t)sys_costs.ms_out.total_land_cost);
        assign("csp.pt.cost.sales_tax.total", (ssc_number_t)sys_costs.ms_out.sales_tax_cost);
        assign("total_indirect_cost", (ssc_number_t)sys_costs.ms_out.total_indirect_cost);
        assign("csp.pt.cost.installed_per_capacity", (ssc_number_t)sys_costs.ms_out.estimated_installed_cost_per_cap);

        // Update construction financing costs, specifically, update: "construction_financing_cost"
        double const_per_interest_rate1 = as_double("const_per_interest_rate1");
        double const_per_interest_rate2 = as_double("const_per_interest_rate2");
        double const_per_interest_rate3 = as_double("const_per_interest_rate3");
        double const_per_interest_rate4 = as_double("const_per_interest_rate4");
        double const_per_interest_rate5 = as_double("const_per_interest_rate5");
        double const_per_months1 = as_double("const_per_months1");
        double const_per_months2 = as_double("const_per_months2");
        double const_per_months3 = as_double("const_per_months3");
        double const_per_months4 = as_double("const_per_months4");
        double const_per_months5 = as_double("const_per_months5");
        double const_per_percent1 = as_double("const_per_percent1");
        double const_per_percent2 = as_double("const_per_percent2");
        double const_per_percent3 = as_double("const_per_percent3");
        double const_per_percent4 = as_double("const_per_percent4");
        double const_per_percent5 = as_double("const_per_percent5");
        double const_per_upfront_rate1 = as_double("const_per_upfront_rate1");
        double const_per_upfront_rate2 = as_double("const_per_upfront_rate2");
        double const_per_upfront_rate3 = as_double("const_per_upfront_rate3");
        double const_per_upfront_rate4 = as_double("const_per_upfront_rate4");
        double const_per_upfront_rate5 = as_double("const_per_upfront_rate5");

        double const_per_principal1, const_per_principal2, const_per_principal3, const_per_principal4, const_per_principal5;
        double const_per_interest1, const_per_interest2, const_per_interest3, const_per_interest4, const_per_interest5;
        double const_per_total1, const_per_total2, const_per_total3, const_per_total4, const_per_total5;
        double const_per_percent_total, const_per_principal_total, const_per_interest_total, construction_financing_cost;

        const_per_principal1 = const_per_principal2 = const_per_principal3 = const_per_principal4 = const_per_principal5 =
            const_per_interest1 = const_per_interest2 = const_per_interest3 = const_per_interest4 = const_per_interest5 =
            const_per_total1 = const_per_total2 = const_per_total3 = const_per_total4 = const_per_total5 =
            const_per_percent_total = const_per_principal_total = const_per_interest_total = construction_financing_cost =
            std::numeric_limits<double>::quiet_NaN();

        N_financial_parameters::construction_financing_total_cost(sys_costs.ms_out.total_installed_cost,
            const_per_interest_rate1, const_per_interest_rate2, const_per_interest_rate3, const_per_interest_rate4, const_per_interest_rate5,
            const_per_months1, const_per_months2, const_per_months3, const_per_months4, const_per_months5,
            const_per_percent1, const_per_percent2, const_per_percent3, const_per_percent4, const_per_percent5,
            const_per_upfront_rate1, const_per_upfront_rate2, const_per_upfront_rate3, const_per_upfront_rate4, const_per_upfront_rate5,
            const_per_principal1, const_per_principal2, const_per_principal3, const_per_principal4, const_per_principal5,
            const_per_interest1, const_per_interest2, const_per_interest3, const_per_interest4, const_per_interest5,
            const_per_total1, const_per_total2, const_per_total3, const_per_total4, const_per_total5,
            const_per_percent_total, const_per_principal_total, const_per_interest_total, construction_financing_cost);

        assign("const_per_principal1", (ssc_number_t)const_per_principal1);
        assign("const_per_principal2", (ssc_number_t)const_per_principal2);
        assign("const_per_principal3", (ssc_number_t)const_per_principal3);
        assign("const_per_principal4", (ssc_number_t)const_per_principal4);
        assign("const_per_principal5", (ssc_number_t)const_per_principal5);
        assign("const_per_interest1", (ssc_number_t)const_per_interest1);
        assign("const_per_interest2", (ssc_number_t)const_per_interest2);
        assign("const_per_interest3", (ssc_number_t)const_per_interest3);
        assign("const_per_interest4", (ssc_number_t)const_per_interest4);
        assign("const_per_interest5", (ssc_number_t)const_per_interest5);
        assign("const_per_total1", (ssc_number_t)const_per_total1);
        assign("const_per_total2", (ssc_number_t)const_per_total2);
        assign("const_per_total3", (ssc_number_t)const_per_total3);
        assign("const_per_total4", (ssc_number_t)const_per_total4);
        assign("const_per_total5", (ssc_number_t)const_per_total5);
        assign("const_per_percent_total", (ssc_number_t)const_per_percent_total);
        assign("const_per_principal_total", (ssc_number_t)const_per_principal_total);
        assign("const_per_interest_total", (ssc_number_t)const_per_interest_total);
        assign("construction_financing_cost", (ssc_number_t)construction_financing_cost);

        // Do unit post-processing here
        double *p_q_pc_startup = allocate("q_pc_startup", n_steps_fixed);
        size_t count_pc_su = 0;
        ssc_number_t *p_q_dot_pc_startup = as_array("q_dot_pc_startup", &count_pc_su);
        if( count_pc_su != n_steps_fixed )
        {
            log("q_dot_pc_startup array is a different length than 'n_steps_fixed'.", SSC_WARNING);
            return;
        }
        for( size_t i = 0; i < n_steps_fixed; i++ )
        {
            p_q_pc_startup[i] = (float)(p_q_dot_pc_startup[i] * (sim_setup.m_report_step / 3600.0));    //[MWh]
        }

        // Convert mass flow rates from [kg/hr] to [kg/s]
        size_t count_m_dot_pc, count_m_dot_rec, count_m_dot_water_pc, count_m_dot_tes_dc, count_m_dot_tes_ch;
        count_m_dot_pc = count_m_dot_rec = count_m_dot_water_pc = count_m_dot_tes_dc = count_m_dot_tes_ch = 0;
        ssc_number_t *p_m_dot_rec = as_array("m_dot_rec", &count_m_dot_rec);
        ssc_number_t *p_m_dot_pc = as_array("m_dot_pc", &count_m_dot_pc);
        ssc_number_t *p_m_dot_water_pc = as_array("m_dot_water_pc", &count_m_dot_water_pc);
        ssc_number_t *p_m_dot_tes_dc = as_array("m_dot_tes_dc", &count_m_dot_tes_dc);
        ssc_number_t *p_m_dot_tes_ch = as_array("m_dot_tes_ch", &count_m_dot_tes_ch);
        if (count_m_dot_rec != n_steps_fixed || count_m_dot_pc != n_steps_fixed || count_m_dot_water_pc != n_steps_fixed
            || count_m_dot_tes_dc != n_steps_fixed || count_m_dot_tes_ch != n_steps_fixed)
        {
            log("At least one m_dot array is a different length than 'n_steps_fixed'.", SSC_WARNING);
            return;
        }
        for (size_t i = 0; i < n_steps_fixed; i++)
        {
            p_m_dot_rec[i] = (ssc_number_t)(p_m_dot_rec[i] / 3600.0);   //[kg/s] convert from kg/hr
            p_m_dot_pc[i] = (ssc_number_t)(p_m_dot_pc[i] / 3600.0);     //[kg/s] convert from kg/hr
            p_m_dot_water_pc[i] = (ssc_number_t)(p_m_dot_water_pc[i] / 3600.0); //[kg/s] convert from kg/hr
            p_m_dot_tes_dc[i] = (ssc_number_t)(p_m_dot_tes_dc[i] / 3600.0);     //[kg/s] convert from kg/hr
            p_m_dot_tes_ch[i] = (ssc_number_t)(p_m_dot_tes_ch[i] / 3600.0);     //[kg/s] convert from kg/hr
        }       

        size_t count;
        ssc_number_t *p_W_dot_net = as_array("P_out_net", &count);
        ssc_number_t *p_time_final_hr = as_array("time_hr", &count);

        // 'adjustment_factors' class stores factors in hourly array, so need to index as such
        adjustment_factors haf(this, "adjust");
        if( !haf.setup() )
            throw exec_error("tcsmolten_salt", "failed to setup adjustment factors: " + haf.error());


        ssc_number_t *p_gen = allocate("gen", count);
        for( size_t i = 0; i < count; i++ )
        {
            size_t hour = (size_t)ceil(p_time_final_hr[i]);
            p_gen[i] = (ssc_number_t)(p_W_dot_net[i] * 1.E3 * haf(hour));           //[kWe]
        }

        accumulate_annual_for_year("gen", "annual_energy", sim_setup.m_report_step / 3600.0, steps_per_hour, 1, n_steps_fixed/steps_per_hour);
        
        accumulate_annual_for_year("P_cycle", "annual_W_cycle_gross", 1000.0*sim_setup.m_report_step / 3600.0, steps_per_hour, 1, n_steps_fixed/steps_per_hour);        //[kWe-hr]
        accumulate_annual_for_year("P_cooling_tower_tot", "annual_W_cooling_tower", 1000.0*sim_setup.m_report_step / 3600.0, steps_per_hour, 1, n_steps_fixed / steps_per_hour);        //[kWe-hr]

        accumulate_annual_for_year("q_dot_rec_inc", "annual_q_rec_inc", sim_setup.m_report_step / 3600.0, steps_per_hour, 1, n_steps_fixed / steps_per_hour);           //[MWt-hr]
        accumulate_annual_for_year("q_thermal_loss", "annual_q_rec_loss", sim_setup.m_report_step / 3600.0, steps_per_hour, 1, n_steps_fixed / steps_per_hour);

        assign("annual_eta_rec_th", (ssc_number_t)(1.0 - as_number("annual_q_rec_loss") / as_number("annual_q_rec_inc")));
        assign("annual_eta_rec_th_incl_refl", (ssc_number_t)(as_number("rec_absorptance")*as_number("annual_eta_rec_th")));

        accumulate_annual_for_year("disp_objective", "disp_objective_ann", 1000.0*sim_setup.m_report_step / 3600.0, steps_per_hour, 1, n_steps_fixed/steps_per_hour);
        accumulate_annual_for_year("disp_solve_iter", "disp_iter_ann", 1000.0*sim_setup.m_report_step / 3600.0, steps_per_hour, 1, n_steps_fixed/steps_per_hour);
        accumulate_annual_for_year("disp_presolve_nconstr", "disp_presolve_nconstr_ann", sim_setup.m_report_step / 3600.0/ as_double("disp_frequency"), steps_per_hour, 1, n_steps_fixed/steps_per_hour);
        accumulate_annual_for_year("disp_presolve_nvar", "disp_presolve_nvar_ann", sim_setup.m_report_step / 3600.0/ as_double("disp_frequency"), steps_per_hour, 1, n_steps_fixed/steps_per_hour);
        accumulate_annual_for_year("disp_solve_time", "disp_solve_time_ann", sim_setup.m_report_step/3600. / as_double("disp_frequency"), steps_per_hour, 1, n_steps_fixed/steps_per_hour );

        // Calculated Outputs
            // First, sum power cycle water consumption timeseries outputs
        accumulate_annual_for_year("m_dot_water_pc", "annual_total_water_use", sim_setup.m_report_step / 1000.0, steps_per_hour, 1, n_steps_fixed/steps_per_hour); //[m^3], convert from kg
            // Then, add water usage from mirror cleaning
        ssc_number_t V_water_cycle = as_number("annual_total_water_use");
        double V_water_mirrors = as_double("water_usage_per_wash") / 1000.0*as_double("A_sf")*as_double("washing_frequency");
        assign("annual_total_water_use", (ssc_number_t)(V_water_cycle + V_water_mirrors));

        ssc_number_t ae = as_number("annual_energy");
        ssc_number_t pg = as_number("annual_W_cycle_gross");
        ssc_number_t convfactor = (pg != 0) ? 100 * ae / pg : (ssc_number_t)0.0;
        assign("conversion_factor", convfactor);

        double kWh_per_kW = 0.0;
        double nameplate = system_capacity;     //[kWe]
        if(nameplate > 0.0)
            kWh_per_kW = ae / nameplate;

        assign("capacity_factor", (ssc_number_t)(kWh_per_kW / ((double)n_steps_fixed / (double)steps_per_hour)*100.));
        assign("kwh_per_kw", (ssc_number_t)kWh_per_kW);
         
        //Single value outputs from radiative cooling system

    }
};

DEFINE_MODULE_ENTRY(tcsmolten_salt, "CSP molten salt power tower with hierarchical controller and dispatch optimization", 1)
