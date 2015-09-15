% create_hawt.m : Reads in parameters from an aero_schedule.dat & hawt_params.m file, builds a HAWT geometry structure
%                 and writes out to a CACTUS-formatted .geom file.
%
%                     Usage : octave create_hawt.m [aero_schedule.dat] [hawt_params_file.m] [output_file.geom]
%
%                 The aerodynamic schedule file (tab-delimited: r/R, c/R, twist, airfoil_id) defining all blades
%                 is specified in aero-schedule.dat.

clc
clear all
close all

%% Read in command line arguments %%
arg_list = argv();

% catch errors
if size(arg_list) < 3
    disp 'Error: Not enough command-line arguments.'
    disp '    Usage : octave create_hawt.m [aero_schedule.dat] [hawt_params_file.m] [output_file.geom]'
    return
end

%% Parameters %%
dat_filename         = arg_list{1};
hawt_params_filename = arg_list{2};
geom_filename        = arg_list{3};

%% Import variables from the hawt_params file
% hawt_params holds three structs:
%       rotor_params
%       blade_params
%       grid_params
run(hawt_params_filename);

%% Generate turbine geometry and write to file
T = hawt_dat_to_geom(dat_filename, geom_filename, rotor_params, blade_params, grid_params);