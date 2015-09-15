%% HAWT Parameter file for CACTUS
% List of variables needed to generate a CACTUS *.geom file for a HAWT.

global CREATEGEOM_PATH

CREATEGEOM_PATH         = './path/to/CreateGeom/';                     % relative path from this file's location to the CreateGeom library
HAWT_DAT_TO_GEOM_PATH   = './path/to/hawt_dat_to_geom/'; % relative path from this file's location to the hawt_dat_to_geom library

addpath(make_absolute_filename(CREATEGEOM_PATH))
addpath(make_absolute_filename(HAWT_DAT_TO_GEOM_PATH))

%% Geometry
% Rotor Parameters
rotor_params.num_blades            = 3;			% Number of blades
rotor_params.radius                = 44.2913;	% Rotor radius (ft)
rotor_params.cone_angle            = 0;			% Cone angle (deg)
rotor_params.tilt_angle            = 0;			% Tilt angle (deg)

% Blade Parameters
blade_params.pitch                 = 0;         % pitch (degrees) - positive is LE into the wind
blade_params.eta                   = 0;         % Blade mount point ratio ((distance behind leading edge of the blade mount point) / (root chord)) 

% Discretization ("grid") Parameters
grid_params.node_distribution_type = 'sin';     % Distribution of grid points ([uniform],given,tanh,sin,custom)
grid_params.num_elems              = 4.0;       % Number of line elements (ignored if node_distribution_type == 'custom' or 'given')
grid_params.r_over_R_start         = 0.00;      % r/R of innermost grid node (ignored if node_distribution_type == 'custom' or 'given') 
grid_params.r_over_R_end           = 1.00;      % r/R of outermost grid node (ignored if node_distribution_type == 'custom' or 'given')

% Specify distribution of grid points (only read if node_distribution_type == 'custom')
grid_params.node_distribution      = [0.03889953, 0.04614853, 0.05474838, 0.06495084, 0.07705454, 0.09141378, 0.10844889, 0.12865853, 0.15263426, 0.18107791, 0.21482209, 0.25485454, 0.30234712, 0.35869001, 0.4255325, 0.5048312, 0.59890735, 0.71051475, 0.84292037, 1.];
