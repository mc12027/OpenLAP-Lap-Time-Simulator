%% OpenLAP Laptime Simulation Project
%
% OpenVEHICLE
%
% Racing vehicle model file creation for use in OpenLAP and OpenDRAG.
% Instructions:
% 1) Select a vehicle excel file containing the vehicles information by
%    assigning the full path to the variable "filename". Use the
%    "OpenVEHICLE tmp.xlsx" file to create a new vehicle excel file.
% 2) Run the script.
% 3) The results will appear on the command window and inside the folder
%    "OpenVEHICLE Vehicles".
%
% More information can be found in the "OpenLAP Laptime Simulator"
% videos on YouTube.
%
% This software is licensed under the GPL V3 Open Source License.
%
% Open Source MATLAB project created by:
%
% Michael Halkiopoulos
% Cranfield University Advanced Motorsport MSc Engineer
% National Technical University of Athens MEng Mechanical Engineer
%
% LinkedIn: https://www.linkedin.com/in/michael-halkiopoulos/
% email: halkiopoulos_michalis@hotmail.com
% MATLAB file exchange: https://uk.mathworks.com/matlabcentral/fileexchange/
% GitHub: https://github.com/mc12027
%
% April 2020.

%% Clearing Memory

clear
clc
close all force
diary('off')
fclose('all') ;

%% Vehicle file selection

filename = 'Formula 1.xlsx' ;

%% Reading vehicle file

info = read_info(filename,'Info') ;
data = read_torque_curve(filename,'Torque Curve') ;

%% Getting variables

% info
name = table2array(info(1,2)) ;
type = table2array(info(2,2)) ;
% index
i = 3 ;
% mass
M = str2double(table2array(info(i,2))) ; i = i+1 ; % [kg]
df = str2double(table2array(info(i,2)))/100 ; i = i+1 ; % [-]
% wheelbase
L = str2double(table2array(info(i,2)))/1000 ; i = i+1 ; % [m]
% steering rack ratio
rack = str2double(table2array(info(i,2))) ; i = i+1 ; % [-]
% aerodynamics
Cl = str2double(table2array(info(i,2))) ; i = i+1 ; % [-]
Cd = str2double(table2array(info(i,2))) ; i = i+1 ; % [-]
factor_Cl = str2double(table2array(info(i,2))) ; i = i+1 ; % [-]
factor_Cd = str2double(table2array(info(i,2))) ; i = i+1 ; % [-]
da = str2double(table2array(info(i,2)))/100 ; i = i+1 ; % [-]
A = str2double(table2array(info(i,2))) ; i = i+1 ; % [m2]
rho = str2double(table2array(info(i,2))) ; i = i+1 ; % [kg/m3]
% brakes
br_disc_d = str2double(table2array(info(i,2)))/1000 ; i = i+1 ; % [m]
br_pad_h = str2double(table2array(info(i,2)))/1000 ; i = i+1 ; % [m]
br_pad_mu = str2double(table2array(info(i,2))) ; i = i+1 ; % [m]
br_nop = str2double(table2array(info(i,2))) ; i = i+1 ; % [m]
br_pist_d = str2double(table2array(info(i,2))) ; i = i+1 ; % [m]
br_mast_d = str2double(table2array(info(i,2))) ; i = i+1 ; % [m]
br_ped_r = str2double(table2array(info(i,2))) ; i = i+1 ; % [m]
% tyres
factor_grip = str2double(table2array(info(i,2))) ; i = i+1 ; % [-]
tyre_radius = str2double(table2array(info(i,2)))/1000 ; i = i+1 ; % [m]
Cr = str2double(table2array(info(i,2))) ; i = i+1 ; % [-]
mu_x = str2double(table2array(info(i,2))) ; i = i+1 ; % [-]
mu_x_M = str2double(table2array(info(i,2))) ; i = i+1 ; % [1/kg]
sens_x = str2double(table2array(info(i,2))) ; i = i+1 ; % [-]
mu_y = str2double(table2array(info(i,2))) ; i = i+1 ; % [-]
mu_y_M = str2double(table2array(info(i,2))) ; i = i+1 ; % [1/kg]
sens_y = str2double(table2array(info(i,2))) ; i = i+1 ; % [-]
CF = str2double(table2array(info(i,2))) ; i = i+1 ; % [N/deg]
CR = str2double(table2array(info(i,2))) ; i = i+1 ; % [N/deg]
% engine
factor_power = str2double(table2array(info(i,2))) ; i = i+1 ;
n_thermal = str2double(table2array(info(i,2))) ; i = i+1 ;
fuel_LHV = str2double(table2array(info(i,2))) ; i = i+1 ; % [J/kg]
% drivetrain
drive = table2array(info(i,2)) ; i = i+1 ;
shift_time = str2double(table2array(info(i,2))) ; i = i+1 ; % [s]
n_primary = str2double(table2array(info(i,2))) ; i = i+1 ;
n_final = str2double(table2array(info(i,2))) ; i = i+1 ;
n_gearbox = str2double(table2array(info(i,2))) ; i = i+1 ;
ratio_primary = str2double(table2array(info(i,2))) ; i = i+1 ;
ratio_final = str2double(table2array(info(i,2))) ; i = i+1 ;
ratio_gearbox = str2double(table2array(info(i:end,2))) ;
nog = length(ratio_gearbox) ;

%% HUD

[folder_status,folder_msg] = mkdir('OpenVEHICLE Vehicles') ;
vehname = "OpenVEHICLE Vehicles/OpenVEHICLE_"+name+"_"+type ;
delete(vehname+".log") ;
diary(vehname+".log") ;
disp([...
    '_______                    ___    ________________  ________________________________';...
    '__  __ \_____________________ |  / /__  ____/__  / / /___  _/_  ____/__  /___  ____/';...
    '_  / / /__  __ \  _ \_  __ \_ | / /__  __/  __  /_/ / __  / _  /    __  / __  __/   ';...
    '/ /_/ /__  /_/ /  __/  / / /_ |/ / _  /___  _  __  / __/ /  / /___  _  /___  /___   ';...
    '\____/ _  .___/\___//_/ /_/_____/  /_____/  /_/ /_/  /___/  \____/  /_____/_____/   ';...
    '       /_/                                                                          '...
    ]) ;
disp('====================================================================================')
disp(filename)
disp('File read successfully')
disp('====================================================================================')
disp("Name: "+name)
disp("Type: "+type)
disp("Date: "+datestr(now,'dd/mm/yyyy'))
disp("Time: "+datestr(now,'HH:MM:SS'))
disp('====================================================================================')
disp('Vehicle generation started.')

%% Brake Model

br_pist_a = br_nop*pi*(br_pist_d/1000)^2/4 ; % [m2]
br_mast_a = pi*(br_mast_d/1000)^2/4 ; % [m2]
beta = tyre_radius/(br_disc_d/2-br_pad_h/2)/br_pist_a/br_pad_mu/4 ; % [Pa/N] per wheel
phi = br_mast_a/br_ped_r*2 ; % [-] for both systems
% HUD
disp('Braking model generated successfully.')

%% Steering Model

a = (1-df)*L ; % distance of front axle from center of mass [mm]
b = -df*L ; % distance of rear axle from center of mass [mm]
C = 2*[CF,CF+CR;CF*a,CF*a+CR*b] ; % steering model matrix
% HUD
disp('Steering model generated successfully.')

%% Driveline Model

% fetching engine curves
en_speed_curve = table2array(data(:,1)) ; % [rpm]
en_torque_curve = table2array(data(:,2)) ; % [N*m]
en_power_curve = en_torque_curve.*en_speed_curve*2*pi/60 ; % [W]
% memory preallocation
% wheel speed per gear for every engine speed value
wheel_speed_gear = zeros(length(en_speed_curve),nog) ;
% vehicle speed per gear for every engine speed value
vehicle_speed_gear = zeros(length(en_speed_curve),nog) ;
% wheel torque per gear for every engine speed value
wheel_torque_gear = zeros(length(en_torque_curve),nog) ;
% calculating values for each gear and engine speed
for i=1:nog
    wheel_speed_gear(:,i) = en_speed_curve/ratio_primary/ratio_gearbox(i)/ratio_final ;
    vehicle_speed_gear(:,i) = wheel_speed_gear(:,i)*2*pi/60*tyre_radius ;
    wheel_torque_gear(:,i) = en_torque_curve*ratio_primary*ratio_gearbox(i)*ratio_final*n_primary*n_gearbox*n_final ;
end
% minimum and maximum vehicle speeds
v_min = min(vehicle_speed_gear,[],'all') ;
v_max = max(vehicle_speed_gear,[],'all') ;
% new speed vector for fine meshing
dv = 0.5/3.6 ;
vehicle_speed = linspace(v_min,v_max,(v_max-v_min)/dv)' ;
% memory preallocation
% gear
gear = zeros(length(vehicle_speed),1) ;
% engine tractive force
fx_engine = zeros(length(vehicle_speed),1) ;
% engine tractive force per gear
fx = zeros(length(vehicle_speed),nog) ;
% optimising gear selection and calculating tractive force
for i=1:length(vehicle_speed)
    % going through the gears
    for j=1:nog
        fx(i,j) = interp1(vehicle_speed_gear(:,j),wheel_torque_gear(:,j)/tyre_radius,vehicle_speed(i),'linear',0) ;
    end
    % getting maximum tractive force and gear
    [fx_engine(i),gear(i)] = max(fx(i,:)) ;
end
% adding values for 0 speed to vectors for interpolation purposes at low speeds
vehicle_speed = [0;vehicle_speed] ;
gear = [gear(1);gear] ;
fx_engine = [fx_engine(1);fx_engine] ;
% final vectors
% engine speed
engine_speed = ratio_final*ratio_gearbox(gear)*ratio_primary.*vehicle_speed/tyre_radius*60/2/pi ;
% wheel torque
wheel_torque = fx_engine*tyre_radius ;
% engine torque
engine_torque = wheel_torque/ratio_final./ratio_gearbox(gear)/ratio_primary/n_primary/n_gearbox/n_final ;
% engine power
engine_power = engine_torque.*engine_speed*2*pi/60 ;
% HUD
disp('Driveline model generated successfully.')

%% Shifting Points and Rev Drops

% finding gear changes
gear_change = diff(gear) ; % gear change will appear as 1
% getting speed right before and after gear change
gear_change = logical([gear_change;0]+[0;gear_change]) ;
% getting engine speed at gear change
engine_speed_gear_change = engine_speed(gear_change) ;
% getting shift points
shift_points = engine_speed_gear_change(1:2:length(engine_speed_gear_change)) ;
% getting arrive points
arrive_points = engine_speed_gear_change(2:2:length(engine_speed_gear_change)) ;
% calculating revdrops
rev_drops = shift_points-arrive_points ;
% creating shifting table
rownames = cell(nog-1,1) ;
for i=1:nog-1
    rownames(i) = {[num2str(i,'%d'),'-',num2str(i+1,'%d')]} ;
end
shifting = table(shift_points,arrive_points,rev_drops,'RowNames',rownames) ;
% HUD
disp('Shift points calculated successfully.')

%% Force model

% gravitational constant
g = 9.81 ;
% drive and aero factors
switch drive
    case 'RWD'
        factor_drive = (1-df) ; % weight distribution
        factor_aero = (1-da) ; % aero distribution
        driven_wheels = 2 ; % number of driven wheels
    case 'FWD'
        factor_drive = df ;
        factor_aero = da ;
        driven_wheels = 2 ;
    otherwise % AWD
        factor_drive = 1 ;
        factor_aero = 1 ;
        driven_wheels = 4 ;
end
% Z axis
fz_mass = -M*g ;
fz_aero = 1/2*rho*factor_Cl*Cl*A*vehicle_speed.^2 ;
fz_total = fz_mass+fz_aero ;
fz_tyre = (factor_drive*fz_mass+factor_aero*fz_aero)/driven_wheels ;
% x axis
fx_aero = 1/2*rho*factor_Cd*Cd*A*vehicle_speed.^2 ;
fx_roll = Cr*abs(fz_total) ;
fx_tyre = driven_wheels*(mu_x+sens_x*(mu_x_M*g-abs(fz_tyre))).*abs(fz_tyre) ;
% HUD
disp('Forces calculated successfully.')

%% GGV Map

% track data
bank = 0 ;
incl = 0 ;
% lateral tyre coefficients
dmy = factor_grip*sens_y ;
muy = factor_grip*mu_y ;
Ny = mu_y_M*g ;
% longitudinal tyre coefficients
dmx = factor_grip*sens_x ;
mux = factor_grip*mu_x ;
Nx = mu_x_M*g ;
% normal load on all wheels
Wz = M*g*cosd(bank)*cosd(incl) ;
% induced weight from banking and inclination
Wy = -M*g*sind(bank) ;
Wx = M*g*sind(incl) ;
% speed map vector
dv = 2 ;
v = (0:dv:v_max)' ;
if v(end)~=v_max
    v = [v;v_max] ;
end
% friction ellipse points
N = 45 ;
% map preallocation
GGV = zeros(length(v),2*N-1,3) ;
for i=1:length(v)
    % aero forces
    Aero_Df = 1/2*rho*factor_Cl*Cl*A*v(i)^2 ;
    Aero_Dr = 1/2*rho*factor_Cd*Cd*A*v(i)^2 ;
    % rolling resistance
    Roll_Dr = Cr*abs(-Aero_Df+Wz) ;
    % normal load on driven wheels
    Wd = (factor_drive*Wz+(-factor_aero*Aero_Df))/driven_wheels ;
    % drag acceleration
    ax_drag = (Aero_Dr+Roll_Dr+Wx)/M ;
    % maximum lat acc available from tyres
    ay_max = 1/M*(muy+dmy*(Ny-(Wz-Aero_Df)/4))*(Wz-Aero_Df) ;
    % max long acc available from tyres
    ax_tyre_max_acc = 1/M*(mux+dmx*(Nx-Wd))*Wd*driven_wheels ;
    % max long acc available from tyres
    ax_tyre_max_dec = -1/M*(mux+dmx*(Nx-(Wz-Aero_Df)/4))*(Wz-Aero_Df) ;
    % getting power limit from engine
    ax_power_limit = 1/M*(interp1(vehicle_speed,factor_power*fx_engine,v(i))) ;
    ax_power_limit = ax_power_limit*ones(N,1) ;
    % lat acc vector
    ay = ay_max*cosd(linspace(0,180,N))' ;
    % long acc vector
    ax_tyre_acc = ax_tyre_max_acc*sqrt(1-(ay/ay_max).^2) ; % friction ellipse
    ax_acc = min(ax_tyre_acc,ax_power_limit)+ax_drag ; % limiting by engine power
    ax_dec = ax_tyre_max_dec*sqrt(1-(ay/ay_max).^2)+ax_drag ; % friction ellipse
    % saving GGV map
    GGV(i,:,1) = [ax_acc',ax_dec(2:end)'] ;
    GGV(i,:,2) = [ay',flipud(ay(2:end))'] ;
    GGV(i,:,3) = v(i)*ones(1,2*N-1) ;
end
% HUD
disp('GGV map generated successfully.')

%% Saving vehicle

% saving
save(vehname+".mat")

%% Plot

% figure
set(0,'units','pixels') ;
SS = get(0,'screensize') ;
H = 900-90 ;
W = 900 ;
Xpos = floor((SS(3)-W)/2) ;
Ypos = floor((SS(4)-H)/2) ;
f = figure('Name','Vehicle Model','Position',[Xpos,Ypos,W,H]) ;
sgtitle(name)

% rows and columns
rows = 4 ;
cols = 2 ;

% engine curves
subplot(rows,cols,1)
hold on
title('Engine Curve')
xlabel('Engine Speed [rpm]')
yyaxis left
plot(en_speed_curve,factor_power*en_torque_curve)
ylabel('Engine Torque [Nm]')
grid on
xlim([en_speed_curve(1),en_speed_curve(end)])
yyaxis right
plot(en_speed_curve,factor_power*en_power_curve/745.7)
ylabel('Engine Power [Hp]')

% gearing
subplot(rows,cols,3)
hold on
title('Gearing')
xlabel('Speed [m/s]')
yyaxis left
plot(vehicle_speed,engine_speed)
ylabel('Engine Speed [rpm]')
grid on
xlim([vehicle_speed(1),vehicle_speed(end)])
yyaxis right
plot(vehicle_speed,gear)
ylabel('Gear [-]')
ylim([gear(1)-1,gear(end)+1])

% traction model
subplot(rows,cols,[5,7])
hold on
title('Traction Model')
plot(vehicle_speed,factor_power*fx_engine,'k','LineWidth',4)
plot(vehicle_speed,min([factor_power*fx_engine';fx_tyre']),'r','LineWidth',2)
plot(vehicle_speed,-fx_aero)
plot(vehicle_speed,-fx_roll)
plot(vehicle_speed,fx_tyre)
for i=1:nog
    plot(vehicle_speed(2:end),fx(:,i),'k--')
end
grid on
xlabel('Speed [m/s]')
ylabel('Force [N]')
xlim([vehicle_speed(1),vehicle_speed(end)])
legend({'Engine tractive force','Final tractive force','Aero drag','Rolling resistance','Max tyre tractive force','Engine tractive force per gear'},'Location','southoutside')

% ggv map
subplot(rows,cols,[2,4,6,8])
hold on
title('GGV Map')
surf(GGV(:,:,2),GGV(:,:,1),GGV(:,:,3))
grid on
xlabel('Lat acc [m/s^2]')
ylabel('Long acc [m/s^2]')
zlabel('Speed [m/s]')
% xlim([-ay_max,ay_max])
% ylim([min(GGV(:,:,1),[],'all'),max(GGV(:,:,1),[],'all')])
view(105,5)
set(gca,'DataAspectRatio',[1 1 0.8])
% cbar = colorbar ;
% set(get(cbar,'Title'),'String','Speed [m/s]')

% saving figure
savefig(vehname+".fig")
% HUD
disp('Plots created and saved.')

%% HUD

% HUD
disp('Vehicle generated successfully.')
% diary
diary('off') ;

%% Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = read_torque_curve(workbookFile,sheetName,startRow,endRow)
    % Input handling
    % If no sheet is specified, read first sheet
    if nargin == 1 || isempty(sheetName)
        sheetName = 1;
    end
    % If row start and end points are not specified, define defaults
    if nargin <= 3
        startRow = 2;
        endRow = 10000;
    end
    % Setup the Import Options
    opts = spreadsheetImportOptions("NumVariables", 2);
    % Specify sheet and range
    opts.Sheet = sheetName;
    opts.DataRange = "A" + startRow(1) + ":B" + endRow(1);
    % Specify column names and types
    opts.VariableNames = ["Engine_Speed_rpm", "Torque_Nm"];
    opts.VariableTypes = ["double", "double"];
    % Setup rules for import
    opts.MissingRule = "omitrow";
    opts = setvaropts(opts, [1, 2], "TreatAsMissing", '');
    % Import the data
    data = readtable(workbookFile, opts, "UseExcel", false);
    for idx = 2:length(startRow)
        opts.DataRange = "A" + startRow(idx) + ":B" + endRow(idx);
        tb = readtable(workbookFile, opts, "UseExcel", false);
        data = [data; tb]; %#ok<AGROW>
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = read_info(workbookFile,sheetName,startRow,endRow)
    % Input handling
    % If no sheet is specified, read first sheet
    if nargin == 1 || isempty(sheetName)
        sheetName = 1;
    end
    % If row start and end points are not specified, define defaults
    if nargin <= 3
        startRow = 2;
        endRow = 10000;
    end
    % Setup the Import Options
    opts = spreadsheetImportOptions("NumVariables", 2);
    % Specify sheet and range
    opts.Sheet = sheetName;
    opts.DataRange = "B" + startRow(1) + ":C" + endRow(1);
    % Specify column names and types
    opts.VariableNames = ["Variable", "Value"];
    opts.VariableTypes = ["string", "string"];
    % Setup rules for import
    opts.MissingRule = "omitrow";
    opts = setvaropts(opts, [1, 2], "TreatAsMissing", '');
    % Import the data
    data = readtable(workbookFile, opts, "UseExcel", false);
    for idx = 2:length(startRow)
        opts.DataRange = "A" + startRow(idx) + ":B" + endRow(idx);
        tb = readtable(workbookFile, opts, "UseExcel", false);
        data = [data; tb]; %#ok<AGROW>
    end
end
