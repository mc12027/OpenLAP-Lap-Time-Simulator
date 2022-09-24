%% OpenLAP Laptime Simulation Project
%
% OpenDRAG
%
% Straight line acceleration and braking simulation using a simple point
% mass model for a racing vehicle.
% Instructions:
% 1) Select a vehicle file created by OpenVEHICLE by assigning the full
%    path to the variable "vehiclefile".
% 2) Run the script.
% 3) The results will appear on the command window and inside the folder
%    "OpenDRAG Sims". You can choose to include the date and time of each
%    simulation in the result file name by changing the
%    "use_date_time_in_name" variable to true.
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

%% Clearing memory

clear
clc
close all force
diary('off')
fclose('all') ;

%% Timer start

% total timer start
total_timer = tic ;

%% Loading vehicle

% filename
vehiclefile = 'OpenVEHICLE Vehicles/OpenVEHICLE_Formula 1_Open Wheel.mat' ;

%% Simulation settings

% date and time in simulation name
use_date_time_in_name = false ;
% time step
dt = 1E-3 ;
% maximum simulation time for memory preallocation
t_max = 60 ;
% acceleration sensitivity for drag limitation
ax_sens = 0.05 ; % [m/s2]
% speed traps
speed_trap = [50;100;150;200;250;300;350]/3.6 ;
% track data
bank = 0 ;
incl = 0 ;

%% Vehicle data preprocessing

% loading file
veh = load(vehiclefile) ;
% mass
M = veh.M ;
% gravity constant
g = 9.81 ;
% longitudinal tyre coefficients
dmx = veh.factor_grip*veh.sens_x ;
mux = veh.factor_grip*veh.mu_x ;
Nx = veh.mu_x_M*g ;
% normal load on all wheels
Wz = M*g*cosd(bank)*cosd(incl) ;
% induced weight from banking and inclination
Wy = M*g*sind(bank) ;
Wx = M*g*sind(incl) ;
% ratios
rf = veh.ratio_final ;
rg = veh.ratio_gearbox ;
rp = veh.ratio_primary ;
% tyre radius
Rt = veh.tyre_radius ;
% drivetrain efficiency
np = veh.n_primary ;
ng = veh.n_gearbox ;
nf = veh.n_final ;
% engine curves
rpm_curve = [0;veh.en_speed_curve] ;
torque_curve = veh.factor_power*[veh.en_torque_curve(1);veh.en_torque_curve] ;
% shift points
shift_points = table2array(veh.shifting(:,1)) ;
shift_points = [shift_points;veh.en_speed_curve(end)] ;

%% Acceleration preprocessing

% memory preallocation
N = t_max/dt ;
T = -ones(N,1) ;
X = -ones(N,1) ;
V = -ones(N,1) ;
A = -ones(N,1) ;
RPM = -ones(N,1) ;
TPS = -ones(N,1) ;
BPS = -ones(N,1) ;
GEAR = -ones(N,1) ;
MODE = -ones(N,1) ;
% initial time
t = 0 ;
t_start = 0 ;
% initial distance
x = 0 ;
x_start = 0 ;
% initial velocity
v = 0 ;
% initial accelerartion
a = 0 ;
% initial gears
gear = 1 ;
gear_prev = 1 ;
% shifting condition
shifting = false ;
% initial rpm
rpm = 0 ;
% initial tps
tps = 0 ;
% initial bps
bps = 0 ;
% initial trap number
trap_number = 1 ;
% speed trap checking condition
check_speed_traps = true ;
% iteration number
i = 1 ;

%% HUD display

%  folder
[folder_status,folder_msg] = mkdir('OpenDRAG Sims') ;
% diary
if use_date_time_in_name
    date_time = "_"+datestr(now,'yyyy_mm_dd')+"_"+datestr(now,'HH_MM_SS') ; %#ok<UNRCH>
else
    date_time = "" ;
end
simname = "OpenDRAG Sims/OpenDRAG_"+veh.name+date_time ;
delete(simname+".log") ;
diary(simname+".log") ;
% HUD
disp([...
    '               _______                    ________________________________';...
    '               __  __ \______________________  __ \__  __ \__    |_  ____/';...
    '               _  / / /__  __ \  _ \_  __ \_  / / /_  /_/ /_  /| |  / __  ';...
    '               / /_/ /__  /_/ /  __/  / / /  /_/ /_  _, _/_  ___ / /_/ /  ';...
    '               \____/ _  .___/\___//_/ /_//_____/ /_/ |_| /_/  |_\____/   ';...
    '                      /_/                                                 '...
    ]) ;
disp('=======================================================================================')
disp(['Vehicle: ',char(veh.name)])
disp("Date:    "+datestr(now,'dd/mm/yyyy'))
disp("Time:    "+datestr(now,'HH:MM:SS'))
disp('=======================================================================================')
disp('Acceleration simulation started:')
disp(['Initial Speed: ',num2str(v*3.6),' km/h'])
disp('|_______Comment________|_Speed_|_Accel_|_EnRPM_|_Gear__|_Tabs__|_Xabs__|_Trel__|_Xrel_|')
disp('|______________________|[km/h]_|__[G]__|_[rpm]_|__[#]__|__[s]__|__[m]__|__[s]__|_[m]__|')

%% Acceleration

% acceleration timer start
acceleration_timer = tic ;
while true
    % saving values
    MODE(i) = 1 ;
    T(i) = t ;
    X(i) = x ;
    V(i) = v ;
    A(i) = a ;
    RPM(i) = rpm ;
    TPS(i) = tps ;
    BPS(i) = 0 ;
    GEAR(i) = gear ;
    % checking if rpm limiter is on or if out of memory
    if v>=veh.v_max
        % HUD
        fprintf('Engine speed limited\t')
        hud(v,a,rpm,gear,t,x,t_start,x_start)
        break
    elseif i==N
        % HUD
        disp(['Did not reach maximum speed at time ',num2str(t),' s'])
        break
    end
    % check if drag limited
    if tps==1 && ax+ax_drag<=ax_sens
        % HUD
        fprintf('Drag limited        \t')
        hud(v,a,rpm,gear,t,x,t_start,x_start)
        break
    end
    % checking speed trap
    if check_speed_traps
        % checking if current speed is above trap speed
        if v>=speed_trap(trap_number)
            fprintf('%s%3d %3d%s ','Speed Trap #',trap_number,round(speed_trap(trap_number)*3.6),'km/h')
            hud(v,a,rpm,gear,t,x,t_start,x_start)
            % next speed trap
            trap_number = trap_number+1 ;
            % checking if speed traps are completed
            if trap_number>length(speed_trap)
                check_speed_traps = false ;
            end
        end
    end
    % aero forces
    Aero_Df = 1/2*veh.rho*veh.factor_Cl*veh.Cl*veh.A*v^2 ;
    Aero_Dr = 1/2*veh.rho*veh.factor_Cd*veh.Cd*veh.A*v^2 ;
    % rolling resistance
    Roll_Dr = veh.Cr*(-Aero_Df+Wz) ;
    % normal load on driven wheels
    Wd = (veh.factor_drive*Wz+(-veh.factor_aero*Aero_Df))/veh.driven_wheels ;
    % drag acceleration
    ax_drag = (Aero_Dr+Roll_Dr+Wx)/M ;
    % rpm calculation
    if gear==0 % shifting gears
        rpm = rf*rg(gear_prev)*rp*v/Rt*60/2/pi ;
        rpm_shift = shift_points(gear_prev) ;
    else % gear change finished
        rpm = rf*rg(gear)*rp*v/Rt*60/2/pi ;
        rpm_shift = shift_points(gear) ;
    end
    % checking for gearshifts
    if rpm>=rpm_shift && ~shifting % need to change gears
        if gear==veh.nog % maximum gear number
            % HUD
            fprintf('Engine speed limited\t')
            hud(v,a,rpm,gear,t,x,t_start,x_start)
            break
        else % higher gear available
            % shifting condition
            shifting = true ;
            % shift initialisation time
            t_shift = t ;
            % zeroing  engine acceleration
            ax = 0 ;
            % saving previous gear
            gear_prev = gear ;
            % setting gear to neutral for duration of gearshift
            gear = 0 ;
        end
    elseif shifting % currently shifting gears
        % zeroing  engine acceleration
        ax = 0 ;
        % checking if gearshift duration has passed
        if t-t_shift>veh.shift_time
            % HUD
            fprintf('%s%2d\t','Shifting to gear #',gear_prev+1)
            hud(v,a,rpm,gear_prev+1,t,x,t_start,x_start)
            % shifting condition
            shifting = false ;
            % next gear
            gear = gear_prev+1 ;
        end
    else % no gearshift
        % max long acc available from tyres
        ax_tyre_max_acc = 1/M*(mux+dmx*(Nx-Wd))*Wd*veh.driven_wheels ;
        % getting power limit from engine
        engine_torque = interp1(rpm_curve,torque_curve,rpm) ;
        wheel_torque = engine_torque*rf*rg(gear)*rp*nf*ng*np ;
        ax_power_limit = 1/M*wheel_torque/Rt ;
        % final long acc
        ax = min([ax_power_limit,ax_tyre_max_acc]) ;
    end
    % tps
    tps = ax/ax_power_limit ;
    % longitudinal acceleration
    a = ax+ax_drag ;
    % new position
    x = x+v*dt+1/2*a*dt^2 ;
    % new velocity
    v = v+a*dt ;
    % new time
    t = t+dt ;
    % next iteration
    i = i+1 ;
end
i_acc = i ; % saving acceleration index
% average acceleration
a_acc_ave = v/t ;
disp(['Average acceleration:    ',num2str(a_acc_ave/9.81,'%6.3f'),' [G]'])
disp(['Peak acceleration   :    ',num2str(max(A)/9.81,'%6.3f'),' [G]'])
% acceleration timer
toc(acceleration_timer)

%% Deceleration preprocessing

% deceleration timer start
deceleration_timer = tic ;
% saving time and position of braking start
t_start = t ;
x_start = x ;
% speed trap condition
check_speed_traps = true ;
% active braking speed traps
speed_trap_decel = speed_trap(speed_trap<=v) ;
trap_number = length(speed_trap_decel) ;

%% HUD display

disp('===============================================================================')
disp('Deceleration simulation started:')
disp(['Initial Speed: ',num2str(v*3.6),' [km/h]'])
disp('|_______Comment________|_Speed_|_Accel_|_EnRPM_|_Gear__|_Tabs__|_Xabs__|_Trel__|_Xrel_|')
disp('|______________________|[km/h]_|__[G]__|_[rpm]_|__[#]__|__[s]__|__[m]__|__[s]__|_[m]__|')

%% Deceleration

while true
    % saving values
    MODE(i) = 2 ;
    T(i) = t ;
    X(i) = x ;
    V(i) = v ;
    A(i) = a ;
    RPM(i) = rpm ;
    TPS(i) = 0 ;
    BPS(i) = bps ;
    GEAR(i) = gear ;
    % checking if stopped or if out of memory
    if v<=0
        % zeroing speed
        v = 0 ;
        % HUD
        fprintf('Stopped             \t')
        hud(v,a,rpm,gear,t,x,t_start,x_start)
        break
    elseif i==N
        % HUD
        disp(['Did not stop at time ',num2str(t),' s'])
        break
    end
    % checking speed trap
    if check_speed_traps
        % checking if current speed is under trap speed
        if v<=speed_trap_decel(trap_number)
            % HUD
            fprintf('%s%3d %3d%s ','Speed Trap #',trap_number,round(speed_trap(trap_number)*3.6),'km/h')
            hud(v,a,rpm,gear,t,x,t_start,x_start)
            % next speed trap
            trap_number = trap_number-1 ;
            % checking if speed traps are completed
            if trap_number<1
                check_speed_traps = false ;
            end
        end
    end
    % aero forces
    Aero_Df = 1/2*veh.rho*veh.factor_Cl*veh.Cl*veh.A*v^2 ;
    Aero_Dr = 1/2*veh.rho*veh.factor_Cd*veh.Cd*veh.A*v^2 ;
    % rolling resistance
    Roll_Dr = veh.Cr*(-Aero_Df+Wz) ;
    % drag acceleration
    ax_drag = (Aero_Dr+Roll_Dr+Wx)/M ;
    % gear
    gear = interp1(veh.vehicle_speed,veh.gear,v) ;
    % rpm
    rpm = interp1(veh.vehicle_speed,veh.engine_speed,v) ;
    % max long dec available from tyres
    ax_tyre_max_dec = -1/M*(mux+dmx*(Nx-(Wz-Aero_Df)/4))*(Wz-Aero_Df) ;
    % final long acc
    ax = ax_tyre_max_dec ;
    % brake pressure
    bps = -veh.beta*veh.M*ax ;
    % longitudinal acceleration
    a = ax+ax_drag ;
    % new position
    x = x+v*dt+1/2*a*dt^2 ;
    % new velocity
    v = v+a*dt ;
    % new time
    t = t+dt ;
    % next iteration
    i = i+1 ;
end
% average deceleration
a_dec_ave = V(i_acc)/(t-t_start) ;
disp(['Average deceleration:    ',num2str(a_dec_ave/9.81,'%6.3f'),' [G]'])
disp(['Peak deceleration   :    ',num2str(-min(A)/9.81,'%6.3f'),' [G]'])
% deceleration timer
toc(deceleration_timer)

%% End of simulation

disp('===============================================================================')
disp('Simulation completed successfully.')% total timer
toc(total_timer)

%% Results compression

% getting values to delete
to_delete = T==-1 ;
% deleting values
T(to_delete) = [] ;
X(to_delete) = [] ;
V(to_delete) = [] ;
A(to_delete) = [] ;
RPM(to_delete) = [] ;
TPS(to_delete) = [] ;
BPS(to_delete) = [] ;
GEAR(to_delete) = [] ;
MODE(to_delete) = [] ;

%% Saving results

save(simname+".mat")
diary('off') ;

%% Plots

% figure window
set(0,'units','pixels') ;
SS = get(0,'screensize') ;
H = 900-90 ;
W = 900 ;
Xpos = floor((SS(3)-W)/2) ;
Ypos = floor((SS(4)-H)/2) ;
f = figure('Name','OpenDRAG Simulation Results','Position',[Xpos,Ypos,W,H]) ;
figtitle = ["OpenDRAG","Vehicle: "+veh.name,"Date & Time: "+datestr(now,'yyyy/mm/dd')+" "+datestr(now,'HH:MM:SS')] ;
sgtitle(figtitle)

% rows & columns
row = 7 ;
col = 2 ;
% plot number
i = 0 ;

% distance
i = i+1 ;
subplot(row,col,i)
hold on
grid on
% title('Traveled Distance')
xlabel('Time [s]')
ylabel('Distance [m]')
plot(T,X)
i = i+1 ;

% speed
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Speed')
xlabel('Time [s]')
ylabel('Speed [km/h]')
plot(T,V*3.6)
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Speed')
xlabel('Distance [m]')
ylabel('Speed [km/h]')
plot(X,V*3.6)

% acceleration
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Acceleration')
xlabel('Time [s]')
ylabel('Acceleration [m/s2]')
plot(T,A)
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Acceleration')
xlabel('Distance [m]')
ylabel('Acceleration [m/s2]')
plot(X,A)

% engine speed
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Engine Speed')
xlabel('Time [s]')
ylabel('Engine Speed [rpm]')
plot(T,RPM)
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Engine Speed')
xlabel('Distance [m]')
ylabel('Engine Speed [rpm]')
plot(X,RPM)

% gear
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Selecetd Gear')
xlabel('Time [s]')
ylabel('Gear [-]')
plot(T,GEAR)
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Selecetd Gear')
xlabel('Distance [m]')
ylabel('Gear [-]')
plot(X,GEAR)

% throttle
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Throttle Position')
xlabel('Time [s]')
ylabel('tps [%]')
plot(T,TPS*100)
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Throttle Position')
xlabel('Distance [m]')
ylabel('tps [%]')
plot(X,TPS*100)

% brake
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Brake Pressure')
xlabel('Time [s]')
ylabel('bps [bar]')
plot(T,BPS/10^5)
i = i+1 ;
subplot(row,col,i)
hold on
grid on
title('Brake Pressure')
xlabel('Distance [m]')
ylabel('bps [bar]')
plot(X,BPS/10^5)

savefig(f,simname+".fig")

%% HUD function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = hud(v,a,rpm,gear,t,x,t_start,x_start)
%     disp(['          Speed         : ',num2str(v*3.6),' km/h'])
%     disp(['          Acceleration  : ',num2str(g/9.81),' G'])
%     disp(['          RPM           : ',num2str(rpm)])
%     disp(['          Gear          : ',num2str(gear)])
%     disp(['          Time          : ',num2str(t-t_start),' s'])
%     disp(['          Distance      : ',num2str(x-x_start),' m'])
%     disp(['          Total Time    : ',num2str(t),' s'])
%     disp(['          Total Distance: ',num2str(x),' m'])
    fprintf('%7.2f\t%7.2f\t%7d\t%7d\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n',v*3.6,a/9.81,round(rpm),gear,t,x,t-t_start,x-x_start) ;
end
