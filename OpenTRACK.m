%% OpenLAP Laptime Simulation Project
%
% OpenTRACK
%
% Track model file creation for use in OpenLAP.
% Instructions:
% 1) Select a track excel file containing the track information by
%    assigning the full path to the variable "filename". To select the
%    input method, change the "mode" variable to 'shape data' or 'logged
%    data'.
%    a) In "shape data" mode use the "OpenTRACK Shape tmp.xlsx" file
%       to create a new track excel file.
%    b) In "logged data" mode, use the "OpenTRACK Logged Data tmp.csv"
%       file to create a new track excel file. Make sure the rows and
%       columns in the functions correspond to the correct channels.
%       Channels needed to generate a usable track are distance, speed and
%       lateral acceleration or yaw velocity. To select between lateral
%       acceleration and yaw velocity for the curvature calculation, set
%       the "log_mode" variable to 'speed & latacc' or 'speed & yaw'.
%       Elevation and banking can be set to 0 everywhere if no data is
%       available. The grip factor should be set to 1 everywhere, and
%       tweaked only increase correlation in specific parts of a track. To
%       filter the logged data, change the duration of the filter from the
%       "filter_dt" variable (a value of 0.5 [s] is recommended).
% 2) Set the meshing size to the desired value (a value of 1 to 5 [m] is
%    recommended).
% 3) Set the track map rotation angle to the desired value in [deg].
%    Zero corresponds to the start of the map pointing towards positive X.
% 4) Run the script.
% 5) The results will appear on the command window and inside the folder
%    "OpenTRACK Tracks".
%
% More information can be found in the "OpenLAP Laptime Simulator"
% videos on YouTube.
%
% This software is licensed under the GPL V3 Open Source License.
%
% Open Source MATLAB project created by:
%
% Michael Chalkiopoulos
% Cranfield University Advanced Motorsport MSc Engineer
% National Technical University of Athens MEng Mechanical Engineer
%
% LinkedIn: https://www.linkedin.com/in/michael-chalkiopoulos/
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

%% Track excel file selection

filename = 'Paul Ricard data.csv' ;

%% Mode

mode = 'logged data' ;
% log_mode = 'speed & yaw' ;
log_mode = 'speed & latacc' ;
% mode = 'shape data' ;

%% Settings

% meshing
mesh_size = 1 ; % [m]
% filtering for logged data mode
filter_dt = 0.5 ; % [s]
% track map rotation angle
rotation = 90+45 ; % [deg]
% track map shape adjuster
lambda = 1 ; % [-]

%% HUD

disp(['Reading track file: ',filename])

%% Reading file
    
if strcmp(mode,'logged data')
    
    %% from logged data
    
    [head,data] = read_logged_data(filename) ;
    info.name = head(2,2) ;
    info.country = head(3,2) ;
    info.city = head(4,2) ;
    info.type = head(5,2) ;
    info.config = head(6,2) ;
    info.direction = head(7,2) ;
    info.mirror = head(8,2) ;
    % channels
    channels = head(11,:) ;
    units = head(12,:) ;
    % frequency
    freq = str2double(head(9,2)) ;
    % data columns
    col_dist = 1 ;
    col_vel = 2 ;
    col_yaw = 3 ;
    col_ay = 4 ;
    col_el = 5 ;
    col_bk = 6 ;
    col_gf = 7 ;
    col_sc = 8 ;
    % extracting data
    x = data(:,col_dist) ;
    v = data(:,col_vel) ;
    w = data(:,col_yaw) ;
    ay = data(:,col_ay) ;
    el = data(:,col_el) ;
    bk = data(:,col_bk) ;
    gf = data(:,col_gf) ;
    sc = data(:,col_sc) ;
    % converting units
    if units(col_dist)~="m"
        switch units(col_dist) % target is [m]
            case 'km'
                x = x*1000 ;
            case 'miles'
                x = x*1609.34 ;
            case 'ft'
                x = x*0.3048 ;
            otherwise
                warning('Check distance units.')
        end
    end
    if units(col_vel)~="m/s"
        switch units(col_vel) % target is [m/s]
            case 'km/h'
                v = v/3.6 ;
            case 'mph'
                v = v*0.44704 ;
            otherwise
                warning('Check speed units.')
        end
    end
    if units(col_yaw)~="rad/s"
        switch units(col_yaw) % target is [rad/s]
            case 'deg/s'
                w = w*2*pi/360 ;
            case 'rpm'
                w = w*2*pi/60 ;
            case 'rps'
                w = w*2*pi ;
            otherwise
                warning('Check yaw velocity units.')
        end
    end
    if units(col_ay)~="m/s/s"
        switch units(col_ay) % target is [m/s2]
            case 'G'
                ay = ay*9.81 ;
            case 'ft/s/s'
                ay = ay*0.3048 ;
            otherwise
                warning('Check lateral acceleration units.')
        end
    end
    if units(col_el)~="m"
        switch units(col_el) % target is [m]
            case 'km'
                el = el*1000 ;
            case 'miles'
                el = el*1609.34 ;
            case 'ft'
                el = el*0.3048 ;
            otherwise
                warning('Check elevation units.')
        end
    end
    if units(col_bk)~="deg"
        switch units(col_bk) % target is [m]
            case 'rad'
                bk = bk/2/pi*360 ;
            otherwise
                warning('Check banking units.')
        end
    end
else
    %% from shape data
    
    [info] = read_info(filename,'Info') ;
    table_shape = read_shape_data(filename,'Shape') ;
    table_el = read_data(filename,'Elevation') ;
    table_bk = read_data(filename,'Banking') ;
    table_gf = read_data(filename,'Grip Factors') ;
    table_sc = read_data(filename,'Sectors') ;
    
end

%% Track name

[folder_status,folder_msg] = mkdir('OpenTRACK Tracks') ;
trackname = "OpenTRACK Tracks/OpenTRACK_"+info.name+"_"+info.config+"_"+info.direction ;
if strcmp(info.mirror,"On")
    trackname = trackname+"_Mirrored" ;
end

%% HUD

delete(trackname+".log") ;
diary(trackname+".log") ;
disp([...
    '_______                    ____________________________________ __';...
    '__  __ \______________________  __/__  __ \__    |_  ____/__  //_/';...
    '_  / / /__  __ \  _ \_  __ \_  /  __  /_/ /_  /| |  /    __  ,<   ';...
    '/ /_/ /__  /_/ /  __/  / / /  /   _  _, _/_  ___ / /___  _  /| |  ';...
    '\____/ _  .___/\___//_/ /_//_/    /_/ |_| /_/  |_\____/  /_/ |_|  ';...
    '       /_/                                                        '...
    ]) ;
disp('==================================================================')
disp(filename)
disp('File read successfully')
disp('==================================================================')
disp("Name:          "+info.name)
disp("City:          "+info.city)
disp("Country:       "+info.country)
disp("Type:          "+info.type)
disp("Configuration: "+info.config)
disp("Direction:     "+info.direction)
disp("Mirror:        "+info.mirror)
disp("Date:          "+datestr(now,'dd/mm/yyyy'))
disp("Time:          "+datestr(now,'HH:MM:SS'))
disp('==================================================================')
disp('Track generation started.')

if strcmp(mode,'logged data')
    %% coarse mesh data
    
    % getting unique points
    [x,rows_to_keep,~] = unique(x) ;
    v = smooth(v(rows_to_keep),round(freq*filter_dt)) ;
    w = smooth(w(rows_to_keep),round(freq*filter_dt)) ;
    ay = smooth(ay(rows_to_keep),round(freq*filter_dt)) ;
    el = smooth(el(rows_to_keep),round(freq*filter_dt)) ;
    bk = smooth(bk(rows_to_keep),round(freq*filter_dt)) ;
    gf = smooth(gf(rows_to_keep),round(freq*filter_dt)) ;
    sc = sc(rows_to_keep) ;
    % curvature
    switch log_mode
        case 'speed & yaw'
            r = w./v ;
        case 'speed & latacc'
            r = ay./v.^2 ;
    end
    r = smooth(r,round(freq*filter_dt)) ;
    % mirroring if needed
    if strcmp(info.mirror,'On')
        r = -r ;
    end
    % old position vector
    xx = x ;
    % track length
    L = x(end) ;
    
    %% Fine mesh data
    
    % new fine position vector
    if floor(L)<L % check for injecting last point
        x = [(0:mesh_size:floor(L))';L] ;
    else
        x = (0:mesh_size:floor(L))' ;
    end
    % distance step vector
    dx = diff(x) ;
    dx = [dx;dx(end)] ;
    % number of mesh points
    n = length(x) ;
    % fine curvature vector
    r = interp1(xx,r,x,'pchip') ;
    % elevation
    Z = interp1(xx,el,x,'linear') ;
    % banking
    bank = interp1(xx,bk,x,'linear') ;
    % sector
    sector = interp1(xx,sc,x,'previous') ;
    % grip factor
    factor_grip = interp1(xx,gf,x,'linear') ;
    % inclination
    incl = -atand((diff(Z)./diff(x))) ;
    incl = [incl;incl(end)] ;
    incl = smooth(incl,round(freq*filter_dt)) ;
    
else
    %% Pre-processing data
    
    % turning radius
    R = table2array(table_shape(:,3)) ;
    R(R==0) = inf ; % correcting straight segment radius
    % segment length
    l = table2array(table_shape(:,2)) ;
    % total length
    L = sum(l) ;
    % segment type
    type_tmp = table2array(table_shape(:,1)) ;
    % memory preallocation
    N = length(l) ;
    type = zeros(N,1) ;
    type(string(type_tmp)=="Straight") = 0 ;
    type(string(type_tmp)=="Left") = 1 ;
    type(string(type_tmp)=="Right") = -1 ;
    if strcmp(info.mirror,'On')
        type = -type ;
    end
    % getting data from tables and ignoring points with x>L
    el = table2array(table_el) ;
    el(el(:,1)>L,:) = [] ;
    bk = table2array(table_bk) ;
    bk(bk(:,1)>L,:) = [] ;
    gf = table2array(table_gf) ;
    gf(gf(:,1)>L,:) = [] ;
    sc = table2array(table_sc) ;
    sc(sc(:,1)>=L,:) = [] ;
    sc = [sc;[L,sc(end,end)]] ;
    % HUD
    disp('Pre-processing completed.')
    
    %% Coarse track meshing
    
    % position vectors
    X = cumsum(l) ; % end position of each segment
    XC = cumsum(l)-l/2 ; % center position of each segment
    j = 1 ; % index
    x = zeros(length(X)+sum(R==inf),1) ; % preallocation
    r = zeros(length(X)+sum(R==inf),1) ; % preallocation
    % editing points
    for i=1:length(X)
        if R(i)==inf % end of straight point injection
            x(j) = X(i)-l(i) ;
            x(j+1) = X(i) ;
            j = j+2 ;
        else % circular segment center
            x(j) = XC(i) ;
            r(j) = type(i)./R(i) ;
            j = j+1 ;
        end
    end
    % saving coarse results
    rr = r ;
    xx = x ;
    % HUD
    disp('Coarse meshing completed.')
    
    %% Fine track meshing
    
    % new fine position vector
    if floor(L)<L % check for injecting last point
        x = [(0:mesh_size:floor(L))';L] ;
    else
        x = (0:mesh_size:floor(L))' ;
    end
    % distance step vector
    dx = diff(x) ;
    dx = [dx;dx(end)] ;
    % number of mesh points
    n = length(x) ;
    % fine curvature vector
    r = interp1(xx,rr,x,'pchip') ;
    % fine turn direction vector
    t = sign(r) ;
    % elevation
    Z = interp1(el(:,1),el(:,2),x,'pchip','extrap') ;
    % banking
    bank = interp1(bk(:,1),bk(:,2),x,'pchip','extrap') ;
    % inclination
    incl = -atand((diff(Z)./diff(x))) ;
    incl = [incl;incl(end)] ;
    % grip factors
    factor_grip = interp1(gf(:,1),gf(:,2),x,'linear','extrap') ;
    % sectors
    sector = interp1(sc(:,1),sc(:,2),x,'previous') ;
    % HUD
    disp("Fine meshing completed with mesh size: "+num2str(mesh_size)+" [m]")
    
end

%% Map

% coordinate vector preallocation
X = zeros(n,1) ;
Y = zeros(n,1) ;
% segment angles
angle_seg = lambda*rad2deg(dx.*r) ;
% heading angles
angle_head = cumsum(angle_seg) ;
if strcmp(info.config,'Closed') % tangency correction for closed track
    dh = [mod(angle_head(end),sign(angle_head(end))*360);angle_head(end)-sign(angle_head(end))*360] ;
    [~,idx] = min(abs(dh)) ;
    dh = dh(idx) ;
    angle_head = angle_head-x/L*dh ;
    angle_seg = [angle_head(1);diff(angle_head)] ;
end
% map points calculation ignoring elevation
for i=2:n
    % previous point
    p = [X(i-1);Y(i-1);0] ;
    % next point
    xyz = rotz(angle_head(i-1))*[dx(i-1);0;0]+p ;
    % saving point coordinates of next point
    X(i) = xyz(1) ;
    Y(i) = xyz(2) ;
end

%% Finding apexes

% finding Apexes
[~,apex] = findpeaks(abs(r)) ;
% correcting corner type
r_apex = r(apex) ;
% HUD
disp('Apex calculation completed.')

%% Map edit

% track direction
if strcmp(info.direction,'Backward')
    x = x(end)-flipud(x) ;
    r = -flipud(r) ;
    apex = length(x)-flipud(apex) ;
    r_apex = -flipud(r_apex) ;
    incl = -flipud(incl) ;
    bank = -flipud(bank) ;
    factor_frip = flipud(factor_grip) ;
    sector = flipud(sector) ;
    X = flipud(X) ;
    Y = flipud(Y) ;
    Z = flipud(Z) ;
end

% track rotation
% rotating track map
xyz = rotz(rotation)*[X';Y';Z'] ;
X = xyz(1,:)' ;
Y = xyz(2,:)' ;
Z = xyz(3,:)' ;
% HUD
disp('Track rotated.')

% distance step vector
dx = diff(x) ;
% closing map if necessary
if strcmp(info.config,'Closed') % closed track
    % HUD
    disp('Closing fine mesh map.')
    % linear correction vectors
    DX = x/L*(X(1)-X(end)) ;
    DY = x/L*(Y(1)-Y(end)) ;
    DZ = x/L*(Z(1)-Z(end)) ;
    db = x/L*(bank(1)-bank(end)) ;
    % adding correction
    X = X+DX ;
    Y = Y+DY ;
    Z = Z+DZ ;
    bank = bank+db ;
    % recalculating inclination
    incl = -atand((diff(Z)./diff(x))) ;
    incl = [incl;(incl(end-1)+incl(1))/2] ;
    % final point injection in distance vector
    dx = [(dx(1)+dx(end))/2;dx] ;
    % HUD
    disp('Fine mesh map closed.')
else % open track
    % final point injection in distance vector
    dx = [dx;dx(end)] ;
end
% smoothing track inclination
incl = smooth(incl) ;
% HUD
disp('Fine mesh map created.')

%% Finish line arrow

% settings
factor_scale = 25 ;
half_angle = 40 ;
% scaling
scale = max([max(X)-min(X);max(Y)-min(Y)])/factor_scale ;
% nondimentional vector from point 2 to point 1
arrow_n = [X(1)-X(2);Y(1)-Y(2);Z(1)-Z(2)]/norm([X(1)-X(2);Y(1)-Y(2);Z(1)-Z(2)]) ;
% first arrow point
arrow_1 = scale*rotz(half_angle)*arrow_n+[X(1);Y(1);Z(1)] ;
% mid arrow point
arrow_c = [X(1);Y(1);Z(1)] ;
% second arrow point
arrow_2 = scale*rotz(-half_angle)*arrow_n+[X(1);Y(1);Z(1)] ;
% arrow vector components
arrow_x = [arrow_1(1);arrow_c(1);arrow_2(1)] ;
arrow_y = [arrow_1(2);arrow_c(2);arrow_2(2)] ;
arrow_z = [arrow_1(3);arrow_c(3);arrow_2(3)] ;
% final arrow matrix
arrow = [arrow_x,arrow_y,arrow_z] ;

%% Ploting Results

% figure
set(0,'units','pixels') ;
SS = get(0,'screensize') ;
H = 900-90 ;
W = 900 ;
Xpos = floor((SS(3)-W)/2) ;
Ypos = floor((SS(4)-H)/2) ;
f = figure('Name',filename,'Position',[Xpos,Ypos,W,H]) ;
figtitle = ["OpenTRACK","Track Name: "+info.name,"Configuration: "+info.config,"Mirror: "+info.mirror,"Date & Time: "+datestr(now,'yyyy/mm/dd')+" "+datestr(now,'HH:MM:SS')] ;
sgtitle(figtitle)

% rows and columns
rows = 5 ;
cols = 2 ;

% 3d map
subplot(rows,cols,[1,3,5,7,9])
title('3D Map')
hold on
grid on
axis equal
axis tight
xlabel('x [m]')
ylabel('y [m]')
scatter3(X,Y,Z,20,sector,'.')
plot3(arrow_x,arrow_y,arrow_z,'k','LineWidth',2)

% curvature
subplot(rows,cols,2)
title('Curvature')
hold on
grid on
xlabel('position [m]')
ylabel('curvature [m^-^1]')
plot(x,r)
scatter(x(apex),r_apex,'.')
xlim([x(1),x(end)])
legend({'curvature','apex'})

% elevation
subplot(rows,cols,4)
title('Elevation')
hold on
grid on
xlabel('position [m]')
ylabel('elevation [m]')
plot(x,Z)
xlim([x(1),x(end)])

% inclination
subplot(rows,cols,6)
title('Inclination')
hold on
grid on
xlabel('position [m]')
ylabel('inclination [deg]')
plot(x,incl)
xlim([x(1),x(end)])

% banking
subplot(rows,cols,8)
title('Banking')
hold on
grid on
xlabel('position [m]')
ylabel('banking [deg]')
plot(x,bank)
xlim([x(1),x(end)])

% grip factors
subplot(rows,cols,10)
title('Grip Factor')
hold on
grid on
xlabel('position [m]')
ylabel('grip factor [-]')
plot(x,factor_grip)
xlim([x(1),x(end)])

% saving plot
savefig(f,trackname+".fig")
% HUD
disp('Plots created and saved.')

%% Saving circuit

% saving
save(trackname+".mat",'info','x','dx','n','r','bank','incl','factor_grip','sector','r_apex','apex','X','Y','Z','arrow')
% HUD
disp('Track generated successfully.')
% diary
diary('off') ;

%% Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [info] = read_info(workbookFile,sheetName,startRow,endRow)
    % Input handling
    % If no sheet is specified, read first sheet
    if nargin == 1 || isempty(sheetName)
        sheetName = 1;
    end
    % If row start and end points are not specified, define defaults
    if nargin <= 3
        startRow = 1;
        endRow = 7;
    end
    % Setup the Import Options
    opts = spreadsheetImportOptions("NumVariables", 2);
    % Specify sheet and range
    opts.Sheet = sheetName;
    opts.DataRange = "A" + startRow(1) + ":B" + endRow(1);
    % Specify column names and types
    opts.VariableNames = ["info", "data"];
    opts.VariableTypes = ["string", "string"];
    opts = setvaropts(opts, [1, 2], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 2], "EmptyFieldRule", "auto");
    % Import the data
    tbl = readtable(workbookFile, opts, "UseExcel", false);
    for idx = 2:length(startRow)
        opts.DataRange = "A" + startRow(idx) + ":B" + endRow(idx);
        tb = readtable(workbookFile, opts, "UseExcel", false);
        tbl = [tbl; tb]; %#ok<AGROW>
    end
    % Convert to output type
    data = tbl.data ;
    info.name = data(1) ;
    info.country = data(2) ;
    info.city = data(3) ;
    info.type = data(4) ;
    info.config = data(5) ;
    info.direction = data(6) ;
    info.mirror = data(7) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tbl] = read_shape_data(workbookFile,sheetName,startRow,endRow)
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
    opts = spreadsheetImportOptions("NumVariables", 3);
    % Specify sheet and range
    opts.Sheet = sheetName;
    opts.DataRange = "A" + startRow(1) + ":C" + endRow(1);
    % Specify column names and types
    opts.VariableNames = ["Type", "SectionLength", "CornerRadius"];
    opts.VariableTypes = ["categorical", "double", "double"];
    opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
    % Setup rules for import
    opts.MissingRule = "omitrow";
    opts = setvaropts(opts, [2, 3], "TreatAsMissing", '');
    % Import the data
    tbl = readtable(workbookFile, opts, "UseExcel", false);
    for idx = 2:length(startRow)
        opts.DataRange = "A" + startRow(idx) + ":C" + endRow(idx);
        tb = readtable(workbookFile, opts, "UseExcel", false);
        tbl = [tbl; tb]; %#ok<AGROW>
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = read_data(workbookFile,sheetName,startRow,endRow)
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
    opts.VariableNames = ["Point", "Data"];
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
function [header,data] = read_logged_data(filename,header_startRow,header_endRow,data_startRow,data_endRow)
    
    % Initialize variables.
    delimiter = ',';
    if nargin<=2
        header_startRow = 1 ;
        header_endRow = 12 ;
        data_startRow = 14 ;
        data_endRow = inf ;
    end

    % Open the text file.
    fileID = fopen(filename,'r');

    % Header array
    % Format for each line of text:
    header_formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
    % Read columns of data according to the format.
    headerArray = textscan(fileID, header_formatSpec, header_endRow(1)-header_startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', header_startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for block=2:length(header_startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, header_formatSpec, header_endRow(block)-header_startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', header_startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for col=1:length(headerArray)
            headerArray{col} = [headerArray{col};dataArrayBlock{col}];
        end
    end
    % Create output variable
    header = [headerArray{1:end-1}];

    % Data array
    % Pointer to start of file
    fseek(fileID,0,'bof') ;
    % Format for each line of text:
    data_formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    % Read columns of data according to the format.
    dataArray = textscan(fileID, data_formatSpec, data_endRow(1)-data_startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', data_startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for block=2:length(data_startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, data_formatSpec, data_endRow(block)-data_startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', data_startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end
    % Create output variable
    data = [dataArray{1:end-1}];

    % Close the text file.
    fclose(fileID);
end