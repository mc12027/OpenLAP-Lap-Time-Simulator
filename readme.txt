_______                    _____________________     ________             ________           _____ 
__  __ \______________________  /___    |__  __ \    ___  __ \__________________(_)____________  /_
_  / / /__  __ \  _ \_  __ \_  / __  /| |_  /_/ /    __  /_/ /_  ___/  __ \____  /_  _ \  ___/  __/
/ /_/ /__  /_/ /  __/  / / /  /___  ___ |  ____/     _  ____/_  /   / /_/ /___  / /  __/ /__ / /_  
\____/ _  .___/\___//_/ /_//_____/_/  |_/_/          /_/     /_/    \____/___  /  \___/\___/ \__/  
       /_/                                                                /___/                    

===================================================================================================

OpenLAP Laptime Simulation Project.

More information can be found in the "OpenLAP Laptime Simulation Project" videos on YouTube.

This software is licensed under the GPL V3 Open Source License.

Open Source MATLAB project created by:

Michael Chalkiopoulos
Cranfield University MSc Advanced Motorsport Engineer
National Technical University of Athens MEng Mechanical Engineer

LinkedIn: https://www.linkedin.com/in/michael-halkiopoulos/
email: halkiopoulos_michalis@hotmail.com
MATLAB file exchange: https://uk.mathworks.com/matlabcentral/fileexchange/
GitHub: https://github.com/mc12027

April 2020.

===================================================================================================

General OpenLAP usage instructions:

1) Create a vehicle with OpenVEHICLE.
2) Create a track with OpenTRACK.
3) Run straight line simulations with OpenDRAG.
   OpenDRAG needs an OpenVEHICLE vehicle file to run.
4) Run lap time simulations with OpenLAP.
   OpenLAP needs both a vehicle file from OpenVEHICLE and a track file from OpenTRACK to run.

===================================================================================================

Instructions on how to use each module/script:
___________________________________________________________________________________________________
OpenVEHICLE:
Racing vehicle model file creation for use in OpenLAP and OpenDRAG.
Instructions:
1) Select a vehicle excel file containing the vehicles information by
   assigning the full path to the variable "filename". Use the
   "OpenVEHICLE tmp.xlsx" file to create a new vehicle excel file.
2) Run the script.
3) The results will appear on the command window and inside the folder
   "OpenVEHICLE Vehicles".
___________________________________________________________________________________________________
OpenTRACK:
Track model file creation for use in OpenLAP.
Instructions:
1) Select a track excel file containing the track information by
   assigning the full path to the variable "filename". To select the
   input method, change the "mode" variable to 'shape data' or 'logged
   data'.
   a) In "shape data" mode use the "OpenTRACK Shape tmp.xlsx" file
      to create a new track excel file.
   b) In "logged data" mode, use the "OpenTRACK Logged Data tmp.csv"
      file to create a new track excel file. Make sure the rows and
      columns in the functions correspond to the correct channels.
      Channels needed to generate a usable track are distance, speed and
      lateral acceleration or yaw velocity. To select between lateral
      acceleration and yaw velocity for the curvature calculation, set
      the "log_mode" variable to 'speed & latacc' or 'speed & yaw'.
      Elevation and banking can be set to 0 everywhere if no data is
      available. The grip factor should be set to 1 everywhere, and
      tweaked only increase correlation in specific parts of a track. To
      filter the logged data, change the duration of the filter from the
      "filter_dt" variable (a value of 0.5 [s] is recommended).
2) Set the meshing size to the desired value (a value of 1 to 5 [m] is
   recommended).
3) Set the track map rotation angle to the desired value in [deg].
   Zero corresponds to the start of the map pointing towards positive X.
4) Run the script.
5) The results will appear on the command window and inside the folder
   "OpenTRACK Tracks".
___________________________________________________________________________________________________
OpenDRAG:
Straight line acceleration and braking simulation using a simple point
mass model for a racing vehicle.
Instructions:
1) Select a vehicle file created by OpenVEHICLE by assigning the full
   path to the variable "vehiclefile".
2) Run the script.
3) The results will appear on the command window and inside the folder
   "OpenDRAG Sims". You can choose to include the date and time of each
   simulation in the result file name by changing the
   "use_date_time_in_name" variable to true.
___________________________________________________________________________________________________
OpenLAP:
Lap time simulation using a simple point mass model for a racing vehicle.
Instructions:
1) Select a vehicle file created by OpenVEHICLE by assigning the full
   path to the variable "vehiclefile".
2) Select a track file created by OpenTRACK by assigning the full path to
   the variable "trackfile".
3) Select an export frequency in [Hz] by setting the variable "freq" to
   the desired value.
4) Run the script.
5) The results will appear on the command window and inside the folder
   "OpenLAP Sims". You can choose to include the date and time of each
   simulation in the result file name by changing the
   "use_date_time_in_name" variable to true.
