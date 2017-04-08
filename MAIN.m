clear SWE Density

run MeasurementLocations.m  %This program determines the easting and northing of transect measurements
run Import_Density.m        %Imports snow density values
run Import_Transect.m       %Imports transect snow depth and measurement location data
run Import_Zigzag.m         %Imports zigzag snow depth and measurement location data
run Import_SWE.m            %Converts to SWE and condences data
run Import_Topo.m           %Get topo params
run Import_SWE2.m           %Averages obs in each DEM cell (also SWE options structure)
