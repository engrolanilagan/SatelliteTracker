% ME505 Antennas in Mobile and Satellite Communications Systems
% Group Assignment - MIT214154 MIT221671 MIT222035
% Design of a satellite tracking system
% with ground station at Melbourne city
%
% The program takes as input 4 orbital parameters:
%   Orbit height over ground Inclination ℎ
%   Inclination 𝑖
%   Right Ascension of the Ascending Node (RAAN)
%   Angle of Perigee

clc;
clear all;
close all;
format shortg;

%% Take orbital parameters and simulation time
orbitHeight = input('Enter orbit height with respect to ground in meters:');
inclination = input('Enter angle of inclination in degrees:');
raanAngle = input('Enter angle of RAAN in degrees:');
perigeeAngle = input('Enter angle of perigee in degrees:');
startTime = input('Enter simulation start date as YYYY-MM-DD (default today): ', 's');
duration = input('Enter simulation duration in hours (default 24h): ');

% Input Checks
if isempty(orbitHeight) || isempty(inclination) || isempty(raanAngle) ||...
    isempty(perigeeAngle)
    error('All orbital parameters required');
end

if isempty(startTime)
    startTime = datetime('today');
end

if isempty(duration)
    duration = 24;
end

if orbitHeight <= 0 || duration <= 0
    error('Orbit Height and simulation duration must be greater than zero.');
end

%% Get Orbit Radius and Period
earthRadius = 6371000; % Radius of Earth in meters
orbitRadius = earthRadius + orbitHeight;

% Display Orbit Period
orbitPeriod = satperiod(orbitHeight);
disp('Satellite Orbit Period (seconds):');
disp(orbitPeriod);

%% Create Ground Station at MIT
gsLat = -37.81021150347838;
gsLon = 144.96094289809363;
gsPos = [gsLat gsLon 0];

%% Create Satellite Scenario
startTime = datetime(startTime);
stopTime = startTime + hours(duration);
sampleTime = 60; % Sampling time of 1 minute
sc = satelliteScenario(startTime,stopTime,sampleTime);

%% Create Satellite
% The Semi-major Axis is the orbital radius for circular orbits
semiMajorAxis = orbitRadius;

% Assume circular orbit
eccentricity = 0;

% No idea what this is yet
trueAnomaly = 0;

% Create Satellite Object using Input Arguments
sat = satellite(sc, semiMajorAxis, eccentricity, inclination, ...
        raanAngle, perigeeAngle, trueAnomaly);


%% Find the latitude φ and longitude λ of the satellite as a function of t
[pos, velocity] = states(sat, 'CoordinateFrame','geographic');
latitude = pos(1,:,1);
longitude = pos(2,:,1);

%% Plot the location over ground on a world map
% Get the satellite's index at -180 degrees
i = 1;
startIndex = 1;
endIndex = length(longitude);

while(i < length(longitude))
    if(sign(longitude(i))==1 &&  sign(longitude(i+1))~=1)
        startIndex = i+1;
        break;
    end
    i = i+1;
end

% Plot the location for the whole simulation period
geoplot(gsLat, gsLon, '^', 'LineWidth', 2);
hold ON;
i = startIndex;
currStart = startIndex;
while(i < length(longitude))
    if(sign(longitude(i))==1 &&  sign(longitude(i+1))~=1)
        geoplot(latitude(currStart:i), longitude(currStart:i),...
        'LineWidth', 2);
        if currStart == startIndex
            endIndex = i;
        end
        currStart = i+1;
    end
    i = i+1;
end

title('Satellite Location over World Map');

%% Plot the location in space over a globe in a 3D view
uif = uifigure;
g = geoglobe(uif);
geoplot3(g, latitude, longitude, orbitHeight, 'm', 'LineWidth',2);

%% Find the elevation and azimuth of one exemplary pass with respect to MIT GS
s = referenceSphere('Earth');
lookIter = 1;
for i=startIndex:endIndex
    [satX, satY, satZ] = geodetic2ecef(s,latitude(i),longitude(i), orbitHeight);
    satPos = [satX, satY, satZ];
    [azimIter, elevIter, visiIter] = lookangles(gsPos, satPos);
    lookAngle(lookIter, :) = [azimIter, elevIter, visiIter, i];
    lookIter = lookIter + 1;
end

% If satellite will never be visibile, end program
if sum(lookAngle(:, 3)) == 0
    error('No visibility between satellite and ground station.');
end

% There should be a maximum of 2 reversals for the visibility plot over
% the world map
azCol = 1;
elCol = 2;
visCol = 3; % Visibility Column in Switch Index matrix
minCol = 4;

previousVis = lookAngle(1, visCol);
switchIter = 1;
for i=2:length(lookAngle)
    if previousVis ~= lookAngle(i, visCol)
        switchIndex(switchIter, :) = [lookAngle(i, minCol) i];
        switchIter = switchIter + 1;
    end
    previousVis = lookAngle(i, visCol);
end

% Plot the Visibility Window on First Orbit Period of Satellite
figure()
geoplot(gsLat, gsLon, '^', 'LineWidth', 2);
hold ON;
startIter = startIndex;
for i=1:(length(switchIndex)+1)

    if i ~= (length(switchIndex)+1)
        endIter = switchIndex(i, 1);
        isVis = lookAngle(switchIndex(i, 2)-1, visCol);
    else
        endIter = endIndex;
        isVis = not(isVis);
        
    end

    if isVis == 0
        geoplot(latitude(startIter:endIter),...
            longitude(startIter:endIter),...
                '-r',...
                'LineWidth', 2);
    else
        geoplot(latitude(startIter:endIter),...
            longitude(startIter:endIter),...
                '-b', ...
                'LineWidth', 2);
    end

    if i <= length(switchIndex)
        startIter = switchIndex(i, 1);
    end
end
title('Satellite Orbit Visibil');

%% Get the satellite visibility timeslots
previousVis = 1;
endTimeIter = 1;
startTimeIter = 1;
started = 0;
for i=1:length(longitude)
    [satX, satY, satZ] = geodetic2ecef(s,latitude(i),longitude(i), orbitHeight);
    satPos = [satX, satY, satZ];
    [azimIter, elevIter, visiIter] = lookangles(gsPos, satPos, 0);
    gsValues(i, :) = [azimIter, elevIter, visiIter];
    if previousVis == 0 && visiIter == 1
        startTimes(startTimeIter) = (startTime) + minutes(i);
        startTimesIndex(startTimeIter) = i;
        startTimeIter = startTimeIter + 1;
        started = 1;
    elseif previousVis == 1 && visiIter == 0 && started == 1
        endTimes(endTimeIter) = (startTime) + minutes(i);
        endTimesIndex(endTimeIter) = i;
        endTimeIter = endTimeIter + 1;
        started = 0;
    end
    previousVis = visiIter;
end

durations = endTimes - startTimes;
rowNames = {'Start Time'; 'End Time'; 'Duration'};
timeTable = table(startTimes', endTimes', durations');
timeTable.Properties.VariableNames = {'Start Time' 'End Time' 'Duration'};
disp('Satellite Visibility Time Slots:');
disp(timeTable);

%% Plot the distance, azimuth and elevation for examplary overpass of satellite
visStartIndex = (startTimes(1)-startTime)/minutes(1);
visEndIndex = visStartIndex + (endTimes(1)-startTimes(1))/minutes(1) - 1;

timeAxis = startTimes(1):minutes(1):(endTimes(1)-minutes(1));

% Plot the Azimuth Angle during Overpass
figure();
plot(timeAxis, gsValues(visStartIndex:visEndIndex, 1), '-b');
grid ON;
xlabel('Time (hh:mm)');
ylabel('Angle of Azimuth');
title('Azimuth Angle of Satellite during Overpass');

% Plot the Elevation Angle during Overpass
figure();
plot(timeAxis, gsValues(visStartIndex:visEndIndex, 2), '-g');
grid ON;
xlabel('Time (hh:mm)');
ylabel('Angle of Elevation');
title('Elevation Angle of Satellite during Overpass');

% Plot the distance of  satellite to the GS during Overpass
[distX, distY, distZ] = ecefOffset(s, latitude(visStartIndex:visEndIndex),...
    longitude(visStartIndex:visEndIndex), orbitHeight, gsLat, gsLon, 0);
for i=1:length(distX)
    distN(i) = norm([distX(i) distY(i) distZ(i)]);
end

figure();
plot(timeAxis, distN, '-r');
grid ON;
xlabel('Time (hh:mm)');
ylabel('Distance (meters');
title('Distance between Satellite and Ground Station during Overpass');


for i=1:length(longitude)
    timestring(i) = string(datetime(startTime + minutes(i)));

end
rowNames = {'Time'; 'Azimuth'; 'Elevation'};
timeAzimElev = table(timestring', gsValues(:, 1), gsValues(:, 2));
timeAzimElev.Properties.VariableNames = {'Time' 'Azimuth' 'Elevation'};
disp('Satellite Azimuth and Elevation:');
disp(timeAzimElev);
