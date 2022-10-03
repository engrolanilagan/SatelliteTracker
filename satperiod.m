function [orbitsecs, orbitmins] = satperiod(h)
% satperiod
% Input: Height of satellite from surface of Earth in meters
% Output 1: Orbital period of satellite in seconds
% Output 2: Orbital period in minutes

    G = 6.67384e-11; % Gravitational Constant
    Me = 5.97219e24; % Mass of Earth in kg
    earth_radius = 6371000+(1.036930945906689e+04); % Radius of Earth in meters
    orbit_radius = earth_radius + h;

    orbitsecs = 2*pi*sqrt((orbit_radius^3)/(G*Me));
    orbitmins = orbitsecs/60;
end