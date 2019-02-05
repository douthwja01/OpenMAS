function [km nmi mi] = haversine(loc1, loc2)
% HAVERSINE     Compute distance between locations using Haversine formula
%   KM = HAVERSINE(LOC1, LOC2) returns the distance KM in km between
%   locations LOC1 and LOC2 using the Haversine formula.  LOC1 and LOC2 are
%   latitude and longitude coordinates that can be expressed as either
%   strings representing degrees, minutes, and seconds (suffixed with
%   N/S/E/W), or numeric arrays representing decimal degrees (where
%   negative indicates West/South).
%
%   [KM, NMI, MI] = HAVERSINE(LOC1, LOC2) returns the computed distance in
%   kilometers (KM), nautical miles (NMI), and miles (MI).
%
%   Examples
%       haversine('53 08 50N, 001 50 58W', '52 12 16N, 000 08 26E') returns
%           170.2547
%       haversine([53.1472 -1.8494], '52 12.16N, 000 08.26E') returns
%           170.2508
%       haversine([53.1472 -1.8494], [52.2044 0.1406]) returns 170.2563
%
%   Inputs
%       LOC must be either a string specifying the location in degrees,
%       minutes and seconds, or a 2-valued numeric array specifying the
%       location in decimal degrees.  If providing a string, the latitude
%       and longitude must be separated by a comma.
%
%       The first element indicates the latitude while the second is the
%       longitude.
%
%   Notes
%       The Haversine formula is used to calculate the great-circle
%       distance between two points, which is the shortest distance over
%       the earth's surface.
%
%       This program was created using equations found on the website
%       http://www.movable-type.co.uk/scripts/latlong.html

% Created by Josiah Renfree
% May 27, 2010


%% Check user inputs

% If two inputs are given, display error
if ~isequal(nargin, 2)
    error('User must supply two location inputs')
    
% If two inputs are given, handle data
else
    
    locs = {loc1 loc2};     % Combine inputs to make checking easier
    
    % Cycle through to check both inputs
    for i = 1:length(locs)
                
        % Check inputs and convert to decimal if needed
        if ischar(locs{i})
            
            % Parse lat and long info from current input
            temp = regexp(locs{i}, ',', 'split');
            lat = temp{1}; lon = temp{2};
            clear temp
            locs{i} = [];           % Remove string to make room for array
            
            % Obtain degrees, minutes, seconds, and hemisphere
            temp = regexp(lat, '(\d+)\D+(\d+)\D+(\d+)(\w?)', 'tokens');
            temp = temp{1};
            
            % Calculate latitude in decimal degrees
            locs{i}(1) = str2double(temp{1}) + str2double(temp{2})/60 + ...
                str2double(temp{3})/3600;
            
            % Make sure hemisphere was given
            if isempty(temp{4})
                error('No hemisphere given')

            % If latitude is south, make decimal negative
            elseif strcmpi(temp{4}, 'S')
                locs{i}(1) = -locs{i}(1);
            end
            
            clear temp

            % Obtain degrees, minutes, seconds, and hemisphere
            temp = regexp(lon, '(\d+)\D+(\d+)\D+(\d+)(\w?)', 'tokens');
            temp = temp{1};
            
            % Calculate longitude in decimal degrees
            locs{i}(2) = str2double(temp{1}) + str2double(temp{2})/60 + ...
                str2double(temp{3})/3600;
            
            % Make sure hemisphere was given
            if isempty(temp{4})
                error('No hemisphere given')
                
            % If longitude is west, make decimal negative
            elseif strcmpi(temp{4}, 'W')
                locs{i}(2) = -locs{i}(2);
            end
            
            clear temp lat lon
        end
    end
end

% Check that both cells are a 2-valued array
if any(cellfun(@(x) ~isequal(length(x),2), locs))
    error('Incorrect number of input coordinates')
end

% Convert all decimal degrees to radians
locs = cellfun(@(x) x .* pi./180, locs, 'UniformOutput', 0);


%% Begin calculation

R = 6371;                                   % Earth's radius in km
delta_lat = locs{2}(1) - locs{1}(1);        % difference in latitude
delta_lon = locs{2}(2) - locs{1}(2);        % difference in longitude
a = sin(delta_lat/2)^2 + cos(locs{1}(1)) * cos(locs{2}(1)) * ...
    sin(delta_lon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
km = R * c;                                 % distance in km


%% Convert result to nautical miles and miles

nmi = km * 0.539956803;                     % nautical miles
mi = km * 0.621371192;                      % miles
end