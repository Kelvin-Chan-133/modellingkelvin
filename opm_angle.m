Dr = 0.46;             % Rim diameter in meters (1.5 ft)
Db = 0.24;            % Ball diameter in meters (0.8 ft)
l  = 7.24;             % Horizontal distance in meters 
                         % (23 ft., 9 in)
playerh = 1.88;           % Height of player (in meters)
releaseh = 1.25*playerh; % Vertical position of center of 
                         % ball at release
basketh  = 3.05;        % Vertical position of the ring 
                         % (in meters), 10ft.
h  = basketh - releaseh; % Vertical distance traversed in 
                         % meters
h  = 0.7;            % Vertical distance traversed in 
                         % meters (2 ft 3.6 in)
g  = -9.8;              % Gravity constant in meters per 
                         % second squared (-32 ft (s)^(-2))
convfactor = 1/0.3048;   % Factor by which measurements in 
                         % m have to be multiplied to get feet.

function [v0, T0] = calcOptv(theta0)
    global l h g   
    % Function to calculate the optimal velocity to hit the center
    % of the basket, given an initial angle of theta0
    v0 = l / cos(theta0 * pi / 180) * sqrt(-g / (2 * (l * tan(theta0 * pi / 180) - h)));
    T0 = l / (cos(theta0 * pi / 180) * v0);
end


function [x,y] = position(v,theta,t);
% Position of the ball at time t, with initial velocity v, 
% and angle theta.
global g
x = v*cos(theta*pi/180)*t;
y = v*sin(theta*pi/180)*t + .5*g*t.*t;
end


function s = distancefromrim(t,v0,theta)
% Calculate the distance from the front of the rim for the ball
% at time t, with initial velocity v, and angle theta. Use this 
% distance to calculate how much space is left between the ball
% and the rim. If this function is negative, the ball hits the 
% rim. If it is zero, the ball scims the rim for this velocity
% and this release angle at time t.

global Dr Db l h 
%global v pos

[x,y] = position(v0,theta,t);
% Calculate the position of the center of the ball at time t,
% with initial velocity v and release angle theta.

part1 = x -(l - Dr/2);
part2 = y - h;

% Calculate the distance from the front of the rim.
s = sqrt(part1.*part1 + part2.*part2)-Db/2;
end


function dist = distancefromback(v,theta)
% Calculate the distance from the back of the rim for the ball
% when its center is level with the rim. The ball is assumed to
% have initial velocity v, and initial angle theta.
global Dr Db l h g 

sintheta = sin(theta*pi/180);
sqrtval  = v*v*sintheta*sintheta+2*g*h;

if (sqrtval < 0)
    dist = -inf;
    return
end
x = v * cos(theta*pi/180)/(-g);
x = x*(v*sintheta + sqrt(sqrtval));
dist = x-l+(Db-Dr)/2;
end


function y = frontfzero(x, v)
    global Dr Db l h g

    sintheta = sin(x * pi / 180);
    % Calculate the relevant time interval.
    timeinterval = [(l - Dr / 2) / (v * cos(x * pi / 180)), -1 / g * (v * sintheta + sqrt(v * v * sintheta * sintheta + 2 * g * h))];
    % If this time interval is not a valid interval, return a negative value.
    % This signals that the ball either goes through the front rim or never reaches it.
    if (timeinterval(1) > timeinterval(2))
        y = -1;
        return;
    end

    % Discretize the time interval
    num_steps = 100; % Number of time steps
    dt = diff(timeinterval) / num_steps; % Time step size

    % Initialize variables for position and velocity
    y = inf; % Initialize minimum distance to infinity
    vy = v * sintheta;
    y_temp = 0; % Temporary variable to store position

    % Perform Euler's method to calculate the position at each time step
    for i = 1:num_steps
        % Update position and velocity
        y_temp = y_temp + vy * dt;
        vy = vy - g * dt;

        % Check if the ball has hit the front rim
        if y_temp <= Db
            y = y_temp;
            return;
        else
            y = min(y, y_temp); % Update the minimum distance
        end
    end
end


function error = calcfronterror(theta0, v0);
% Routine to calculate the maximum allowed error in the angle
% with respect to the front of the rim when the ball is thrown
% with given velocity.

% Calculate the angle when the ball just skims the front rim, 
% and still goes in.
ang  = fzero(@(x) frontfzero(x,v0), theta0);
ang2 = fzero(@(x) frontfzero(x,v0), theta0+1);
% Calculate the error allowed.
error = min(abs(theta0-ang),abs(theta0-ang2));
end


function [error, ang] = calcbackerror(theta0, v0);
% Routine to calculate the maximum allowed error in the angle
% with respect to the back of the rim when the ball is thrown
% with given velocity.
% Joerg Gablonsky, 06/07/2005
global Dr Db l h g 
% Find the maximum distance from the back the ball thrown 
% with this velocity, and varying angles can have. Note that 
% we have to multiply the distance with (-1) to use the Matlab
% fminsearch function to find a maximum.
[ang, val] = fminsearch(@(x) -distancefromback(v0,x), theta0);

% Negate the value to reverse the multiplication with (-1) that
% was necessary to do a maximization.
val = -val
% If this maximum is small enough, the ball cannnot reach the
% back of the rim. Therefore the error can be infinite.
if (val < 0)
    error = NaN;
else
    if (val > Dr-Db/2)
        error = NaN;
    else
        % The ball can reach the back of the rim. Find the 
        % angle when the ball just skims the back of the rim 
        % by minimizing the negative distance from the rim.
        ang = fzero(@(x) -distancefromback(v0,x), theta0);
        % Calculate the error in angle allowed.
        error = abs(theta0-ang);
    end
end
end


function error = calcerror(y)
% Set all variables.
theta0 = y(1);
v0     = y(2);

% Check if the angle is in the valid range, that is, between 
% 10 and 85 degrees.
if ((theta0 < 10) || (theta0 > 85))
    sprintf(...
    'Angle theta0 = %5.2f is either too small or too large.'
    ,theta0)
    error = NaN;  % The ball is too far from the back of the
                  % rim to still be in the rim.
    return
end
    
% Calculate this value since it is used several times below.
sintheta = sin(theta0*pi/180);
global g h l Dr Db
% Calculate the horizontal position of the ball as it comes 
% back down to the basket height.
x = v0 * cos(theta0*pi/180)/(-g);
x = x*(v0*sintheta + sqrt(v0*v0*sintheta*sintheta+2*g*h));

if (x < l - Dr/2+Db/2)
    sprintf('The ball does not reach the basket.')
    error = NaN;  % The ball is too far from the back of the 
                  % rim to still be in the rim.
    return
end
if (x > l + Dr/2 - Db/2)
    sprintf('The ball goes too far.')
    error = NaN;  % The ball is too close to the back of the 
                  % rim or behind the back of the rim.
    return
end

% Check to see if the ball hits the front of the rim. This is
% done by calculating the minimum distance the ball has from 
% the front of the rim. 
% If this distance is below 0, the ball hits the front rim.
interm = frontfzero(theta0, v0);
if (interm < 0)  % The ball hits the front rim.
    sprintf('The ball hits the front of the rim.')
    error = NaN;
    return
end
fronterror = calcfronterror(theta0, v0);
backerror  = calcbackerror(theta0, v0);
error = -min(fronterror, backerror);
end


function error = calcerrorvopt(theta0)
    % First calculate the optimal velocity for this angle.
    v0 = calcOptv(theta0);
    % Use the general routine to calculate the maximum allowable 
    % error given an initial velocity and angle.
    error = -calcerror([theta0, v0]);
end 

thetas = [40:0.4:48, 48.15:0.001:48.2, 48.6:0.4:56];
errors = thetas;
for i = 1:size(thetas, 2)
    errors(i) = calcerrorvopt(thetas(i));
end


plot(thetas, errors, '.-')
grid on
ylabel('e(theta0)');
xlabel('theta0');

[maxError, maxIndex] = max(errors);
bestTheta_deg = thetas(maxIndex);

fprintf('Best shooting angle: %.2f degrees\n', bestTheta_deg);
fprintf('Maximum error: %.2f\n', maxError);