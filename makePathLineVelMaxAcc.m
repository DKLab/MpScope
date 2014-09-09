function [scanPath, scanPathObjSubNum] = makePathLineVelMaxAcc(startPoint,endPoint,scanVel,nLines,maxAcc,dt)

% function takes start and stop coordinates, and creates a line scan path
% if there is more than one line, they are joined by splines and
% scanned back and fourth
% version  "makePathLineVel" takes a scan velocity, which is in "natural units", 
% that is, distance per time step.  For instance, a velocity of 0.1 
% between 0 and 1 will generate points [0 0.1 0.2 0.3], etc.

%% parameters
ds = scanVel;                       % in natural units, these are the same
length = norm(startPoint-endPoint);

%% parse inputs

if nargin < 5
    nStepsTurn = 10;              % need a better number for this ...
end
if nargin < 4
    nLines = 1;                   % one line is the default
end
if nargin < 3
    error 'too few parameters passed function makePathLineVel'
end

%% find path for a line
% this works by making two lines, the odd numbered lines (which go throught
% the orignal points in order), and the even number lines (which go through
% the original points in reverse.  then, these points are connected with
% splines

% generate the line
nSteps = round(length / ds);      % how many steps to take along the length of the vessel


%There is a quirk in MpScope3 that screws up the display (but not the data)
%if the number of pixels per line is not a multple of 4.

% if mod(nSteps,4) == 1,
    nSteps = ceil(nSteps/4) * 4; %force number of points to be multiple of 4 (to circumvent a quirk in MpScope3 code for display) %pst - added 2/3/14
% end%if mod(nSteps,2) == 1,

    
    
% initial scan line
scanPathOdd = [];                    % will hold x and y scan info
for n = 0:nSteps-1
    point = startPoint + n/(nSteps-1) * (endPoint-startPoint);
    scanPathOdd = [scanPathOdd; point];
end

% stop if only one line is to be generated
if nLines == 1
    scanPath = scanPathOdd;
    scanPathObjSubNum = ones(size(scanPath,1),1);
    return
end

scanPathObjSubNum = ones(size(scanPathOdd,1),1);

%% generate more than one line, and connect them by splines
% generate the even (return) lines, and the splines to connect to two

scanPathEven = makePathLineVel(endPoint,startPoint,scanVel);  % flipped line

% create the spline functions to connect the two line types
% (connect even odd lines with even lines)
scanPathSplineOddToEven = splineFuncMaxAcc(scanPathOdd,scanPathEven,maxAcc,dt);

% (connect even lines with odd lines)                 
scanPathSplineEvenToOdd = splineFuncMaxAcc(scanPathEven,scanPathOdd,maxAcc,dt);

scanPath = scanPathOdd;    % start with the first line                      
for li = 2:nLines          % iterate to add the rest
    % add the turn, and line
    if mod(li,2)            % li is odd
        scanPath = [scanPath; scanPathSplineEvenToOdd; scanPathOdd];
        scanPathObjSubNum = [scanPathObjSubNum; zeros(size(scanPathSplineEvenToOdd,1),1); li*ones(size(scanPathOdd,1),1)];
    else
        scanPath = [scanPath; scanPathSplineOddToEven; scanPathEven];
        scanPathObjSubNum = [scanPathObjSubNum; zeros(size(scanPathSplineOddToEven,1),1); li*ones(size(scanPathEven,1),1)];
    end
end
                      



