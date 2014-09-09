function path = splineFuncMaxAcc(line1,line2,maxAcc,dt)
% version of the spline finding function which calculates the number of points in the turn
% based on the maxAcc that the program is passed

%
% Part I - calculate the length (in time => samples) of the turn,
% so that the spline (mirrors) run just to the maxAcc threshold, and do not
% accelerate faster
%
% For this portion of the code, it is easier to deal with time explicitly, rather than
% samples.  The time for the turn is converted into a number of points for the second
% part of the function, whith uses units where the spline starts at t=0 and end at t=1
%
%   The governing equations are:
%
%     y = a + b*t + c*t^2 + d*t^3       or
%     y = xi + vi*t +  ai*t*2 + j*t^3
%
%   The acc is a linear function of time, and the maxAcc will occur either at 
%   t=0 or t=tau (at the end of the path).
%
%   The value of tau which will give maxAcc (m) at t=tau is given by the solution to:
%
%     0 = -m*tau^2 + (4*vf + 2*vi)*tau + (6*xi - 6*xf)
%
%   note in terms of finding roots to the quadratic equation,
%   this polynomial is represented in matlab (as a function of tau) as
%  
%   [-m; (4*vf + 2*vi); (6*xi - 6*xf)]

if nargin < 2
    error 'error, too few parameters passed to splineFunc'
end


xi = line1(end,:);
vi = (line1(end,:)-line1(end-1,:)) / dt;

xf = line2(1,:);
vf = (line2(2,:)-line2(1,:)) / dt;



% calculate the end time, at which the acceleration would be maxAcc
% since there are two scan mirrors, need to look at two seperate sets of roots ...
% note of the positive and negative acceleation values, roots could be imaginary,
% for instance, if the acctual acceleration for the path is always positive, there will be
% no point at which the negative maxAcc is reached

%%%% for the case where the max acc occurs at the beginning:
% what values of tau result in a positive maxAcc, at the beginning of the path?
tauPosBegAxis1 = roots( [maxAcc, 4*vi(1) + 2*vf(1), 6*xi(1) - 6*xf(1)] );   % x portion  note vi and vf switched
tauPosBegAxis2 = roots( [maxAcc, 4*vi(2) + 2*vf(2), 6*xi(2) - 6*xf(2)] );   % y portion note vi and vf switched

% what values of tau result in a negative maxAcc, at the beginning of the path?
tauNegBegAxis1 = roots( [-maxAcc, 4*vi(1) + 2*vf(1), 6*xi(1) - 6*xf(1)] );  % x portion note vi and vf switched
tauNegBegAxis2 = roots( [-maxAcc, 4*vi(2) + 2*vf(2), 6*xi(2) - 6*xf(2)] );  % y portion note vi and vf switched

%%%% for the case where the max acc occurs at the end:
% note, the only difference is that vi and and vf are swapped!

% what values of tau result in a positive maxAcc, at the end of the path?
tauPosEndAxis1 = roots( [-maxAcc, 4*vf(1) + 2*vi(1), 6*xi(1) - 6*xf(1)] );  % x portion
tauPosEndAxis2 = roots( [-maxAcc, 4*vf(2) + 2*vi(2), 6*xi(2) - 6*xf(2)] );  % y portion

% what values of tau result in a negative maxAcc, at the end of the path?
tauNegEndAxis1 = roots( [maxAcc, 4*vf(1) + 2*vi(1), 6*xi(1) - 6*xf(1)] );  % x portion
tauNegEndAxis2 = roots( [maxAcc, 4*vf(2) + 2*vi(2), 6*xi(2) - 6*xf(2)] );  % y portion

%%%% figure out which tau to actually use 
% collect together possible tau values for each axis, discarding imaginary value

% for the beginning
tausBegAxis1 = [isreal(tauPosBegAxis1(1))*tauPosBegAxis1(1), ...
                isreal(tauPosBegAxis1(2))*tauPosBegAxis1(2), ...
                isreal(tauNegBegAxis1(1))*tauNegBegAxis1(1), ...
                isreal(tauNegBegAxis1(2))*tauNegBegAxis1(2) ];
tausBegAxis1 = tausBegAxis1( tausBegAxis1 > 0);                  % keep only positive values

tausBegAxis2 = [isreal(tauPosBegAxis2(1))*tauPosBegAxis2(1), ...
                isreal(tauPosBegAxis2(2))*tauPosBegAxis2(2), ...
                isreal(tauNegBegAxis2(1))*tauNegBegAxis2(1), ...
                isreal(tauNegBegAxis2(2))*tauNegBegAxis2(2) ];
tausBegAxis2 = tausBegAxis2( tausBegAxis2 > 0);                  % keep only positive values

% for the end
tausEndAxis1 = [isreal(tauPosEndAxis1(1))*tauPosEndAxis1(1), ...
                isreal(tauPosEndAxis1(2))*tauPosEndAxis1(2), ...
                isreal(tauNegEndAxis1(1))*tauNegEndAxis1(1), ...
                isreal(tauNegEndAxis1(2))*tauNegEndAxis1(2) ];
tausEndAxis1 = tausEndAxis1( tausEndAxis1 > 0);                  % keep only positive values

tausEndAxis2 = [isreal(tauPosEndAxis2(1))*tauPosEndAxis2(1), ...
                isreal(tauPosEndAxis2(2))*tauPosEndAxis2(2), ...
                isreal(tauNegEndAxis2(1))*tauNegEndAxis2(1), ...
                isreal(tauNegEndAxis2(2))*tauNegEndAxis2(2) ];
tausEndAxis2 = tausEndAxis2( tausEndAxis2 > 0);                  % keep only positive values

tauArray = [tausBegAxis1 tausEndAxis1 tausBegAxis2 tausEndAxis2]; % array of possible taus


% now,  need to pick the shortest tau with acc < 1.01 maxAcc
% extrema in accelerations happen at the end, so just check the end points

%tauArray = sort(tauArray);    % short from shortest to longest turn time

acceptableTau = [];
for ti = 1:length(tauArray)
    tau = tauArray(ti);

    accInit = (6*xf - 6*xi)*tau^-2 - (4*vi + 2*vf)*tau^-1;      % initial acceleration
    accFinal = (6*xf - 6*xi)*tau^-2 - (4*vf + 2*vi)*tau^-1;     % final acceleration
    
    if max(abs([accInit accFinal])) < maxAcc * 1.01             % can be 1% higher
        acceptableTau = [acceptableTau tau];
    end
end

%disp(['num possible choices: ' num2str(length(acceptableTau))])

if length(acceptableTau) == 0
    % this should never happen ...
    beep
    disp 'acceptable Tau not found, using longest tau'
    tau = max(tauArray);
else
    tau = min(acceptableTau);    % choose the minimum tau
end

%%%% see if a smaller value of tau would work as well ...
% tauStep = tau / 1000;
% for tauShorter = tau:-tauStep:tauStep
%     accInit = (6*xf - 6*xi)*tau^-2 - (4*vi + 2*vf)*tauShorter^-1;      % initial acceleration
%     accFinal = (6*xf - 6*xi)*tau^-2 - (4*vf + 2*vi)*tauShorter^-1;     % final acceleration
%     
%     % loop until acceleration is too high
%     if max(abs([accInit accFinal])) > maxAcc
%         break;
%     end
% end

% tauShorter = tauShorter + tauStep;   % last point was too fast, go back one
% %disp(['tau ' num2str(tau) ' tau shorter ' num2str(tauShorter)])
% tau = tauShorter;                                                  % set to shortest possible length

%disp( mat2str( max( abs( [accInit accFinal]))))

% convert tau to nPoints
nPoints = ceil(tau / dt);





%There is a quirk in MpScope3 that screws up the display (but not the data)
%if the number of pixels per line is not a multple of 4.
% if mod(nPoints,4) == 1,
    nPoints = ceil(nPoints/4) * 4  + 1;%force number of points to be multiple of 4 after trimming(to circumvent a quirk in MpScope3 code for display) %pst - added 2/3/14
% end%if mod(nSteps,2) == 1,





%%%% Part II - calculate the spline path, using the nPoints calculated above

% calculates a spline path through the passed in points, in any number of dimensions
% generate a cubic spline to connect two points alpha and beta
%
% spline has the equation:
%   P = a + b*t + c*t^2 + d*t^3
%
% and the derivative equation 
%   P' = b + 2*c*t + 3*d*t^2
%
% the constranits are that the positions and the derivatives must match
% at the end of the spline points ... the algebra becomes much simpler if
% the start and the end of the path are taken to be t = 0 and t = 1
% 
% this version of the program takes two lines, where the points are listed
% as columns, and creates the spline between them
%
% note, that the first and last points are cut out of the spline ...
% these points are the same as the ends of the lines that they are
% connecting

% parse endpoints and slopes from input lines
point1 = line1(end,:);
point1slope = line1(end,:)-line1(end-1,:);

point2 = line2(1,:);
point2slope = line2(2,:)-line2(1,:);

% calculate spline
alpha = point1;
beta = point2;

% the velocities must be scaled by the step size, so that they are in the
% same units
alphaPrime = nPoints * point1slope;
betaPrime = nPoints * point2slope;

a = alpha;
b = alphaPrime;
c = - 3*alpha - 2*alphaPrime + 3*beta - betaPrime;
d = 2*alpha + alphaPrime - 2*beta + betaPrime;

%Plength = b + c + d

tStep = 0:1/nPoints:1 ;

path = zeros(length(tStep),length(alpha));  % vector to hold the path 

for tIter= 1:length(tStep)
    t = tStep(tIter);
    path(tIter,:) = a + b*t + c*t^2 + d*t^3;
end

% cut out two end points, since these are already part of the line
path = path(2:end-1,:);