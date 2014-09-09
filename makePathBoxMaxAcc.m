function [path, pathObjSubNum] = makePathBoxMaxAcc(startPoint,endPoint,scanVel,nLines,orientation,maxAcc,dt)
% makes a scan box
% new version takes an additional parameter, the scanVelocity ... this is
% in "natural" units, (distance / dt), so is equivalent to ds per time
% point

% orientation is as follows:
%   1 - path comes in/out, oriented left-right (default)
%   2 - path comes in/out, oriented up-down

if nargin < 4
    error 'too few parameters in makePathBox'
end

if nargin < 5
    orientation = 1;     % does not need to be set by the user
end

if orientation == 1
    % horizontal lines
    spacing = linspace(startPoint(2),endPoint(2),nLines);  % parallel spacings of lines

    firstLine = makePathLineVel(startPoint,[endPoint(1) spacing(1)],scanVel);

    path = firstLine;
    pathObjSubNum = ones(size(firstLine,1),1);    % pathObjSubNum is labels the line number, or 0 for a turn

    for i = 2:length(spacing)
        if mod(i,2)
            nextline = makePathLineVel([startPoint(1) spacing(i)],[endPoint(1) spacing(i)],scanVel);
        else
            nextline = makePathLineVel([endPoint(1) spacing(i)],[startPoint(1) spacing(i)],scanVel);
        end
    
        % calculate the spline to connect the previous path with the new line
        pathTurn = splineFuncMaxAcc(path,nextline,maxAcc,dt);    % note - n turn points is hard-coded here
        path = [path; pathTurn; nextline];
        pathObjSubNum = [ pathObjSubNum; zeros(size(pathTurn,1),1); i*ones(size(nextline,1),1) ];
    end
else
    % vertical lines
    spacing = linspace(startPoint(1),endPoint(1),nLines);  % parallel spacings of lines

    firstLine = makePathLineVel(startPoint,[spacing(1),endPoint(2)],scanVel);

    nPointsTurn = round(length(firstLine)/3); %use a third of a line length to turn, maybe a better solution

    path = firstLine;
    pathObjSubNum = ones(size(firstLine,1),1);    % pathObjSubNum is labels the line number, or 0 for a turn

    for i = 2:length(spacing)
        if mod(i,2)
            nextline = makePathLineVel([spacing(i) startPoint(2)],[spacing(i) endPoint(2)],scanVel);
        else
            nextline = makePathLineVel([spacing(i) endPoint(2)],[spacing(i) startPoint(2)],scanVel);
        end
    
        % calculate the spline to connect the previous path with the new line
        pathTurn = splineFunc(path,nextline,maxAcc,dt);  
        path = [path; pathTurn; nextline];
        pathObjSubNum = [ pathObjSubNum; zeros(size(pathTurn,1),1); i*ones(size(nextline,1),1) ];
    end
    
end
    
