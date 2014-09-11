function varargout = MpArbScanGUI(varargin)
% GUI to define, draw, and scan paths
% new version, for MpScope 3 (started 2012 June)
% Jonathan Driscoll & Philbert Tsai, David Kleinfeld Lab of Neurophysics, UCSD

% change log:
%01/28/14 - Changed "scanVelocity" variable to "scanStepSize" to be more accurate [volts (or millivolts) per pixel].
% Last Modified by GUIDE v2.5 09-Sep-2014 10:21:03

%% INITIALIZATION - MATLAB and user initialzation code, and exit
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MpArbScanGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MpArbScanGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before MpArbScanGUI is made visible.
function MpArbScanGUI_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.output = hObject;   % Choose default command line output for MpArbScanGUI
    
    % user initialization code here ...
    
    % if MpArbScanGUI is called with the line 'Save', open autosave function and exit
    if numel(varargin) && ischar(varargin{1}) && strcmpi(varargin{1},'save')
        disp 'Save function called on MpArbScanGUI ...'
        %handles
        figure(handles.pathGUIfig)     % make sure window pops up, instead of opening in background
        pushButtonSaveData_Callback(hObject, eventdata, handles)
        return   % do not exectute the rest of the initializeation code, which would overwrite handles
    end
    
    % create initial handles fields
    handles.path = [];            % initialize a blank path
    handles.returnedPath = [];    % initialize a blank path (from ADC)
    handles.lineIter = 1;         % iterator to keep track of lines
    handles.boxIter = 1;          % iterator to keep track of boxes
    
    handles.pathObjNum = [];      % for each point in the path, which item does it correspond to? (0 for turn)
    handles.pathObjSubNum = [];   % which sub item? (i.e., which line of a multi-line box?), 0 for a turn
    
    handles.im = [];              % blank image 
    handles.axisLimRow = [-2,2];  % image axis limits
    handles.axisLimCol = [-2,2];  % image axis limits
    
    handles.scanStepSize = 1e-3 * str2num(get(handles.editScanStepSize,'String'));  %entered in mV/pixel, convert to V/pixel

    handles.fileDirectory = '.\';     % initial file directory 
    
    handles.drawingLine = false;      % used by draw line callback
    handles.lineHandle = 0;           % handle to line being drawn
    
    handles.pathGUIfig = gcf;     % store the current figure
    
    % NEW: the struct array handles.group will contain a scanCoords struct
    % for every group that the user creates. A group is a set of lines or
    % boxes that that can be repeated independently of other groups.
    handles = createGroup(handles, 1);
    % createGroup will create the fields:
        % handles.group   and
        % handles.currentGroupIndex
    
    
    guidata(hObject, handles);    % Update handles structure
    
    pushButtonClearPath_Callback(hObject, eventdata, handles)  % set up blank path and image
    % prepare the scan coords listbox
    set(handles.listboxScanCoords, 'String', {' '});
    set(handles.listboxScanCoords, 'Value', 1);
    
    % no path, so assign the length zero (so Delpi knows), may not be necessary
    assignin('base','SCAN_path_len',int32(0) );
    
    set(gcf,'name','MpArbScanGUI v0.10');   
    
    % uiwait(handles.figure1);   % UIWAIT makes MpArbScanGUI wait for user response (see UIRESUME)
    
% --- Outputs from this function are returned to the command line.
function varargout = MpArbScanGUI_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;  % Get default command line output from handles structure


%% stand-alone functions 
    
function [ groupIndex, scanCoordsIndex, listboxIndex ] = getListboxItem( handles )
    % retrieve the currently selected item in the scan coords listbox.
    % If the item can't be found in handles.group, then 0 is returned for
    % the groupIndex and scanCoordsIndex.
    % If the item is a group and not a line or a box, then groupIndex will
    % have the correct number and scanCoordsIndex will be 0.
    
    listboxIndex = get(handles.listboxScanCoords,'Value');
    listboxStringList = get(handles.listboxScanCoords, 'String');
    
    itemName = listboxStringList{listboxIndex};

    % remove any HTML from the item name
    itemName = regexprep(itemName, '<[^>]*>', '');
    % remove any &nbsp; (HTML non blanking space)
    itemName = strrep(itemName, '&nbsp;', '');
    
    % check if the user clicked on a group
    if ~isempty( strfind(itemName, 'Group') )
        % this is a group, not a line or box
        % get the group number -- one or more digts followed by an open paranthesis '('
        groupIndex = str2double( regexp(itemName, '\d+(?=\()', 'match') );
        scanCoordsIndex = 0;
        
        if isnan(groupIndex)
            groupIndex = 0;
        end
 
    else
        % this is a line or a box, find the two indicies
        [ groupIndex, scanCoordsIndex ] = findIndex( handles, itemName );
    end

function [ groupIndex, scanCoordsIndex ] = findIndex( handles, objectName )
    % look through each scanCoords struct array within the handle.group
    % struct array and return the index that the item was found at, or 0 if
    % the item doesn't exist.
    
    groupIndex = 0;
    scanCoordsIndex = 0;
    
    for searchGroupIndex = 1 : length( handles.group )
       scanCoords = handles.group(searchGroupIndex).scanCoords;
       
       for searchScanCoordsIndex = 1 : length( scanCoords ) 
            if strcmp( objectName, scanCoords(searchScanCoordsIndex).name );
                % item has been found
                groupIndex = searchGroupIndex;
                scanCoordsIndex = searchScanCoordsIndex;
            end
       end
       
    end
    
function handles = moveListboxItem(handles, moveUpFlag)
    % moves an item up or down -- changing its order in the scan coords
    % listbox and changing the order that the path is drawn in
    
    % if moveUpFlag is true, the item will be moved up, if its false it
    % will be moved down
    
     [groupIndex, scIndex, listboxIndex] = getListboxItem(handles);
    
    % if scIndex is 1 then this item can't be moved up any further, and if
    % its 0 then this item is a group
    if moveUpFlag
        inBounds = scIndex > 1;
        increment = -1;
    else
        inBounds = scIndex < length(handles.group(groupIndex).scanCoords);
        increment = 1;
    end
    
    
    if groupIndex > 0 && inBounds;
        
        swapItem = handles.group(groupIndex).scanCoords(scIndex);
        handles.group(groupIndex).scanCoords(scIndex) = ...
                    handles.group(groupIndex).scanCoords(scIndex + increment);
                
        handles.group(groupIndex).scanCoords(scIndex + increment) = swapItem;
        
        % move the selected item as well
        set(handles.listboxScanCoords,'Value', listboxIndex + increment );
          
        % the hObject being sent to updatePath is only used for guidata
        % calls -- so any object in the figure can be used as 'hObject'
        updatePath(handles.listboxScanCoords, [], handles);  % update the scan path
    end        
        
        
function handles = createGroup(handles, optional_GroupIndex)
    if exist('optional_GroupIndex', 'var')
        newIndex = optional_GroupIndex;
    else
        newIndex = length(handles.group) + 1;
    end
    
    % dont ask for number of cycles if this is the first group (just assume
    % 1 cycle)
    if newIndex > 1
    
        prompt = sprintf('How many times should Group %d be repeated before moving on to the next group?',...
                            newIndex);
        dlg_title = 'Group Cycles';
        defAns = {'1'};
        cyclesStr = inputdlg( prompt, dlg_title, 1, defAns);

        if isempty(cyclesStr)
            % user cancled dialog box -- don't create a new group afterall
            return;
        end
        
        
        cycles = str2double(cyclesStr{1});

        if isnan(cycles)
            cycles = 1;
        end
    else
        cycles = 1;
    end
    

    newScanCoords = struct('scanShape', 'blank', ...
                                'startPoint', [], ...
                                'endPoint', [], ...
                                'nLines', [], ...
                                'orientation', [], ...
                                'name','<blank>');

    
    handles.group(newIndex) = struct(...
        'scanCoords', newScanCoords,...
        'cycles', cycles );
    
    handles.currentGroupIndex = newIndex;
    
    refreshListbox(handles);
    
    % select the last item in the listbox -- 
    % this will be where the new group is
    stringList = get( handles.listboxScanCoords, 'String' );
    set( handles.listboxScanCoords, 'Value', length(stringList) );
    

function handles = deleteGroup(handles, groupIndex)

    % user can delete any group except group 1
    if groupIndex > 1
        % set the scan coords listbox to the first element (this is
        % guaranteed to still exist)
        set(handles.listboxScanCoords, 'Value', 1);
        
        handles.group(groupIndex) = [];
        
        % and reset the current group index to the last group
        handles.currentGroupIndex = length(handles.group);
    elseif groupIndex == 1
        disp('Cannot delete Group 1');
    end
    
function mouseDown_Callback(hObject, ~)
    % Any number of graphics objects may have this callback registered as
    % the ButtonDownFcn.
    % On mouse down, check to see if an imline needs to be completed by
    % asking the user if they want to add a new line, finish the path, or
    % create a new group to add lines to.
    
    handles = guidata(hObject);
    groupIndex = handles.currentGroupIndex;
    
    if handles.drawingLine

        qstring = sprintf('Would you like to add a new line to Group %d, start a new Group, or finish placing lines?', ...
                           groupIndex);
        title = 'Line Placed';
        newLine = 'New Line';
        newGroup = 'New Group';
        finishLine = 'Done Placing Lines';
        default = newLine;
        
        answer = questdlg( qstring, title, newGroup, newLine, finishLine, default); 

        switch answer
            case newLine
                handles = StopLine(hObject, [], handles);
                handles = StartNewLine(hObject, [], handles);
                
            case newGroup
                handles = StopLine(hObject, [], handles);
                handles = createGroup(handles);
                handles = StartNewLine(hObject, [], handles);
                
            case finishLine
                handles = StopLine(hObject, [], handles);
                
            otherwise
                % disregard this mouse click
                return;
        end
    else
        % Nothing is being drawn -- ask the user if they want to create a
        % new group, or add a line or a box to the current group
        qstring = sprintf('Would you like to create a new Group, add a line to Group %d, or add a box to Group %d?', ...
                            groupIndex, groupIndex);
        title = 'Add New...';
        newLine = 'Add Line';
        newGroup = 'Create New Group';
        newBox = 'Add Box';
        default = newGroup;
        
        answer = questdlg( qstring, title, newGroup, newLine, newBox, default);
        
        switch answer
            case newLine
                handles = StartNewLine(hObject, [], handles);
            case newGroup
                handles = createGroup(handles);
                guidata(hObject, handles);
            case newBox
                handles = addBox(handles);
                guidata(hObject, handles);
                updatePath(hObject, [], handles);
            otherwise
                return;
        end
    end
    
    
    
function refreshListbox(handles)
    % make a cell structure of existing names
    for groupIndex = 1 : length(handles.group)
        if groupIndex == 1
            % handle group 1 seperately -- this is always the first item in
            % the list
            strmat = sprintf('<HTML><b>Group 1</b>&nbsp;&nbsp;(x%d)</HTML>',...
                        handles.group(groupIndex).cycles);
        else
            strmat = char(strmat, sprintf('<HTML><b>Group %d</b>&nbsp;&nbsp;(x%d)</HTML>', ...
                        groupIndex, handles.group(groupIndex).cycles));
        end

        for scIndex = 1:length(handles.group(groupIndex).scanCoords)
            if ~strcmp( handles.group(groupIndex).scanCoords(scIndex).scanShape, 'blank') 
                strmat = char( strmat, ...
                    sprintf('<HTML>&nbsp;&nbsp;&nbsp;&nbsp;%s</HTML>',...
                        handles.group(groupIndex).scanCoords(scIndex).name));
            end
        end
    end
    set(handles.listboxScanCoords,'String',cellstr(strmat));
    

function handles = updatePath(hObject, eventdata, handles,extra)
    % funtion is called whenever a new point is added to the scanCoords
    
 
    fast = (nargin == 4) && strcmp(extra, 'fast');

    %fast = false; %jd - don't use this for now

    assignin('base', 'scanCoords', handles.group(1).scanCoords);
    
    % create a list of scanCoords that includes the scanCoords from every
    % group
    nGroups = length(handles.group);
    scanCoords = [];        % the struct array that will hold the full scan path
                    % with all group cycles included
                    
    for groupIndex = 1 : nGroups
        cycles = handles.group(groupIndex).cycles;
        
        newScanCoords = handles.group(groupIndex).scanCoords;
 
        % handle repetitions (the first cycle is already accounted for)
        while cycles > 1
            cycles = cycles - 1;
            newScanCoords = [ newScanCoords, ...
                                    handles.group(groupIndex).scanCoords ];      
        end
        
        if isempty(scanCoords)
            scanCoords = newScanCoords;
        else
            % (preallocating scanCoords may be tricky -- would need a seperate loop
            % before this one to calculate the total length that sc will be)
            scanCoords = [ scanCoords, newScanCoords];
        end 
    end
    
    if ~fast
        % plot the start and endpoints on the graph, and place text
        for i = 1:length(scanCoords)
            sc = scanCoords(i);
            
            if ~strcmp(sc.scanShape,'blank')

                % mark start and end point
                plot(sc.startPoint(1),sc.startPoint(2),'g*')
                plot(sc.endPoint(1),sc.endPoint(2),'r*')

                % find a point to place text
                placePoint = sc.startPoint + .1*(sc.endPoint-sc.startPoint);
                text(placePoint(1),placePoint(2),sc.name,'color','red')

                % draw a line or box (depending on data structure type)
                if strcmp(sc.scanShape,'line')
                    line([sc.startPoint(1) sc.endPoint(1)],[sc.startPoint(2) sc.endPoint(2)])
                elseif strcmp(sc.scanShape,'box')
                    % width and height must be > 0 to draw a box
                    boxXmin = min([sc.startPoint(1),sc.endPoint(1)]);
                    boxXmax = max([sc.startPoint(1),sc.endPoint(1)]);
                    boxYmin = min([sc.startPoint(2),sc.endPoint(2)]);
                    boxYmax = max([sc.startPoint(2),sc.endPoint(2)]);

                    rectangle('Position',[boxXmin,boxYmin, ...
                        boxXmax-boxXmin,boxYmax-boxYmin], ...
                        'EdgeColor','green');
                end
            end
        end
    end % fast 
    
    if ~fast
        % update the listbox
        refreshListbox(handles);
    end % fast
        
    % update the actual path
    handles.path = [];          % clear previously set path
     

    if isempty(scanCoords) || strcmp(scanCoords(1).scanShape,'blank')
        guidata(hObject, handles);
        return;   % no path to find
    end
    
    % actual scan path is created here!
    % pathObjNum and pathObjSumNum are also created here ...
    
    % TODO: Once again, will need to iterate over groups and handle
    % transitions between groups
    
    for i = 1:length(scanCoords)
        sc = scanCoords(i);     % copy to a structure, to make it easier to access
        
        if strcmp(sc.scanShape,'line')
            % create a line here
            [pathLine objSubNum] = makePathLineVelMaxAcc(sc.startPoint,sc.endPoint,handles.scanStepSize,sc.nLines,handles.maxAcc,handles.pixelDwellTime);

            % add a spline to connect to previous path, if previous path is non-zero
            if isempty(handles.path)
                handles.path = pathLine;                           % this is the first line, just add
                handles.pathObjNum = ones(size(pathLine,1),1);     % first line
                handles.pathObjSubNum = objSubNum;
            else
                pathSpline = splineFuncMaxAcc(handles.path,pathLine,handles.maxAcc,handles.pixelDwellTime);      % splineFunc determines turn points
                handles.path = [handles.path; pathSpline; pathLine]; 
                handles.pathObjNum = [handles.pathObjNum; zeros(size(pathSpline,1),1); i*ones(size(pathLine,1),1)]; 
                handles.pathObjSubNum = [handles.pathObjSubNum; zeros(size(pathSpline,1),1); objSubNum]; 
            end
            
            
        elseif strcmp(scanCoords(i).scanShape,'box')
            % create a box here
            [pathBox objSubNum] = makePathBoxMaxAcc(sc.startPoint,sc.endPoint,handles.scanStepSize,sc.nLines,1,handles.maxAcc,handles.pixelDwellTime);
 
            % add a spline to connect to previous path, if previous path is non-zero
            if isempty(handles.path)
                % no need need to connect to previous line
                handles.path = pathBox;                
                handles.pathObjNum = i*ones(size(pathBox,1),1);    % first line, object num is 1
                handles.pathObjSubNum = objSubNum;
            else
                % add previous connecting spline, and box
                %pathSpline = splineFunc(handles.path,pathBox);  % path spline determines nPoints
                pathSpline = splineFuncMaxAcc(handles.path,pathBox,handles.maxAcc,handles.pixelDwellTime);  % path spline determines nPoints
                
                handles.path = [handles.path; pathSpline; pathBox];
                handles.pathObjNum = [handles.pathObjNum; zeros(size(pathSpline,1),1); i*ones(size(pathBox,1),1) ]; 
                handles.pathObjSubNum = [handles.pathObjSubNum; zeros(size(pathSpline,1),1); objSubNum];
            end


        else 
            warndlg('scanType set incorrectly')    % undefined scan type
        end
    end % loop over scan objects
    
    % close the path (loop back to beginning)
    %closeSpline = splineFunc(handles.path,handles.path);
    closeSpline = splineFuncMaxAcc(handles.path,handles.path,handles.maxAcc,handles.pixelDwellTime);
    
    

    
    
    handles.path = [handles.path; closeSpline];
    handles.pathObjNum = [handles.pathObjNum; zeros(size(closeSpline,1),1)];
    handles.pathObjSubNum = [handles.pathObjSubNum; zeros(size(closeSpline,1),1)];
    
    % make sure paths were added correctly ...
    if size(handles.path,1) ~= size(handles.pathObjNum,1)
        warndlg 'oops, path lengths are not the same ...'
    end
    
    guidata(hObject, handles);   
    
    % write these variables to the base space, so that Delphi can access them
    scanN = size(handles.path,1);  % number of points in scan

%pst 03/11/14 - removed negative flip on x and added pixel shift offset
% %     negX = -handles.path(:,1);
% %     assignin('base','SCAN_path_x',negX );                % needs to be flipped for mirrors this! %pst 03/11/2014 - removed negative flip!
% %     assignin('base','SCAN_path_x',handles.path(:,1) );
% %     assignin('base','SCAN_path_y',handles.path(:,2) );   

    pixelsToShiftX = int16(handles.pixelsToShift(1));
    pixelsToShiftY = int16(handles.pixelsToShift(1));
    shiftedPathVoltageX = circshift(handles.path(:,1),double([-1*pixelsToShiftX,0]));
    shiftedPathVoltageY = circshift(handles.path(:,2),double([-1*pixelsToShiftY,0])); 
% %     assignin('base','SCAN_path_x',shiftedPathVoltageX );
% %     assignin('base','SCAN_path_y',shiftedPathVoltageY );
% %     assignin('base','SCAN_path_len',int32(scanN) );  %removed pixel shift - taken care of in MpScope3 (processFrame).
    
    assignin('base','SCAN_path_x',handles.path(:,1) );
    assignin('base','SCAN_path_y',handles.path(:,2) );
    assignin('base','SCAN_path_len',int32(scanN) );
    
    
    
    % handle Z data, if it exists
    if size(handles.path,2) == 3
        % Z is assigned
        assignin('base','SCAN_path_z',handles.path(:,2) );
        assignin('base','SCAN_path_z_len',int32(scanN) );
    else
        assignin('base','SCAN_path_z',0 );
        assignin('base','SCAN_path_len_z',int32(0) );
    end

    
    
    
    
function handles = StartNewLine(hObject, ~, handles) 
    useImline = true;          % use imline to draw line, if available
    set(handles.toggleButtonAddLine,'Visible','off');
    set(handles.pushbuttonLineDone,'Visible','on');
    set(handles.pushbuttonNextLine,'Visible','on');
    
    %disp 'starting line'
    % start by getting point from use
    % this is necessary, since imline(gca,[]) doesn't work in a gui -
    % must specify the starting points
    startPoint = ginput(1);
    plot(startPoint(1),startPoint(2),'g*')
    endPoint = ginput(1);
    plot(endPoint(1),endPoint(2),'r*')

    
    %if exist('imline') == 2 && useImline
    if useImline
        % draw a placeable line, using imline
        handles.drawingLine = true; 
        handles.lineHandle = imline(gca,[startPoint(1) endPoint(1)],[startPoint(2) endPoint(2)]);
        guidata(hObject, handles);
        return;
    end


function handles = StopLine(hObject, eventdata, handles) 
    set(handles.toggleButtonAddLine,'Visible','on');
    set(handles.pushbuttonLineDone,'Visible','off');
    set(handles.pushbuttonNextLine,'Visible','off');   
    
        API = iptgetapi(handles.lineHandle);      
        position= API.getPosition();              % get the line's position
        API.delete();                             % delete the line
            
        startPoint = position(1,:);
        endPoint = position(2,:);
            
        %guidata(hObject, handles);     % there's already a call to guidata
                                        % later in this function


    %%%
    
    % store the line
    lineName = ['line ' num2str(handles.lineIter)];
    handles.lineIter = handles.lineIter + 1;
    handles.drawingLine = false;
    
    sc = struct('scanShape', 'line', ...
                'startPoint', startPoint, ...
                'endPoint', endPoint, ...
                'nLines', handles.nLines, ...
                'orientation', 1, ...
                'name', lineName );
    
    groupIndex = handles.currentGroupIndex;
            
    if isempty(handles.group(groupIndex).scanCoords) || ...
            strcmp(handles.group(groupIndex).scanCoords(1).scanShape, 'blank')
        
        handles.group(groupIndex).scanCoords = sc;                             % this is the first element
    else
        handles.group(groupIndex).scanCoords( end + 1 ) = sc;  % append to end
    end
    
    guidata(hObject, handles);   % Update handles structure   
    handles = updatePath(hObject, eventdata, handles);  % update the scan path
    
    
    
function handles = addBox(handles)
    startPoint = ginput(1);
    plot(startPoint(1),startPoint(2),'g*')

    endPoint = ginput(1);
    plot(endPoint(1),endPoint(2),'r*')  
    
    boxName = ['box ' num2str(handles.boxIter)];
    handles.boxIter = handles.boxIter + 1;

    sc = struct('scanShape', 'box', ...
                'startPoint', startPoint, ...
                'endPoint', endPoint, ...
                'nLines', handles.numBoxLines, ...
                'orientation', 1, ...
                'name', boxName);
    
    groupIndex = handles.currentGroupIndex;
    
    if isempty(handles.group(groupIndex).scanCoords) || ...
            strcmp(handles.group(groupIndex).scanCoords(1).scanShape, 'blank')
        
        handles.group(groupIndex).scanCoords = sc;                            % this is the first element
    else
        handles.group(groupIndex).scanCoords( end + 1 ) = sc; % append to end
    end



    
    
    
    
%% BUTTONS    


% --- BUTTON - Add Line
function toggleButtonAddLine_Callback(hObject, eventdata, handles)

    StartNewLine(hObject, eventdata, handles);
    
    
    
% --- BUTTON - Add Box
function pushButtonAddBox_Callback(hObject, eventdata, handles)
    handles = addBox(handles);
    
    guidata(hObject, handles);   % Update handles structure 
    updatePath(hObject, eventdata, handles);
    
    
% --- BUTTON - Draw Path
function pushButtonDrawPath_Callback(hObject, eventdata, handles)
    
    if handles.drawingLine == true,
        handles = StopLine(hObject, eventdata, handles);
    end
    
    RefreshAllParameters(hObject, eventdata, handles);
    pushButtonClearGraph_Callback(hObject, eventdata, handles) % clear the graph
    updatePath(hObject, eventdata, handles);        % make sure path is updated, and draw boxes / lines
    handles = guidata(hObject);

    
    % to draw a path that can be saved, put a breakpoint in the program, and type
    %   figure
    %   imagesc(handles.axisLimCol,handles.axisLimRow,handles.im);colormap('gray');axis image; axis off
    %   hold on;plot(handles.path(:,1),handles.path(:,2),'w.');hold off
    
    %{
    groupIndex = handles.currentGroupIndex;
    if strcmp(handles.group(groupIndex).scanCoords(1).scanShape,'blank')
        % first element is blank (no path to draw!)
        helpdlg('no path to draw ...')
        return
    end
    %}
    
    nPoints = size(handles.path,1);

    comet(handles.path(:,1),handles.path(:,2));
    %return
    
% no need for this, comet is very fast
     color = 'red';    
     plotFast = 1;
     
     dInd=1:handles.drawEveryPoints:nPoints;   % select points to draw
     
     curfig = gcf;    
 
     if plotFast
         plot(handles.path(dInd,1),handles.path(dInd,2),'.','color',color)   
         drawnow
     else
         for d=dInd
             figure(curfig);
             plot(handles.path(d,1),handles.path(d,2),'.','color',color)
             drawnow
         end
     end
                
    %color = hsv2rgb([(1:nPoints)'/nPoints,ones(nPoints,1),ones(nPoints,1)]);  
 
    updatePath(hObject, eventdata, handles);  % draw the coordinates oven the path
    handles = guidata(hObject);
    pathLength_pixels = size(handles.path,1);
    pathDuration_ms = pathLength_pixels  * handles.pixelDwellTime * 1000;
    set(handles.editPathLength,'String',num2str(pathLength_pixels,'%0.0f'));
    set(handles.editPathDuration,'String',num2str(pathDuration_ms,'%0.3f'));
    
    


    
% --- BUTTON - Clear Path
function pushButtonClearPath_Callback(hObject, eventdata, handles)
    % sets the scanCoords back to <blank> value, and drawn image
%     if handles.drawingLine == true,
%         handles = StopLine(hObject, eventdata, handles);
%     end
    
    % delete everything except the first group
    handles.currentGroupIndex = 1;
    handles.group( 2 : end ) = [];
    

    handles.group(1).scanCoords = struct('scanShape', 'blank', ...
                                    'startPoint', [], ...
                                    'endPoint', [], ...
                                    'nLines', [], ...
                                    'orientation', [], ...
                                    'name','<blank>');
                                     % will become the array of the scanstructure,
                                     % note the first element is 'blank'
  
    set(handles.listboxScanCoords, 'Value', 1);
    guidata(hObject, handles);   % Update handles structure
    updatePath(hObject, eventdata, handles);
    
    pushButtonClearGraph_Callback(hObject, eventdata, handles);


% --- BUTTON - Clear Graph
function pushButtonClearGraph_Callback(hObject, eventdata, handles)
    % draws an empty graph, or graph with image if image is set
%     if handles.drawingLine == true,
%         handles = StopLine(hObject, eventdata, handles);
%     end
    
    if isempty( handles.im )
        % plot a blank image
        hold off
        %image(handles.axisLimCol,handles.axisLimRow,100*ones(100,100));
        image([-2 2],[-2 2],100*ones(100,100));
        colormap('gray')
        hold on
    else
        % plot the actual graph
        cla
        testImage = log(max(handles.im,1));
        imageHandle = imagesc(handles.axisLimCol,handles.axisLimRow,handles.im);
%        imagesc(handles.axisLimCol,handles.axisLimRow,testImage);
       
        set(imageHandle, 'ButtonDownFcn', @mouseDown_Callback);

        % just in case the user clicks somewhere other than the image:
        set(gca, 'ButtonDownFcn', @mouseDown_Callback);
        
        axis on
        axis tight
        axis equal
        colormap('gray')
    end
    
    updatePath(hObject, eventdata, handles);


% --- BUTTON - Save Data
function pushButtonSaveData_Callback(hObject, eventdata, handles)
    % load data from handles to a scanData structure, and save it
    % handles
%     if handles.drawingLine == true,
%         handles = StopLine(hObject, eventdata, handles);
%     end
    
    scanData.pixelDwellTime = handles.pixelDwellTime;
%     scanData.fs = handles.fs;
    scanData.scanStepSize = handles.scanStepSize;   %NeedToDo:  changed from ScanVelocity - check downstream.
    scanData.path = handles.path;
    scanData.returnedPath = handles.returnedPath; 
    
    scanData.maxAcc = handles.maxAcc;
    scanData.im = handles.im;
    scanData.axisLimRow = handles.axisLimRow;      % limits in the up-down direction (the rows)
    scanData.axisLimCol = handles.axisLimCol;      % limits in the left right direction (the cols)
    
    %scanData.scanCoords = handles.scanCoords;  % this line is being
    %replaced by handles.group
    scanData.group = handles.group;
    
    %jd - this doesn't make sense, meant to save to scanData.?
    handles.imPath = handles.fileDirectory;                  %  path of loaded data
    
    scanData.pathObjNum =  handles.pathObjNum;      % for each point in the path, which item does it correspond to? (0 for turn)
    scanData.pathObjSubNum = handles.pathObjSubNum; 
         
    [fileName,handles.fileDirectory] = uiputfile(handles.fileDirectory,'save file');
    
    if isequal(fileName,0)              % check to make sure a file was selected
        return
    end
    
    %saveFilename = 'data/savedDate'
    %save([pathName fileName ],'scanData');
    save([handles.fileDirectory fileName ],'scanData');
    %disp([pathName fileName]);
    
    updatePath(hObject, eventdata, handles);

    
    
% --- Executes on button press in pushbuttonLoadData.
function pushbuttonLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [fileName,handles.fileDirectory] = uigetfile(handles.fileDirectory,'load file');
    load([handles.fileDirectory fileName ],'scanData');
    handles.pixelDwellTime = scanData.pixelDwellTime;
    handles.scanStepSize = scanData.scanStepSize;   %NeedToDo:  changed from ScanVelocity - check downstream.
    handles.path = scanData.path;
    handles.returnedPath = scanData.returnedPath; 
    
    handles.maxAcc = scanData.maxAcc;
    handles.im = scanData.im;
    handles.axisLimRow = scanData.axisLimRow;      % limits in the up-down direction (the rows)
    handles.axisLimCol = scanData.axisLimCol;      % limits in the left right direction (the cols)
    
    %handles.scanCoords = scanData.scanCoords;  % just like in the save
    %function, replacing this line with .group
    handles.group = scanData.group;
    
    %jd - this doesn't make sense, meant to save to scanData.?
    handles.imPath = handles.fileDirectory;                  %  path of loaded data
    
    handles.pathObjNum =  scanData.pathObjNum;      % for each point in the path, which item does it correspond to? (0 for turn)
    handles.pathObjSubNum = scanData.pathObjSubNum; 
         
    
    
    if isequal(fileName,0)              % check to make sure a file was selected
        return
    end
    

    
    updatePath(hObject, eventdata, handles);    
    
    
    
% --- Button - Display Path Data
function pushButtonPathAnalysis_Callback(hObject, eventdata, handles)
%     if handles.drawingLine == true,
%         handles = StopLine(hObject, eventdata, handles);
%     end
    pushButtonDisplayPathData_Callback(hObject, eventdata, handles)

    
% --- Button - Display Path Data
function pushButtonDisplayPathData_Callback(hObject, eventdata, handles)
    % calculates information about the scan path, and displays to screen
%     if handles.drawingLine == true,
%         handles = StopLine(hObject, eventdata, handles);
%     end
    path = handles.path;
    
    if isempty(path)
        disp 'path has not been set ...'
        return
    end
    
    pathVelocity = diff(path) ./ handles.pixelDwellTime;
    pathAcceleration = diff(pathVelocity) ./ handles.pixelDwellTime;
    pathPeriod = size(path,1) * handles.pixelDwellTime;
    pathFrequency = 1 / pathPeriod;

    disp(['*** path data'])
    disp(['    position (min/max): ' ...
               num2str(min(handles.path(:))) '  ' num2str(max(handles.path(:)))])
    disp(['    speed kV/s (min/max): ' ...
               num2str(min(abs(pathVelocity(:)))/1e3) '  ' num2str(max(abs(pathVelocity(:)))/1e3)])         
    disp(['    abs accelation MV/s^2 (min/max): ' ...
               num2str(min(abs(pathAcceleration(:)))/1e6) '  ' num2str(max(abs(pathAcceleration(:)))/1e6)])
    disp(['    path Period (s): ' num2str(pathPeriod) ' frequency (Hz): ' num2str(pathFrequency) ])           
    disp(['    nChannels in path: ' num2str(size(handles.path,2))])
    disp(['    nPoints in Path: ' num2str(size(handles.path,1))])
    disp(['    pixelDwellTime for path: ' num2str(handles.pixelDwellTime)]);
    disp(['    percent ROI: ' num2str( 100*sum(handles.pathObjNum>0)/length(handles.pathObjNum)) ' %']);
    
    % TODO: Need a new way of calculating this:
    %disp(['    Average time spent over each ROI (ms): ' num2str( 1/length(handles.scanCoords)*1000*pathPeriod*sum(handles.pathObjNum>0)/length(handles.pathObjNum))]);
    
    % plot 
    figure
    subplot(3,2,1:2)
    plot(pathVelocity ./1e3)
    ylabel 'velocity (kV/s)'
    
    subplot(3,2,3:4)
    plot(pathAcceleration ./ 1e6)
    ylabel 'acceleration (MV/s^2)'
    
    subplot(3,2,5:6)
    plot([handles.pathObjNum handles.pathObjSubNum])
    ylabel 'object (sub) number'
        
% --- BUTTON Delete Path Element.
function pushDeletePathElement_Callback(hObject, eventdata, handles)
    % clears the selected item out of the scanCoords
%     if handles.drawingLine == true,
%         handles = StopLine(hObject, eventdata, handles);
%     end

   [groupIndex, scIndex, listboxIndex ] = getListboxItem(handles);
   
   if groupIndex > 0
       if scIndex == 0
           % this is a group
           handles = deleteGroup(handles, groupIndex);
       else
           % this is a line or a box
            set(handles.listboxScanCoords, 'Value', listboxIndex - 1);
            
            handles.group(groupIndex).scanCoords(scIndex) = [];
        end
   end
   
   % update here, otherwise, scanCoords gets overwritten when deleting first element
    guidata(hObject, handles);   % Update handles structure
    updatePath(hObject, eventdata, handles);
    pushButtonClearGraph_Callback(hObject, eventdata, handles);

% --- BUTTON pushButtonLoadImage
function pushButtonLoadImage_Callback(hObject, eventdata, handles)
    % loads in a data file, and sets:
    %   handles.im          (image for graph)
    %   handles.axisLimRow  (axis limits, row)
    %   handles.axisLimCol  (axis limits, column)
%     if handles.drawingLine == true,
%         handles = StopLine(hObject, eventdata, handles);
%     end
    
    % set(handles.editDwellTime,'String',num2str(1e6*d.Header.PixelClockSecs));   % set the string (sets in microseconds)    

    % guidata(hObject, handles);                                                  % Update handles structure before sending
    % frameHeight = str2double(d.Header.Frame_Height);
    % frameWidth = str2double(d.Header.Frame_Width);
    % xFrameOffset = str2double(d.Header.X_Frame_Offset);
    % yFrameOffset = str2double(d.Header.Y_Frame_Offset); 
    % magStr = d.Header.Magnification;

    % rotation is not saved here, because it is applied in MpScope3
                               
    
%         if handles.imageCh == 1
%             im = d.Ch1;
%         elseif handles.imageCh == 2
%             im = d.Ch2;
%         elseif handles.imageCh == 3
%             im = d.Ch3;
%         elseif handles.imageCh == 4
%             im = d.Ch4;
%         end
    
    if ~(evalin('base','exist(''MpData'',''var'')')),
        msgbox('No MpData available');
        return
    end%if ~(evalin('base','exist(''MpData'',''var'')')),
    
    % assigned data to handles structure, and update
    imageChStr = get(handles.editImageCh,'String');
    cmdStr = ['mean(MpData.framesCh',imageChStr,',3);'];
    handles.im = evalin('base',cmdStr);   % load the average of all frames  
%     handles.im = evalin('base','mean(MpData.frames,3)');   % load the average of all frames
    
    
    handles.axisLimRow = [evalin('base','MpData.startVoltageY'),evalin('base','MpData.stopVoltageY')];
    handles.axisLimCol = [evalin('base','MpData.startVoltageX'),evalin('base','MpData.stopVoltageX')];
    handles.pixelDwellTime =  evalin('base','MpData.pixelDwellTime');%seconds
    handles.pixelsToShift = [evalin('base','MpData.pixelsToShiftX'),evalin('base','MpData.pixelsToShiftY')];
    set(handles.editPixelDwellTime,'String',num2str(handles.pixelDwellTime*1e6,'%0.3f'));
    
    guidata(hObject, handles);             % Update handles structure (save the image)
    
    pushButtonClearGraph_Callback(hObject, eventdata, handles) % updates graph (draws image)
    
% --- BUTTON - Run Mask Path
function pushButtonRunMaskPath_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonRunMaskPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Ilya's code goes here
    % all relevant information should be in the handles data structure
    %               im:
    %               axisLimRow:
    %               axisLimCol:
    %               scanCoords:
    homefig=gcf;
    curdir=pwd;
    
    groupIndex = handles.currentGroupIndex;
    
    try
        cd([pwd '\..\mask_generation']);
        handles.group(groupIndex).scanCoords=maskwrapper(handles);
    catch 
		e=lasterror();
        disp(e.message)
        disp(e.stack)
        disp('Could not make mask!');
        cd(curdir);
        figure(homefig);
        rethrow(e);
    end
    cd(curdir);
    figure(homefig);
    pushButtonClearGraph_Callback(hObject, eventdata, handles); 
    guidata(hObject, handles);         % Update handles structure (save the image)
    updatePath(hObject, eventdata, handles);

% --- BUTTON - Move Up (moves scan item up list)
function pushButtonMoveUp_Callback(hObject, eventdata, handles)
    handles = moveListboxItem(handles, true);
    guidata(hObject, handles);
    
% --- BUTTON - Move Down (moves scan item down list)
function pushButtonMoveDown_Callback(hObject, eventdata, handles)
    handles = moveListboxItem(handles, false);         
    guidata(hObject, handles);
    
%% Other objects (LISTBOX)
    
% --- Executes on selection change in listboxScanCoords.
function listboxScanCoords_Callback(hObject, eventdata, handles)
% hObject    handle to listboxScanCoords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listboxScanCoords contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxScanCoords


% --- Executes during object creation, after setting all properties.
function listboxScanCoords_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%% TEXT BOXES - reads text into data structures on initialization and change


%--- EDIT (enter) - Image Ch
function editImageCh_Callback(hObject, eventdata, handles)
    handles.imageCh = str2double(get(hObject,'String'));
    guidata(hObject, handles);   % Update handles structure

% --- EDIT (creation) - Image Ch
function editImageCh_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    editImageCh_Callback(hObject, eventdata, handles)   % execute, to read initial value


% --- EDIT (enter) - Num Lines
function editNumLines_Callback(hObject, eventdata, handles)
    handles.nLines = str2double(get(hObject,'String'));
    guidata(hObject, handles);   % Update handles structure

% --- EDIT (creation) - Num Lines
function editNumLines_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end   
    editNumLines_Callback(hObject, eventdata, handles)   % execute, to read initial value

    
% --- EDIT (enter) - Num Spline Points
function editMaxAcc_Callback(hObject, eventdata, handles)
    handles.maxAcc = 1e6 * str2double(get(hObject,'String'));  % translate to V / s^2
    guidata(hObject, handles);   % Update handles structure

% --- EDIT (creation) - Num Spline Points
function editMaxAcc_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    editMaxAcc_Callback(hObject, eventdata, handles)   % execute, to read initial value
    
% --- EDIT (enter) - Num Box Lines
function editNumBoxLines_Callback(hObject, eventdata, handles)
    handles.numBoxLines = str2double(get(hObject,'String'));
    guidata(hObject, handles);   % Update handles structure
    
% --- EDIT (create) - Num Box Lines
function editNumBoxLines_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    editNumBoxLines_Callback(hObject, eventdata, handles)   % execute, to read initial value

% --- EDIT (enter) - Scan Velocity    
function editScanStepSize_Callback(hObject, eventdata, handles)
    %jd ... need to think about this ...
    handles.scanStepSize = 1e-3 * str2double(get(hObject,'String')); % entered in mV, convert to V
    guidata(hObject, handles);   % Update handles structure
    
% --- EDIT (create) - Scan Velocity
function editScanStepSize_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    editScanStepSize_Callback(hObject, eventdata, handles)   % execute, to read initial value
    
% --- EDIT (enter) - Draw Every Points
function editDrawEveryPoints_Callback(hObject, eventdata, handles)
    handles.drawEveryPoints = str2double(get(hObject,'String')); % entered in mV, convert to V
    guidata(hObject, handles);   % Update handles structure

% --- EDIT (create) - Draw Every Points
function editDrawEveryPoints_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    editDrawEveryPoints_Callback(hObject, eventdata, handles)   % execute, to read initial value
   

% --- BUTTON - Optimize Path

function pushButtonOptimize_Callback(hObject, eventdata, handles)
      % function uses random iterations to find a more correct path ...
    % different algorithms could be used, this one should be fairly
    % simple, fast, and transparent (and close to optimal
    disp('Computing distance matrix...');
    
    groupIndex = handles.currentGroupIndex;
    nScanCoords = length(handles.group(groupIndex).scanCoords);   % how many elements to optimize?
    
    if nScanCoords < 3
        % to short to optimize
        helpdlg('path too short to optimize!')
        return
    end

    % create a matrix, pathMath, which keeps track of distance from one object to another
    pathMat = zeros(nScanCoords,nScanCoords);   % matrix, holds distance from one point to another
    
    for i = 1:nScanCoords
      
       % find the path corresponding to object i
        
       % find indices of the start object
       allIndicesObjectI = find(handles.pathObjNum == i);
       firstIndexObjectI = allIndicesObjectI(1);
       lastIndexObjectI = allIndicesObjectI(end);      
        
       % cut out this portion of the path
       pathI = handles.path(firstIndexObjectI:lastIndexObjectI,:);
       centerI=allIndicesObjectI(ceil(end/2):(ceil(end/2)+1));
       
        for j = 1:nScanCoords
            if i == j 
                continue;    % path from a point to itself is irrelevant, leave at zero
            end
            
            % find the path corresponding to object j
        
            % find indices of the start object
            allIndicesObjectJ = find(handles.pathObjNum == j);
            firstIndexObjectJ = allIndicesObjectJ(1);
            lastIndexObjectJ = allIndicesObjectJ(end);      
        
            % cut out this portion of the path
            pathJ = handles.path(firstIndexObjectJ:lastIndexObjectJ,:);
            centerJ=allIndicesObjectJ(ceil(end/2):(ceil(end/2)+1));

            % have both paths now ... compute the distance!
             pathMat(i,j) = length(splineFuncMaxAcc(pathI,pathJ,handles.maxAcc,handles.pixelDwellTime));

        end
    end

    % run optimizations
    disp('Starting ant colony...');
  
    pathLengthPrevious = size(handles.path,1);

    pathMat(pathMat==0)=Inf;
    
    alpha=1;
    beta=5;
    rho=0.5;
    MaxITime=100;
    AntNum=length(pathMat);
    [GBTour] = ...
    AS(pathMat,AntNum,MaxITime,alpha,beta,rho); 
    opt_rte=GBTour(1:end-1);
    % [pospath ]=tsp_ant(pathMat,100,100);
    %  [opt_rte ]=tsp_ga([(1:length(pathMat))' (1:length(pathMat))'],pathMat,1000,1000,1:length(pathMat));
    %  figure(hFig);
    handles.group(groupIndex).scanCoords = handles.group(groupIndex).scanCoords(opt_rte);
    
    % done, genetic algorithm everything
    guidata(hObject, handles);                  % Update handles structure before sending
    updatePath(hObject, eventdata, handles);    % update the path to find the length
    handles=guidata(gca);
	gaPath=(size(handles.path,1));
    
    %starting orientation checking
    
    %temp=input('Do you want to optimize rotation  (y/n): ','s');  
    %if temp=='y',
    
    % changed this to a dialog box ...
    dlgAns = questdlg('Optimize rotations (ROI orientations) as well?','','yes','no','no');
    if strcmp(dlgAns,'yes')
        disp('Starting orientation checker...');
        handles.group(groupIndex).scanCoords = rotateoptimize(hObject,handles);
    end
    
    guidata(hObject, handles);                  % Update handles structure before sending
    updatePath(hObject, eventdata, handles);
    handles=guidata(gca);
    
    pathLengthNew=(size(handles.path,1));
	disp('Done...');
	disp(['Previous path = ' num2str(pathLengthPrevious)]);
    disp(['GA path = ' num2str(gaPath)]);
	disp(['New path = ' num2str(pathLengthNew)]);


% --- Executes on button press in pushButtonRandomizePath.
function pushButtonRandomizePath_Callback(hObject, eventdata, handles)
% hObject    handle to pushButtonRandomizePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    groupIndex = handles.currentGroupIndex;
    
    nScanCoords = length(handles.group(groupIndex).scanCoords);
    
    if nScanCoords < 3
        return           % nothing to shuffle
    end

    % make an array the same length as nScanCoords, and shuffle
    
    newOrder = [1 randperm(nScanCoords-1)+1];   % will start with 1 but be random after that
    handles.group(groupIndex).scanCoords = handles.group(groupIndex).scanCoords(newOrder); 
    updatePath(hObject,eventdata, handles);  



function editPixelDwellTime_Callback(hObject, eventdata, handles)
% hObject    handle to editPixelDwellTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPixelDwellTime as text
%        str2double(get(hObject,'String')) returns contents of editPixelDwellTime as a double
handles.pixelDwellTime = str2num(get(handles.editPixelDwellTime,'String'))*1e-6;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function editPixelDwellTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPixelDwellTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonNextLine.
function pushbuttonNextLine_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonNextLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = StopLine(hObject, eventdata, handles);
StartNewLine(hObject, eventdata, handles);
   
    

% --- Executes on button press in pushbuttonLineDone.
function pushbuttonLineDone_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLineDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = StopLine(hObject, eventdata, handles);




% --- Executes on button press in buttonLoadFakeData.
function buttonLoadFakeData_Callback(hObject, eventdata, handles)
% hObject    handle to buttonLoadFakeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    evalin('base','MpData.framesCh1 = zeros(500,500);');
    evalin('base','MpData.startVoltageX = -2;');
    evalin('base','MpData.stopVoltageX = +2;');
    evalin('base','MpData.startVoltageY = -2;');
    evalin('base','MpData.stopVoltageY = +2;');
    evalin('base','MpData.pixelDwellTime = 1e-6;');
    evalin('base','MpData.pixelsToShiftX = 77;');
    evalin('base','MpData.pixelsToShiftY = 77;');
    
    % assigned data to handles structure, and update
    handles.im = evalin('base','MpData.framesCh1;');   % load the fake image data
    handles.axisLimRow = [evalin('base','MpData.startVoltageX'),evalin('base','MpData.stopVoltageX')];
    handles.axisLimCol = [evalin('base','MpData.startVoltageY'),evalin('base','MpData.stopVoltageY')];
    handles.pixelDwellTime =  evalin('base','MpData.pixelDwellTime');%seconds
    handles.pixelsToShift = [evalin('base','MpData.pixelsToShiftX'),evalin('base','MpData.pixelsToShiftY')];
    set(handles.editPixelDwellTime,'String',num2str(handles.pixelDwellTime*1e6,'%0.3f'));
    guidata(hObject, handles);             % Update handles structure (save the image)
    pushButtonClearGraph_Callback(hObject, eventdata, handles); % updates graph (draws image)


    
function RefreshAllParameters(hObject, eventdata, handles)

handles.pixelDwellTime = str2num(get(handles.editPixelDwellTime,'String'))*1e-6; %translate to seconds
handles.scanStepSize = str2num(get(handles.editScanStepSize,'String'))*1e-3; % translate to Volts
handles.maxAcc = str2num(get(handles.editMaxAcc,'String'))*1e6; % translate to Volts / sec^2

guidata(hObject,handles);

    

function editPathLength_Callback(hObject, eventdata, handles)
% hObject    handle to editPathLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPathLength as text
%        str2double(get(hObject,'String')) returns contents of editPathLength as a double


% --- Executes during object creation, after setting all properties.
function editPathLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPathLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPathDuration_Callback(hObject, eventdata, handles)
% hObject    handle to editPathDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPathDuration as text
%        str2double(get(hObject,'String')) returns contents of editPathDuration as a double


% --- Executes during object creation, after setting all properties.
function editPathDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPathDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in textDiagnosticHidden.
function buttonDiagnosticHidden_Callback(hObject, eventdata, handles)
% hObject    handle to textDiagnosticHidden (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over textDiagnosticHidden.
function textDiagnosticHidden_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to textDiagnosticHidden (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
diagnosticWindowState = get(handles.panelHiddenDiagnostic,'Visible');
if isequal(diagnosticWindowState,'off'),
    set(handles.panelHiddenDiagnostic,'Visible','on');
else
    set(handles.panelHiddenDiagnostic,'Visible','off');
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
