function [B]=func_bubblevelocity(B, tstep, minbubbledia_vel, ylim1, ylim2, cylgeometry, lagrangetracking, diaratio, dmax, tolerance)

% ----------------------------------------------------------------------
% this function adds velocity components as columns to matrix B
% in order to avoid ambiguity because of bubble coalescence, splitting and eruption
% a bubble in frame t is linked with a bubble in frame t+1 iff: 
% 1. distance between bubbles < dmax
% 2. change in bubble dia < diaratio
% 3. y-location of bubble in t+1 > y-location of bubble in t - tolerance 
% for lagrangetracking = 0 case (not recommended), addition filter of 
% equal number of bubbles in frames t & t+1 is implemented 

% To improve linking, remove very small (fictitious?) bubbles by increasing minbubbledia_vel 
% ----------------------------------------------------------------------

Bcopy = B;
m = length(B(:,1)); B = [(1:m)',B]; 
% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2]

% discard (a) small bubbles and (b) bubbles outside domain of interest   
TF = B(:,9)>ylim1 & B(:,10)<ylim2 & B(:,6)>minbubbledia_vel; B = B(TF,:);

% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2]
if ~lagrangetracking; B = func_globaltracking(B, tstep, cylgeometry, diaratio, dmax, tolerance);  
else; [bubbletrace, B] = func_lagrangetracking(B, tstep, cylgeometry, diaratio, dmax, tolerance); 
end
% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2, vx, vy, vz]
% bubbletrace is matrix of tracked bubbles; each column represents a trace
% (or series of bubble numbers linked) 
% properties of ith bubble in jth trace can be accesssed using B(bubbletrace(i,j),:) 

% Overlay B on Bcopy 
m = length(Bcopy(:,1)); Bcopy = [(1:m)', Bcopy, zeros(m,1), zeros(m,1), zeros(m,1)]; 
[Lia,Locb] = ismember(B(:,1),Bcopy(:,1)); 
Bcopy(Locb,:) = B; 

B = Bcopy; B(:,1) = []; 
% B = [frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2, vx, vy, vz]
end



function [B]=func_globaltracking(B, tstep, cylgeometry, diaratio, dmax, tolerance)

% ----------------------------------------------------------------------
% this function links bubbles based on global indexing 
% ----------------------------------------------------------------------

nframe1 = min(B(:,2)); B(:,2) = B(:,2) - nframe1 + 1;               % modifying starting frame to begin from 1
nframes = max(B(:,2)); 
frame_nbubbles = histc(B(:,2),1:nframes);                           % # of bubbles in frame  
frame_linked = frame_nbubbles == circshift(frame_nbubbles, -1);     % equal # bubbles in consecutive frames 

% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2]
m = length(B(:,1)); s1 = zeros(m,1); 
B = [linspace(1,m,m)', B(:,1:6), s1, s1, s1, s1, s1, B(:,7:14)];          
B(:,11) = frame_linked(B(:,3));             % B(:,11) = 1 if frame linked to next frame 
B(:,12) = B(:,1)+frame_nbubbles(B(:,3));    % note: we don't care about non-linked frames

% B = [new-bubble#, orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, vx, vy, vz, linkedframe, linkedbubble, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2]
for i=1:length(B(:,1))  
    if B(i,11)>0 && B(i,3) ~= nframes 
        if B(B(i,12),7)/B(i,7)>=1/diaratio && B(B(i,12),7)/B(i,7)<=diaratio     % diameter based condition 
            B(i,8:10) = (B(B(i,12),4:6)-B(i,4:6))/tstep; 
        end
    end
end
B(:,11:12) = []; B(:,1) = [];

% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, vx, vy, vz, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2] 
TF = sqrt((B(:,7)*tstep).^2+(B(:,8)*tstep).^2+(B(:,9)*tstep).^2) > dmax & B(:,7)<= -tolerance*tstep;    % distance based condition 
B(TF,7:9) = 0;

if cylgeometry==1
    % convert to vr, vy, vtheta (cylindrical coord not used before to avoid 0/2pi issues)
    rmean = sqrt(B(:,3).^2 + B(:,5).^2); 
    thetamean = atan(B(:,5)./B(:,3));
    TF = sign(B(:,3)) == -1;                            % modifying for original theta in [pi/2,3pi/2]
    thetamean(TF) = thetamean(TF)+pi;
    TF = sign(thetamean) == -1;                         % modifying for original theta in [3*pi/2,2pi]
    thetamean(TF) = thetamean(TF)+2*pi;
    vr = B(:,7).*cos(thetamean) + B(:,9).*sin(thetamean); 
    vtheta = (-B(:,7).*sin(thetamean) + B(:,9).*cos(thetamean)); 
    B(:,7) = vr; B(:,9) = vtheta; 
end 

B = [B(:,1:6), B(:,10:17), B(:,7:9)];
B(:,2) = B(:,2) + nframe1 - 1;  
% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2, vx, vy, vz] 
end

function [bubbletrace, B] = func_lagrangetracking(B, tstep, cylgeometry, diaratio, dmax, tolerance)

% ----------------------------------------------------------------------
% this function tracks bubbles one by one 
% ----------------------------------------------------------------------

% send orig-bubble# to last column to avoid changing convention in func_bubbletrace
m = zeros(length(B(:,1)),1); 
B = [B(:,2:14), m, m, m, B(:,1)]; 
% B = [frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2, vx, vy, vz, orig-bubble#]
diacol = 5;             % bubble-dia column in bubbleproperties: 4 for 2D case, 5 for 3D case
coordloc = 2:4;         % coordinate columns in bubbleproperties: 2:3 for 2D case, 2:4 for 3D case

[B, bubbletrace] = func_bubbletrace(B, diacol, coordloc, diaratio, dmax, tolerance);

[m n] = size(bubbletrace); 
for i=1:length(B(:,1))
    if sum(sum(bubbletrace == i)) == 1        
        [row, col] = ind2sub([m,n], find(bubbletrace == i));    % if bubble exists in bubble trace 
        if row ~= m
            linkedbubble = bubbletrace(row+1,col);              % if bubble not at end of its trace
            if linkedbubble ~= 0                                % if bubble not at end of its trace
                B(i,14:16) = (B(linkedbubble,2:4)-B(i,2:4))/tstep; 
            end        
        end
    end
end    

% modify trace to have actual bubble numbers
% keep separate from velocity loop; max to min to avoid numbering conflicts
for i=length(B(:,1)):-1:1
    if sum(sum(bubbletrace == i)) == 1        
        [row, col] = ind2sub([m,n], find(bubbletrace == i)); 
        bubbletrace(row,col) = B(i,17); 
    end
end

B = [B(:,17), B(:,1:5), B(:,14:16), B(:,6:13)];
% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, vx, vy, vz, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2] 

if cylgeometry==1
    % convert to vr, vy, vtheta (cylindrical coord not used before to avoid 0/2pi issues)
    rmean = sqrt(B(:,3).^2 + B(:,5).^2); 
    thetamean = atan(B(:,5)./B(:,3));
    TF = sign(B(:,3)) == -1;                            % modifying for original theta in [pi/2,3pi/2]
    thetamean(TF) = thetamean(TF)+pi;
    TF = sign(thetamean) == -1;                         % modifying for original theta in [3*pi/2,2pi]
    thetamean(TF) = thetamean(TF)+2*pi;
    vr = B(:,7).*cos(thetamean) + B(:,9).*sin(thetamean); 
    vtheta = (-B(:,7).*sin(thetamean) + B(:,9).*cos(thetamean)); 
    B(:,7) = vr; B(:,9) = vtheta; 
end 

B = [B(:,1:6), B(:,10:17), B(:,7:9)];
% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2, vx, vy, vz]
end


function [B, bubbletrace] = func_bubbletrace(B, diacol, coordloc, diaratio, dmax, tolerance)

% ----------------------------------------------------------------------
% this function tracks bubbles across consecutive frames 
% output B = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR1, vx, vy, orig-bubble#]
% output bubbletrace = every column is a trace of the bubble as it moves through 
% NOTE: trace ends if bubble has split, coalesced or erupted 
% ----------------------------------------------------------------------

% B = [frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2, vx, vy, vz, orig-bubble#]

% remove isolated frames where linking is not possible  
Bframes = unique(B(:,1)); 
TF1 = Bframes(:,1)-circshift(Bframes(:,1),1)  == 1;        % bubbles in frame before
TF2 = circshift(Bframes(:,1),-1)-Bframes(:,1) == 1;        % bubbles in frame after 
TF = TF1 | TF2; 
Bframes = Bframes(find(TF)); 
TF = ismember(B(:,1),Bframes); B = B(TF,:); 

% starting frame numbers from 1. We do not need these for bubbletrace, so it is okay!
nframes1 = min(B(:,1)); B(:,1) = B(:,1)-nframes1+1; 

% there are two cases: 
% 1. no bubbles in next frame: stop linking; reset framehistory = 0; 
% 2. bubbles in next frame: framehistory = min(framehistory + 1,3); business as usual 

Bframes = unique(B(:,1));                                   % set of frame numbers 
TF = circshift(Bframes(:,1),-1)-Bframes(:,1) == 1;           
Bframes = [Bframes, TF];                                    % [frame number, if bubbles in next frame] 
nframes = length(Bframes(:,1)); 

% initiate bubbletrace matrix 
framehistory = 0;
frame1 = find(B(:,1) == Bframes(1,1)); nbubbles = length(frame1);
bubbletrace = zeros(nframes,nbubbles);
bubbletrace(1,:) = frame1; 
framehistory = framehistory + 1; 

% since isolated frames have been removed, frame 2 necessarily has bubbles 
for i = 2:nframes
    
    if framehistory == 0 
        frame = find(B(:,1) == Bframes(i,1)); 
        bubbletrace = func_bubbletrace0(bubbletrace, i, frame);             % add new trace 
        if Bframes(i,2)== 0; framehistory = 0; 
        elseif framehistory < 3; framehistory = framehistory + 1;           % if framehistory = 3, keep unchanged
        end       
        
    elseif framehistory == 1
        bubblesconsider = bubbletrace(i-1,:) > 0;
        frame1 = bubbletrace(i-1,bubblesconsider);                          % #s of bubbles in previous frame
        frame2 = find(B(:,1) == Bframes(i,1));                              % #s of bubbles in current frame 
        
        TF1 = bubbletrace(i-1,:) > 0;                         
        predictcoord = zeros(length(frame1),coordloc(end)-1);               % predict coordinates based on extrapolation 
        predictcoord(TF1(bubblesconsider),:) = func_bubbleextrapolate(B, bubbletrace(i-1:i-1,TF1), 1, coordloc);
        [bubbletrace, frame2] = func_linkbubbles (bubbletrace, i, B, frame1, frame2, predictcoord, coordloc, diaratio, diacol, dmax, tolerance); 
        bubbletrace           = func_bubbletrace0(bubbletrace, i, frame2);  % add new trace with bubbles not linked 
        if Bframes(i,2)== 0; framehistory = 0; 
        elseif framehistory < 3; framehistory = framehistory + 1;           % if framehistory = 3, keep unchanged 
        end        
        
    elseif framehistory == 2
        bubblesconsider = bubbletrace(i-1,:) > 0;
        frame1 = bubbletrace(i-1,bubblesconsider);                          % #s of bubbles in previous frame
        frame2 = find(B(:,1) == Bframes(i,1));                              % #s of bubbles in current frame

        TF1 = bubbletrace(i-2,:) == 0 & bubbletrace(i-1,:) > 0;                          
        TF2 = bubbletrace(i-2,:) >  0 & bubbletrace(i-1,:) > 0; 
        predictcoord = zeros(length(frame1),coordloc(end)-1);               % predict coordinates based on extrapolation 
        predictcoord(TF1(bubblesconsider),:) = func_bubbleextrapolate(B, bubbletrace(i-1:i-1,TF1), 1, coordloc);
        predictcoord(TF2(bubblesconsider),:) = func_bubbleextrapolate(B, bubbletrace(i-2:i-1,TF2), 2, coordloc);
        [bubbletrace, frame2] = func_linkbubbles (bubbletrace, i, B, frame1, frame2, predictcoord, coordloc, diaratio, diacol, dmax, tolerance); 
        bubbletrace           = func_bubbletrace0(bubbletrace, i, frame2);  % add new trace with bubbles not linked 
        if Bframes(i,2)== 0; framehistory = 0; 
        elseif framehistory < 3; framehistory = framehistory + 1;           % if framehistory = 3, keep unchanged 
        end        
    
    elseif framehistory == 3 
        bubblesconsider = bubbletrace(i-1,:) > 0;
        frame1 = bubbletrace(i-1,bubblesconsider);                          % #s of bubbles in previous frame
        frame2 = find(B(:,1) == Bframes(i,1));                              % #s of bubbles in current frame
        
        TF1 = bubbletrace(i-2,:) == 0 & bubbletrace(i-1,:) > 0;             % predict coordinates based on extrapolation                      
        TF2 = bubbletrace(i-3,:) == 0 & bubbletrace(i-2,:) > 0 & bubbletrace(i-1,:) > 0; 
        TF3 = bubbletrace(i-3,:) >  0 & bubbletrace(i-2,:) > 0 & bubbletrace(i-1,:) > 0; 
        predictcoord = zeros(length(frame1),coordloc(end)-1);               % add new trace with bubbles not linked 
        predictcoord(TF1(bubblesconsider),:) = func_bubbleextrapolate(B, bubbletrace(i-1:i-1,TF1), 1, coordloc);
        predictcoord(TF2(bubblesconsider),:) = func_bubbleextrapolate(B, bubbletrace(i-2:i-1,TF2), 2, coordloc);
        predictcoord(TF3(bubblesconsider),:) = func_bubbleextrapolate(B, bubbletrace(i-3:i-1,TF3), 3, coordloc);
        [bubbletrace, frame2] = func_linkbubbles (bubbletrace, i, B, frame1, frame2, predictcoord, coordloc, diaratio, diacol, dmax, tolerance);
        bubbletrace           = func_bubbletrace0(bubbletrace, i, frame2);
        if Bframes(i,2)== 0; framehistory = 0; 
        elseif framehistory < 3; framehistory = framehistory + 1;           % if framehistory = 3, keep unchanged 
        end  
    end
end

B(:,1) = B(:,1)+nframes1-1;
bubbletrace = func_bubbletracecondense(bubbletrace);                        % remove unlinked bubbles in bubbletrace 
end

function bubbletrace = func_bubbletrace0 (bubbletrace, current, frame2)

% ----------------------------------------------------------------------
% this function adds new trace with bubbles which were unable to be linked with previous frame 
% ----------------------------------------------------------------------

n = length(frame2);
if n > 0
    [x y] = size(bubbletrace);
    bubbletrace = [bubbletrace zeros(x,n)];
    bubbletrace(current,(y+1):(y+n)) = frame2;                          
end
end

function newlocations = func_bubbleextrapolate (B, prevlocatons, history, coordloc)

% ----------------------------------------------------------------------
% this function predicts location of bubbles in frame t+1 based on
% locations in frames t, t-1 and/or t-2 (depending on history available) 
% prevlocatons = bubble numbers in bubbletrace matrix
% newlocatons = array of extrapolated x, y (and z) locations of the bubbles
% ----------------------------------------------------------------------

if history == 1
    % No extrapolation
    newlocations = B(prevlocatons(1,:),coordloc); 
end
if history == 2
    % Linear extrapolation x = at + b. (0,x1), (1,x2) are solutions
    % => a = (x1-2x2+x3)/2, b = (-3x1+4x2-x3)/2, c=x1 
    a = B(prevlocatons(2,:),coordloc) - B(prevlocatons(1,:),coordloc); 
    b = B(prevlocatons(1,:),coordloc); 
    newlocations = 2*a + b;
end
if history == 3
    % Quadtratic extrapolation x = at^2 + bt + c. (0,x1), (1,x2), (2,x3) are solutions
    % => a = (x1-2x2+x3)/2, b = (-3x1+4x2-x3)/2, c=x1 
    a =  0.5*(   B(prevlocatons(1,:),coordloc) - 2*B(prevlocatons(2,:),coordloc) + B(prevlocatons(3,:),coordloc));
    b =  0.5*(-3*B(prevlocatons(1,:),coordloc) + 4*B(prevlocatons(2,:),coordloc) - B(prevlocatons(3,:),coordloc));
    c = B(prevlocatons(1,:),coordloc);
    newlocations = 9*a+3*b+c; 
end 
end

function [bubbletrace, frame2] = func_linkbubbles (bubbletrace, current, B, frame1, frame2, predictcoord, coordloc, diaratio, diacol, dmax, tolerance)

% ----------------------------------------------------------------------
% this function links bubbles in frame t and t+1 
% linking is based on predicting locations of bubbles in t and comparing with 
% actual locations of bubbles in t+1 
% bubbletrace = array updated with bubbles linked from current frame
% frame2 = bubbles that could not be linked with bubbles in frame1
% ----------------------------------------------------------------------

% distance between predicted and actual bubbles in frame2     
m = length(frame1); n = length(frame2); distance = zeros(m,n);
for ivar = 1:m
    for jvar = 1:n
        distance(ivar,jvar) = sqrt(sum((predictcoord(ivar,:)-B(frame2(jvar),coordloc)).^2));
    end
end

% link bubbles (going in icreasing order of distance)  
flag = true;
while flag   
    [M, I] = min(distance(:)); [row, col] = ind2sub([m,n], I);    
    if M >= dmax; flag = false; 
    else
        tempdiaratio = B(frame2(col),diacol)/B(frame1(row),diacol); 
        if tempdiaratio < 1; tempdiaratio = 1/tempdiaratio; end
        if frame2(col) > 0 & tempdiaratio <= diaratio & B(frame2(col),3)-B(frame1(row),3)>=-tolerance  % account for grid resolution
        % check that bubble not linked yet, diameter within diaratio, positive y location
        index = find(bubbletrace(current-1,:) == frame1(row));
        bubbletrace(current,index) = frame2(col); frame2(col) = 0;
        end
    end   
    % remove linked bubbles from search
    distance(row,:) = dmax; distance(:,col) = dmax;
end

% remove bubbles that have been linked
TF = frame2 == 0; frame2(TF) = [];
end

function bubbletrace = func_bubbletracecondense (bubbletrace)

% ----------------------------------------------------------------------
% this function condenses bubbletrace array by removing unlinked bubbles 
% ----------------------------------------------------------------------

TF = sum(bubbletrace ~= 0) == 1;   bubbletrace(:,TF) = [];                  % no links- remove columns with only one element

[m, n] = size(bubbletrace);
for i = 1:n
    TF = find(bubbletrace(:,i));                                            % making first entry in column non-zero
    bubbletrace(:,i) = circshift(bubbletrace(:,i),1-TF(1));                 
end

TF = sum(transpose(bubbletrace)) == 0; bubbletrace(TF,:) = [];              % remove empty rows
end
