function [B]=func_bubblevelocity(B, tstep, R, Z, minbubbledia_vel, ylim1, ylim2, cylgeometry)

% ----------------------------------------------------------------------
% this function adds velocity components to the matrix B
% note: bubble velocity in frame i is 0 if
% (a) total number of bubbles in frame i is not equal to total number of bubbles in frame i+1 (coalescence/splitting/eruption)
% (b) computed vx and vy are physically unreasonable 

% Smaller bubbles distort numbering; to improve linking, use larger values for minbubbledia
% increasing minbubledia_vel and choosing [ylim1, ylim2] to exclude small bubbles may improve linking 

% max bubble velocity constraints are based on observations in visualization
% the defult setting for vxmax, vxmax is that a bubble can travel atmost Width/10 in one time-step
% vymax based on bubble diameter
% -----------------------------------------------------------------------

Btemp = B;

% B = [frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2]
nframe1 = min(B(:,1)); B(:,1) = B(:,1) - nframe1 + 1;               % modifying starting frame to begin from 1
nframes = max(B(:,1)); 
m = length(B(:,1)); B = [linspace(1,m,m)',B]; 

% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2]
TF = B(:,9)>ylim1 & B(:,10)<ylim2 & B(:,6)>minbubbledia_vel; 
B = B(TF,:);                                                        % discard small bubbles to improve linking  
frame_nbubbles = histc(B(:,2),1:nframes);                           % # of bubbles in frame  
frame_linked = frame_nbubbles == circshift(frame_nbubbles, -1);     % equal # bubbles in consecutive frames 

B(:,7:14) = [];                             % B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia] 
m = length(B(:,1)); s1 = zeros(m,1); 
B = [linspace(1,m,m)', B, s1, s1, s1];      % B = [new-bubble#, orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, vx, vy, vz] 
B(:,11) = frame_linked(B(:,3));             % B(:,11) = 1 if frame linked to next frame
B = [B, s1]; 
B(:,12) = B(:,1)+frame_nbubbles(B(:,3));    % note: we don't care about non-linked frames

% B = [new-bubble#, orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, vx, vy, vz, linkedframe, linkedbubble] 
for i=1:length(B(:,1))  
    if B(i,11)>0 && B(i,3) ~= nframes 
        B(i,8) = (B(B(i,12),4)-B(i,4))/tstep; 
        B(i,9) = (B(B(i,12),5)-B(i,5))/tstep; 
        B(i,10)= (B(B(i,12),6)-B(i,6))/tstep; 
    end
end
B(:,1) = [];

% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, vx, vy, vz, linkedframe, linkedbubble] 
% max possible velocities are empirical and related to frequency of data sampling
if cylgeometry==1; vxmax = (2*R/10)/tstep; vzmax = vxmax; 
else; vxmax = (R/10)/tstep; vzmax = (Z/10)/tstep;  
end
vymax = 5*0.71*sqrt(9.81.*B(:,6));

TF = abs(B(:,7))>vxmax | abs(B(:,9))>vzmax | B(:,8)<0 | B(:,8)>vymax(:); B(TF,7:9) = 0;  
TF = B(:,7)==0 & B(:,8)==0 & B(:,9)==0; B(TF,:) = [];   % keep only non-zero elements

% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, vx, vy, vz, linkedframe, linkedbubble] 

if cylgeometry==1
    % convert to vr, vy, vtheta (cylindrical coord not used before to avoid 0/2pi issues 
    rmean = sqrt(B(:,3).^2 + B(:,5).^2); 
    thetamean = atan(B(:,5)./B(:,3));
    TF = sign(B(:,3)) == -1;                            % modifying for original theta in [pi/2,3pi/2]
    thetamean(TF) = thetamean(TF)+pi;
    TF = sign(thetamean) == -1;                         % modifying for original theta in [3*pi/2,2pi]
    thetamean(TF) = thetamean(TF)+2*pi;
    vr = B(:,7).*cos(thetamean) + B(:,9).*sin(thetamean); 
    vtheta = (-B(:,7).*sin(thetamean) + B(:,9).*cos(thetamean)); 
    B(:,7) = abs(vr); B(:,9) = abs(vtheta); 
end 

% B = [orig-bubble#, frame#, xmean, ymean, zmean, bubble-dia, vx, vy, vz, linkedframe, linkedbubble] 
% Btemp = [frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2]

% combine matrices Btemp and B 
m=length(Btemp(:,1)); Btemp = [Btemp, zeros(m,1), zeros(m,1), zeros(m,1)]; 
for i=1:length(B(:,1))
    Btemp(B(i,1),14) = B(i,7); Btemp(B(i,1),15) = B(i,8); Btemp(B(i,1),16) = B(i,9); 
end
B = Btemp;  
% B = [frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2, vr, vy, vtheta]

end



