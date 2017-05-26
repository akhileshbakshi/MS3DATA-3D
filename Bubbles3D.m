% variables are tailored for 'bubblestats2D.txt' and 'geometry.xlsx' provided 
% Major Edit 1: option lagrangetracking has been added where bubbles are individually tracked through frames 

clear all; clc; 

% 1. Mfix file properties 
nframes = 0;            % 0 to read from void fraction data 
                        % non-zero to specify exact number of frames 
epgcutoff = 0.65;       % consider only cells so that ep_g>epgcutoff
                        % note: epgcutoff<epgbubble 
ycutoff2   = 0.79;      % maximum domain y-extremity (ycutoff2 must be < ymax from simulation data) 
ycutoff1 =  0.003;      % minimum domain y-extremity
                        % domain extremeties are required to remove ambiguous bubbles (touching freeboard and distributor)                    
                        
% 2. Modify Geometry.xlsx and enter other simulation data here
cylgeometry = 1;        % set 1 if cylindrical reactor, 0 if Cartesian reactor 
cylcoord = 1;           % set 1 if cylindrical coordinates, 0 if cartesian coordinates 
                        % NOTE: cylcoord = 1 only possible if cylgeometry = 1                        
R = 0.075;              % R = bed radius if cylcoord = 1. R = bed diameter/width if cylcord = 0
Z = 0.0;                % Z = bed depth if cylgeometry = 0. leave unspecified otherwise
tstep = 0.01;           % time step of data sampling 

% 3. input/output files names 
% sample file provided has data corresponding to 150 frames 
bubblefile = 'bubblestats3D.txt';
printfile  = 'bubblestats3D'; 

% 4. criteria for bubble detection     
epgbubble = 0.7;        % threshold voidage for bubble (interphase) detection 
mincordlength = 0.01;   % discard small bubbles 
minCSlength = 0.01;     % discard small bubbles (infinite AR)
minbubbledia = 0.01;    % discard small bubbles      
ysmooth = 1;            % y-grid refinement 
xsmooth = 2;            % x-grid refinement. use xsmooth >= 2  if cylcoord = 1
zsmooth = 2;            % z-grid refinement. use zsmooth >= 2  if cylcoord = 1 
    
% 5. criteria for postprocessing of detected bubbles
ylim1 = 0.0;            % min y for postprocessing 
ylim2 = 0.6;            % max y for postprocessing 
rlim1 = 0;              % min r/x for postprocessing (if cylgeometry = 1, r=0 is centerline)
rlim2 = R;              % max r/x for postprocessing 
zlim1 = 0;              % min z for postprocessing (if cylgeometry = 1, leave unspecified)
zlim2 = 0;              % max z for postprocessing (if cylgeometry = 1, leave unspecified)
minbubbledia_vel = 0.01;% discard very small bubbles for bubble linking
diaratio = 1.1;         % maximum permissible ratio of bubble dia for linking  
dmax = 0.05;            % maximum permissible distance traveled by bubble in one time-step  
tolerance  = 0.0;       % minimum permissible bubble y-velocity = -tolerance x time-step 
lagrangetracking = 1;   % (recommended) set 1 to turn on lagrangian tracking of bubbles 
                        % linking is affected by bubble activity- splitting, coalescence and eruption
                        % for best linking results write data at high frequency 
                        % if globaltracking (lagrangetracking=0), consider increasing minbubledia_vel
                        % and choosing [ylim1, ylim2] to exclude small bubbles may improve linking 

% 6. Statistics for average computations 
nbinsax = 10;           % # bins for axial statistics between [ylim1, ylim2]
nbinsrad = 4;           % # bins for radial statistics between [rlim1, rlim2]

% ----------------------------------------------------------------
if cylgeometry == 0 && cylcoord == 1; error('cylcoord = 1 only possible if cylgeometry = 1'); end 
    
[nframes, bubblepropertiestotal] = func_bubbledetection(bubblefile, xsmooth, ysmooth, zsmooth, epgcutoff, epgbubble, mincordlength, minCSlength, minbubbledia, nframes, ycutoff1, ycutoff2, cylgeometry, cylcoord);
% bubblepropertiestotal = [frame#, xmean, ymean, zmean, bubbledia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2]
% if cylgeometry=1, xmean and zmean are in range [-radius,radius] 

bubblepropertiestotal = func_bubblevelocity(bubblepropertiestotal, tstep, minbubbledia_vel, ylim1, ylim2, cylgeometry, lagrangetracking, diaratio, dmax, tolerance); 
% bubblepropertiestotal_1 = [frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2, vx, vy, vz]
% if cylgeometry=1, vx = radial velocity, vz = theta velocity 

% ----------------------------------------------------------------
% sample for computing average statistics 
[bubblestats_2D, bubblestats_ax, bubblestats_rad]=func_bubblestatistics(bubblepropertiestotal_1, nbinsax, nbinsrad, ylim1, ylim2, rlim1, rlim2, zlim1, zlim2, cylgeometry);
% bubblestats_2D = [binr, biny, nb, vol-dia, CSmax, cord, AR1, AR2, nbubbles_linked, abs(vx), vy, abs(vz)]; 
% bubblestats_ax = [biny, nb_y, vol-dia, CSmax, cord, AR1, AR2, nbubbles_linked, abs(vx), vy, abs(vz)]; 
% bubblestats_rad= [binr, nb_r, vol-dia, CSmax, cord, AR1, AR2, nbubbles_linked, abs(vx), vy, abs(vz)]; 
% if cylgeometry=1, vx = radial velocity, vz = theta velocity 

% ----------------------------------------------------------------
% % sample for writing to files 
% filename = strcat(printfile,'_BubbleStats_Ax.txt');
% dlmwrite(filename,bubblestats_ax,'delimiter',' ','precision',4); 



    

      






