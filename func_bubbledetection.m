function [nframes, bubblepropertiestotal] = func_bubbledetection(bubblefile, xsmooth, ysmooth, zsmooth, epgcutoff, epgbubble, mincordlength, minCSlength, minbubbledia, nframes, ycutoff1, ycutoff2, cylgeometry, cylcoord)
  
bubblepropertiestotal = [0 0 0 0 0 0 0 0 0 0 0 0 0];

fileID = fopen(bubblefile);
Anet = textscan(fileID, '%f %f %f');
fclose(fileID);
Anet = [Anet{1,1} Anet{1,2} Anet{1,3}];

frame1 = Anet(1,1);                                     % begin frame numbering from 1 
Anet(:,1) = Anet(:,1) - frame1 + 1; 

% script to compute nframes based on bubblefile 
% (executed only if nframes = 0 set by user)
Anet = [Anet; 100000, 100000, 100000];
[s1,s2] = ismember(unique(Anet(:,1)),Anet(:,1));
if nframes==0
    nframes = length(unique(Anet(:,1)))-1;              % to account for the extra 100000 at the end 
end
frameloc = s2;                                          % location of beginning of framei 

[R,nx,H,ny,Z,nz,coarsegridglobal] = func_readgeometry(cylcoord);    
% coarsegridglobal = [cell#, x, y, z];

% ----------------------------------------------------------------
% time loop begins 

parfor framei = 1:nframes

    framei
   
A1 = Anet(frameloc(framei):frameloc(framei+1)-1,2);     % cell-indices  
A2 = Anet(frameloc(framei):frameloc(framei+1)-1,3);     % ep_g 
epgcoarse = epgcutoff*ones(nx*ny*nz,1); 
epgcoarse(A1) = A2;                                     % replacing ep_g where ep_g>epgcutoff

coarsegrid = coarsegridglobal; 
xcoarse = coarsegrid(:,2); 
ycoarse = coarsegrid(:,3); 
zcoarse = coarsegrid(:,4);

% Converting (xcoarse,zcoarse) from cylindrical to cartesian coordinates
if cylcoord == 1
    xtemp = xcoarse.*cos(zcoarse); ztemp = xcoarse.*sin(zcoarse);  
    xcoarse = xtemp; zcoarse = ztemp; 
    xlimits_new = [-R,R]; 
    zlimits_new = [-R,R]; 
    nz1 = nx;
else
    xlimits_new = [0,R]; 
    zlimits_new = [0,Z]; 
    nz1 = nz;    
end
          
xgrid = linspace(xlimits_new(1),xlimits_new(2),xsmooth*nx); 
ygrid = linspace(0,H,ysmooth*ny);
zgrid = linspace(zlimits_new(1),zlimits_new(2),zsmooth*nz1);                               

[xgrid,ygrid, zgrid] = meshgrid(xgrid,ygrid,zgrid);
xgrid = reshape(xgrid,[],1);
ygrid = reshape(ygrid,[],1);
zgrid = reshape(zgrid,[],1);
 
epggrid = griddata(xcoarse,ycoarse,zcoarse,epgcoarse,xgrid,ygrid,zgrid);
% interpolating data to cartesian and/or finer grid 
% can be improved for better code efficiency 

% ----------------------------------------------------------------
% restructuring matrix

deltax = (xlimits_new(2)-xlimits_new(1))/(xsmooth*nx-1);          % in fine grid 
deltay = H/(ysmooth*ny-1); 
deltaz = (zlimits_new(2)-zlimits_new(1))/(zsmooth*nz1-1); 

B = [xgrid ,ygrid ,zgrid, epggrid]; 
B(isnan(B)) = 0;        % boundary cells have NaN (no interpolation) 
B = sortrows(B,2);      % sortrows based on axial location 


% ----------------------------------------------------------------
% checking neighbouring cells for linking 
% B = [# x y z epg]

m = length(B(:,1)); 
v = zeros(m,1);          
B = [linspace(1,m,m)',B, v, v, v];   % numbering & 2 columns to store neighbours

B(:,6) = B(:,5)>epgbubble & circshift(B(:,5),-1)>epgbubble;                         % right neighbour
B(:,6) = (B(:,1)+1).* B(:,6);

B(:,7) = B(:,5)>epgbubble & circshift(B(:,5),-xsmooth*nx)>epgbubble;                % back neighbour
B(:,7) = (B(:,1)+xsmooth*nx).* B(:,7);

B(:,8) = B(:,5)>epgbubble & circshift(B(:,5),-(xsmooth*nx*zsmooth*nz1))>epgbubble;   % top neighbour
B(:,8) = (B(:,1)+((xsmooth*nx*zsmooth*nz1))).* B(:,8);

TF = B(:,7)==0;             % fill B(:,6) before B(:,7) and before B(:,8)
s = TF.*B(:,8); 
B(:,7)=B(:,7)+s; 
B(:,8)=B(:,8)-s; 

TF = B(:,6)==0;             % fill B(:,6) before B(:,7) and before B(:,8)
s = TF.*B(:,7); 
B(:,6)=B(:,6)+s; 
B(:,7)=B(:,7)-s; 

TF = B(:,7)==0;             % repeat for the initial cases 0,non-0, non-0
s = TF.*B(:,8); 
B(:,7)=B(:,7)+s; 
B(:,8)=B(:,8)-s; 

% ----------------------------------------------------------------
% checking if bubble is at interface (1=interface 0=interior)
% B = [# x y z epg neighbour#1 neighbour#2 neighbour#3]

v = zeros(length(B(:,1)),1); 
B = [B, v]; 

s = B(:,5);
TF1 = s>epgbubble;            
TF2 = circshift(s,1)>epgbubble;                         % left neighbour
TF3 = circshift(s,-1)>epgbubble;                        % right neighbour
TF4 = circshift(s,xsmooth*nx)>epgbubble;                % back neighbour
TF5 = circshift(s,-xsmooth*nx)>epgbubble;               % front neighbour 
TF6 = circshift(s,xsmooth*nx*zsmooth*nz1)>epgbubble;    % bottom neighbour
TF7 = circshift(s,-xsmooth*nx*zsmooth*nz1)>epgbubble;   % top neighbour 

B(:,9) = TF1 & TF2 & TF3 & TF4 & TF5 & TF6 & TF7; 
B(:,9) = 1-B(:,9); 

TF = B(:,5)>epgbubble;                        % delete rows with epg<epgbubble
B = B(TF,:); 


% ----------------------------------------------------------------
% renumbering indices of points 
% B = [# x y z epg neighbour#1 neighbour#2 neighbour#3 interface]

[s1,s2] = ismember(B(:,6),B(:,1)); B(:,6) = s2; 
[s1,s2] = ismember(B(:,7),B(:,1)); B(:,7) = s2; 
[s1,s2] = ismember(B(:,8),B(:,1)); B(:,8) = s2; 

B(:,1) = []; B(:,4) = [];   

% ----------------------------------------------------------------
% linking bubbles
% B = [x y z neighbour#1 neighbour#2 neighbour#3 interface]

v = zeros(length(B(:,1)),1); B = [B, v, v, v];           
s = B(:,7); B(:,7)=[]; B = [B, s];
% B = [x y z neighbour#1 neighbour#2 neigbour#3 bubble#1 bubble#2 bubble#3 interface]                         

% start numbering bubbles 
bubblenumctr = 0;           % keeps note of latest bubblenum

for i=1:length(B(:,1))   
    if B(i,7)>0 
        if B(i,4)>0 
            if B(B(i,4),7) ==0 
                B(B(i,4),7)= B(i,7); 
            elseif B(B(i,4),7)~=B(i,7) && B(B(i,4),8)==0
                B(B(i,4),8)=B(i,7); 
            elseif B(B(i,4),7)~=B(i,7) && B(B(i,4),8)~=B(i,7) && B(B(i,4),9)==0
                B(B(i,4),9)=B(i,7); 
            end
               
            if B(i,5)>0                     % because B(i,5)>0 only possible if B(i,4)>0
                if B(B(i,5),7) ==0 
                    B(B(i,5),7)= B(i,7); 
                elseif B(B(i,5),7)~=B(i,7) && B(B(i,5),8)==0
                    B(B(i,5),8)=B(i,7); 
                elseif B(B(i,5),7)~=B(i,7) && B(B(i,5),8)~=B(i,7) && B(B(i,5),9)==0
                    B(B(i,5),9)=B(i,7); 
                end
                
                if B(i,6)>0
                    B(B(i,6),7)= B(i,7); 
                end
            end
        end
    
    elseif B(i,4)>0 && B(B(i,4),7)>0        
        B(i,7) = B(B(i,4),7); 
        if B(i,5)>0                     % because B(i,5)>0 only possible if B(i,4)>0
           if B(B(i,5),7) ==0 
                B(B(i,5),7)= B(i,7); 
           elseif B(B(i,5),7)~=B(i,7) && B(B(i,5),8)==0
                B(B(i,5),8)=B(i,7); 
           elseif B(B(i,5),7)~=B(i,7) && B(B(i,5),8)~=B(i,7)  && B(B(i,5),9)==0
                B(B(i,5),9)=B(i,7); 
           end

            if B(i,6)>0
                B(B(i,6),7)= B(i,7); 
            end
        end

    elseif B(i,5)>0 && B(B(i,5),7)>0  
        B(i,7) = B(B(i,5),7);
        B(B(i,4),7)= B(i,7);                % because B(i,4)>0 && B(B(i,4),7)==0 by loop construction 
        if B(i,6)>0
            B(B(i,6),7)= B(i,7); 
        end
                     
    else                                    % B(i,5)=0 and (B(i,3)=0 || B(B(i,3),5)=0)
        bubblenumctr = bubblenumctr+1; 
        B(i,7) = bubblenumctr; 
        if B(i,4)>0      
            B(B(i,4),7)= B(i,7);
            if B(i,5)>0
                B(B(i,5),7)= B(i,7); 
                if B(i,6)>0
                    B(B(i,6),7)= B(i,7);
                end
            end
        end
    end 
end



% ----------------------------------------------------------------
% find unique dispute combinations (i.e. B(i,7)~=B(i,8)~=B(i,9))
% B = [x y z neighbour#1 neighbour#2 neigbour#3 bubble#1 bubble#2 bubble#3 interface]                         

dispute = [B(:,7:9), B(:,8) > 0]; 
% dispute = [bubble#1 bubble#2 bubble#3 dispute]
TF = dispute(:,end)<1; dispute(TF,:)=[];    % keep only dispute cases

% get unqiue ID of combination (assume max 9999 disputes allowed) 
dispute(:,4) = []; 
dispute = sort(dispute,2,'descend'); 

% remove the third column and add combinations at the end 
TF = dispute(:,3) ~= 0;
s1= dispute(TF,2);  
s2= dispute(TF,3);
dispute(:,end) = []; 
dispute = [dispute; s1, s2]; 

dispute(:,3) = dispute(:,1)*10^4 + dispute(:,2);

disputecombination = zeros(length(unique(dispute(:,3))),3); 
disputecombination(:,1) = mod(unique(dispute(:,3)),10^4);
disputecombination(:,2) = (unique(dispute(:,3))-disputecombination(:,1))/10^4; 
disputecombination=sortrows(disputecombination,2); 

% checking for possibility that c->a and c->b 

[m n] = size(disputecombination);
if m>1          % consider cases with >1 disputes only 
    % add a fictitious row in dipute combination for cases where only two
    % entries [a b; a c] 
    disputecombination = [zeros(1,3); disputecombination];    
    s = disputecombination(:,2); 
    % disputecombination = [disputecombination, zeros(length(s),1)];
    disputecombination(:,3) = s == circshift(s,1);
    conflicts = sum(disputecombination(:,3));
    disputecombination(1,:) = []; 
    
    while conflicts>0  
        [s1, s2] = ismember(1,disputecombination(:,3));
        disputecombination(s2,2)= max(disputecombination(s2,1),disputecombination(s2-1,1));
        disputecombination(s2,1)= min(disputecombination(s2,1),disputecombination(s2-1,1));

        s1 = disputecombination(:,1)*10^4 + disputecombination(:,2);

        disputecombinationtemp = zeros(length(disputecombination(:,1)),3); 
        disputecombinationtemp(1:length(unique(s1)),2) = mod(unique(s1),10^4);
        disputecombinationtemp(1:length(unique(s1)),1) = (unique(s1)-disputecombinationtemp(1:length(unique(s1)),2))/10^4; 
        if (disputecombinationtemp(end,1)==0) 
            disputecombinationtemp(end,:)=[];
        end           
        disputecombinationtemp=sortrows(disputecombinationtemp,2); 

        s = disputecombinationtemp(:,2); 
        disputecombinationtemp(:,3) = s == circshift(s,1);
        conflicts = sum(disputecombinationtemp(:,3));
        disputecombination = disputecombinationtemp; 
        disputecombinationtemp = zeros(length(disputecombination(:,1)),3);
    end
end

% ----------------------------------------------------------------
% modify bubble# for disputed cases
% B = [x y z neighbour#1 neighbour#2 neighbor#3 bubble#1 bubble#2 bubble#3 interface]
% dispute combination = [newvalue oldvalue conflict] 

% change B from highest dispute case e.g. 19->10 & 10->5 => 19->5 
 
for i=length(disputecombination(:,1)):-1:1
    s1 = B(:,7);
    s1(s1==disputecombination(i,2))=disputecombination(i,1);
    numbermod = find(s1>disputecombination(i,2));            % renumbering the higher bubble    
    s1(numbermod) = s1(numbermod)-1; 
    B(:,7) = s1;  
end 

B(:,8:9)=[];          % bubble#2 and bubble#3 deleted      
B = sortrows(B,7); B(:,4:6)=[];        


% ----------------------------------------------------------------
% compute bubble properties 
% B = [x y z bubble# interface]
bubbleproperties = func_bubbleproperties(B(:,1:4), deltax, deltay, deltaz); 
% bubbleproperties = [bubble#, xmean, ymean, zmean, volume, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2]
bubbleproperties(:,5) = (6*bubbleproperties(:,5)/pi).^(1/3); 
% filters to eliminate small bubbles 
TF1 = bubbleproperties(:,9)-bubbleproperties(:,8)<mincordlength; 
TF2 = bubbleproperties(:,7)-bubbleproperties(:,6)<minCSlength;
TF3 = bubbleproperties(:,11)-bubbleproperties(:,10)<minCSlength;
TF4 = bubbleproperties(:,5)<minbubbledia;
bubbleproperties(TF1 | TF2 | TF3 | TF4,:) = []; 

% the first column now gets frame number 
bubbleproperties(:,1) = framei;
bubblepropertiestotal = [bubblepropertiestotal; bubbleproperties];
% bubblepropertiestotal = [frame#, xmean, ymean, zmean, bubble-dia, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2]
end

bubblepropertiestotal(1,:) = [];            % remove first row of zeros 
TF = bubblepropertiestotal(:,9) > ycutoff2 | bubblepropertiestotal(:,8) < ycutoff1 ; 
bubblepropertiestotal(TF,:) = [];           % remove bubbles touching top and bottom of bed

% add frame1 to recover original frame numbers 
bubblepropertiestotal(:,1) = bubblepropertiestotal(:,1) + frame1 - 1; 
% modify coordinates for cylgeometry 
if cylgeometry==1 && cylcoord == 0
    bubblepropertiestotal(:,2) = bubblepropertiestotal(:,2)-R/2; % for this case, R = bed diameter 
    bubblepropertiestotal(:,4) = bubblepropertiestotal(:,4)-R/2; % for this case, R = bed diameter 
end

end





function [R,nr,H,ny,Z,nz,A] = func_readgeometry(cylcoord)

File = 'Geometry.xlsx';
R = xlsread(File,'Sheet1','C3');
nr = xlsread(File,'Sheet1','C4');
H = xlsread(File,'Sheet1','E3');
ny = xlsread(File,'Sheet1','E4');
Z = xlsread(File,'Sheet1','G3');
nz = xlsread(File,'Sheet1','G4');

% x direction 
[num,txt,raw] = xlsread(File,'Sheet1','C5');
if strcmp(txt,'no')
    drrange = strcat('C7:C',num2str(7+nr-1));
    dr = xlsread(File,'Sheet1',drrange);
else dr = (R/nr)*ones(nr,1);
end

% y direction 
[num,txt,raw] = xlsread(File,'Sheet1','E5');
if strcmp(txt,'no')
    dyrange = strcat('C7:C',num2str(7+ny-1));
    dy = xlsread(File,'Sheet1',drrange);
else dy = (H/ny)*ones(ny,1);
end

% z direction 
if cylcoord==1; Z = 2*pi; dz = (Z/nz)*ones(nz,1); 
else % cylcoord==0
    [num,txt,raw] = xlsread(File,'Sheet1','G5');
    if strcmp(txt,'no')
        dzrange = strcat('G7:G',num2str(7+nz-1));
        dz = xlsread(File,'Sheet1',drrange);
    else dz = (Z/nz)*ones(nz,1);
    end
end 
    
%--------------------------------------------------------
% setup coarsegrid to replicate simulation grid 

% setup coarse grid in x direction 
xtemp(1)=dr(1);         % to get wall centers 
for i=2:size(dr)
    xtemp(i)=0; for k=1:i; xtemp(i)=xtemp(i)+dr(k); end
end
xtemp = [0 xtemp];      % to get cell centers 
for i=1:size(dr); xcell(i)=0.5*(xtemp(i+1)+xtemp(i)); end   

% setup coarse grid in y direction 
ytemp(1)=dy(1);         
for i=2:size(dy)
    ytemp(i)=0; for k=1:i; ytemp(i)=ytemp(i)+dy(k); end
end
ytemp = [0 ytemp];
for i=1:size(dy); ycell(i)=0.5*(ytemp(i+1)+ytemp(i)); end

% setup coarse grid in z direction 
ztemp(1)=dz(1);         
for i=2:size(dz) 
    ztemp(i)=0; for k=1:i; ztemp(i)=ztemp(i)+dz(k); end
end
ztemp = [0 ztemp];
for i=1:size(dz); zcell(i)=0.5*(ztemp(i+1)+ztemp(i)); end

xgeom = zeros(nr*ny*nz,1);
ygeom = zeros(nr*ny*nz,1);
zgeom = zeros(nr*ny*nz,1);
cellgeom = linspace(1,nr*ny*nz,nr*ny*nz); 

% cell order must replicate order in which data is written
% assume r/x first, then y, and then z/theta 

for i=1:nr*ny*nz
    if mod(i,nr*ny)>0
        zgeom(i) = zcell(1+(i-mod(i,nr*ny))/(nr*ny));
        zgeomctr = 1+(i-mod(i,nr*ny))/(nr*ny);
    else
        zgeom(i) = zcell(i/(nr*ny));
        zgeomctr = i/(nr*ny); 
    end        
    index = i-(zgeomctr-1)*nr*ny;   
    if mod(index,nr)>0
        ygeom(i)= ycell(1+(index-mod(index,nr))/nr);
        ygeomctr = 1+(index-mod(index,nr))/nr;
    else
        ygeom(i) = ycell(index/nr);
        ygeomctr = index/nr;
    end    
    xgeom(i) = xcell(index-(ygeomctr-1)*nr);       
end

cellgeom = cellgeom'; 
A = [cellgeom, xgeom, ygeom, zgeom];
clear xtemp; clear ytemp; clear xcell; clear ycell; clear cellgeom; 

end





function [C] = func_bubbleproperties(B, deltax, deltay, deltaz)

% B = [x y z bubble#]

% Adding fictitious row in B for ease of indices later 
B = [B; 10000, 10000, 10000, 10000]; 
[s1,s2] = ismember(unique(B(:,4)),B(:,4));
C = [unique(B(:,4)), s2];           % bubble number and first cell 

% last row in C will be fictitious (due to fictitious last row in B) 
for i=1:size(C(:,1))-1
    C(i,3) = mean(B(C(i,2):C(i+1,2)-1,1));          % xmean 
    C(i,4) = mean(B(C(i,2):C(i+1,2)-1,2));          % ymean
    C(i,5) = mean(B(C(i,2):C(i+1,2)-1,3));          % zmean    
    C(i,6) = (C(i+1,2)-C(i,2))*deltax*deltay*deltaz;% Volume    
    C(i,7) = min(B(C(i,2):C(i+1,2)-1,1));           % xmin
    C(i,8) = max(B(C(i,2):C(i+1,2)-1,1));           % xmax
    C(i,9) = min(B(C(i,2):C(i+1,2)-1,2));           % ymin
    C(i,10) = max(B(C(i,2):C(i+1,2)-1,2));          % ymax
    C(i,11) = min(B(C(i,2):C(i+1,2)-1,3));          % zmin
    C(i,12) = max(B(C(i,2):C(i+1,2)-1,3));          % zmax
    C(i,13) =(C(i,10)-C(i,9))/(C(i,8)-C(i,7));      % AR1
    C(i,14) =(C(i,10)-C(i,9))/(C(i,12)-C(i,11));    % AR2
end

B(end,:) = []; C(end,:) = []; 
% sortrows based on ylocation 
C = sortrows(C,4); C(:,2) = []; 
% renumber bubbles 
C(:,1) = linspace(1,length(C(:,1)),length(C(:,1)))';
% C = [bubble#, xmean, ymean, zmean, Volume, xmin, xmax, ymin, ymax, zmin, zmax, AR1, AR2]

end




