function [K,R] = Ripleys_Kc(Data,Type,Eps,R)

% perform Ripleys Kc as described by Protazio (2007)

% Assume Data is X, Y, Radius
% Type has to be 'rect' or 'circle'  (If rectangle we assume that it is
% aligned (approximately) with the axes of the Cartesian plane
% You can either specify R or let it calculate automatically

% Other inputs:
% Eps is the width of the ring used in the ring pattern calculation.
%     It works best if it is quite small relative to the radii of the data
%     and the sample radii R.  It will be automatically calculated if not
%     provided.
% Specifying Eps = 0 will coerce the code to perform regular point-pattern
% Ripley's K calculation, ignoring the radii
% 
% R is the range of sample radii to be considered.  Should not exceed 50%
% of the sample plot scale.  Will be calculated automatically if not `
% provided.

% Output: 
% K is the value of Kc (or K) calculated as each spatial scale R
% Values of K > 1 indicate clustering; values < 1 indicate regularity
% Note that the smallest value of R is often much greater than 1 just
% because it lies interior to some of the data objects.


if ~any(strcmp(Type,{'circle','rect'}))
    error('Type MUST BE EITHER circle OR rectangle')
end

% If the data are all points, must set Eps = 0 so it calculates
% point-pattern Ripley's K
if max(Data(:,3))==0 
    Eps = 0;
end
  
n = size(Data,1);

% Center the data
Data(:,1:2) = Data(:,1:2) - repmat(mean(Data(:,1:2)),[length(Data(:,1)),1]);

switch Type
    case 'circle'
    maxR = max((Data(:,1).^2 + Data(:,2).^2).^0.5)+max(Data(:,3)); % max radius of center points + largest radius (note: scalar)
    case 'rect'
    maxR = max(abs(Data(:,1:2))) + max(Data(:,3)); % max X or Y + largest radius (note: 1x2 vector) 
    % Vertex points of the rectangle
    Vertices = [-maxR(1), maxR(2);...
                maxR(1), maxR(2);...
                maxR(1), -maxR(2);...
                -maxR(1), -maxR(2)];
        
end % end switch Type

if ~exist('R','var')
% Scales of influence
R = linspace(quantile(Data(:,3),0.75),max((maxR/2)),20); 
end

if ~exist('Eps','var')
% Ring width
Eps = 0.001*min(R);
end

% If calculating point-pattern Ripley's K, ensure that we ignore radii
if Eps == 0
    Data(:,3) = 0;
end


% Calculate mu: proportional area occupied by the objects in Data
switch Type
    case 'circle'
    % Circle areas:
    Ca = sum(Data(:,3).^2); % area of each circle
    totalA = pi*maxR.^2; % total area of the study plot (used later for point-pattern K)
    Mu = Ca./(maxR.^2);     % area of the study area (leave out pi's bc it cancels out)
    case 'rect'
    Ca = pi*sum(Data(:,3).^2); % area of each circle
    totalA = prod(2*maxR);
    Mu = Ca ./ prod(2*maxR);   % area of the total data rectangle
end % end switch Type

A = nan(n,n,length(R));
Ai = nan(n,length(R));
wi = nan(n,1);
wo = nan(n,1);

for r = 1:length(R) % Loop over each scale of clustering
for i = 1:n         % Loop over each data point
            
    %---------------------------------------------------------------------
    % Calculate weighting factors wo & wi
    % (proportion of the study area that is within the circle of influence
    
    switch Type
        case 'circle'
    % If circle lies entirely within study area
    OL1 = is_overlap(Data(i,:),[0,0,maxR],R(r)+Eps/2);
    OL2 = is_overlap(Data(i,:),[0,0,maxR],R(r)-Eps/2);
    
    if OL1
    wo(i) = 1;
    else
    wo(i) = Afxn(Data(i,:),[0,0,maxR],R(r)+Eps/2)./(pi*(R(r)+Eps/2)^2);
    end
    
    if OL2
    wi(i) = 1;
    else
    wi(i) = Afxn(Data(i,:),[0,0,maxR],R(r)-Eps/2)./(pi*(R(r)-Eps/2)^2);
    end    
        case 'rect'
    % If circle lies entirely within rectangular study area
    OL1 = is_overlap_rect(Data(i,:),Vertices,R(r)+Eps/2);
    OL2 = is_overlap_rect(Data(i,:),Vertices,R(r)-Eps/2);    
        
    if OL1
    wo(i) = 1;
    else
    wo(i) = Afxn_rect(Data(i,:),Vertices,R(r)+Eps/2)./(pi*(R(r)+Eps/2)^2);
    end
    
    if OL2
    wi(i) = 1;
    else
    wi(i) = Afxn_rect(Data(i,:),Vertices,R(r)-Eps/2)./(pi*(R(r)-Eps/2)^2);
    end   
    
    end % end switch Type
    % End calculation of wo & wi
    %---------------------------------------------------------------------
    
    for j = 1:n % loop over every potential neighbor point:
        
        if j == i % if not the same point
           A(i,j,r) = 0;
        else
        
        d  = ((Data(i,1)-Data(j,1)).^2 + (Data(i,2)-Data(j,2)).^2).^0.5;
        rj = Data(j,3);
        
        % If just doing regular Ripley's K (not concentric rings):
        if Eps == 0
            
            if R(r) > d;
                A(i,j,r) = wo(i)*1/n; % if the two points are within the radius.  Scale by proportional overlap, wo (wo == wi in this case)
            else
                A(i,j,r) = 0;
            end
            
        else % otherwise do Ripley's Kc:
      
        % If i is totally inside j:
        if rj >= d + R(r) + Eps/2 
            A(i,j,r) = pi*( (R(r)+Eps/2).^2 - (R(r)-Eps/2).^2 );% total area of i
            
            
        elseif rj < d + R(r) + Eps/2 && rj > d + R(r) - Eps/2; % only the inner one is totally enveloped
            Outer =  Afxn(Data(i,:),Data(j,:),R(r)+Eps/2);
            Inner = pi*(R(r)-Eps/2).^2; % total area of i
            A(i,j,r) = wo(i)*Outer-wi(i)*Inner;
            
            
        else % if there is some opportunity for intersection in both rings:
            
        
        if R(r)+Eps/2 > d - rj %%% if they actually intersect - otherwise don't calculate
         
           if  R(r) + Eps/2 > rj + d % complete envelopment:
               Outer = pi.*rj^2;
           else
           Outer =  Afxn(Data(i,:),Data(j,:),R(r)+Eps/2);
           end
           
           if R(r) - Eps/2 > ((Data(i,1)-Data(j,1))^2 + (Data(i,2)-Data(j,2))^2)^0.5 - Data(j,3) % if the inner ring also intersects
           
               if  R(r) - Eps/2 > rj + d % complete envelopment:
               Inner = pi.*rj^2;
               else
               Inner = Afxn(Data(i,:),Data(j,:),R(r)-Eps/2);
               end
               
           else
           Inner = 0;
           end
            
        A(i,j,r) = wo(i)*Outer - wi(i)*Inner;
        
        else
        A(i,j,r) = 0;
        end
        end
        
        end % end if Eps == 0
        end % end if j == i
    end % end j
    
    Ai(i,r) = sum(A(i,:,r)); % sum of all j overlaps
    
    
end % end i
end % end r

if Eps > 0
K = sum(Ai,1)./(2.*pi.*R.*Eps)./n./Mu; 
% normalized ring object function, as in Protazio Eq. 4.10
else
K = sum(Ai,1)./(pi.*R.^2)./(n./totalA);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = Afxn(Di,Dj,r)

% radii
rj = Dj(3); % radius of target circle

% distance
d = ( (Di(1)-Dj(1)).^2 + (Di(2)-Dj(2)).^2).^0.5;

if d == 0
   % error('Two objects cannot have the same center point')
    
    A = pi*( min(r,rj))^2;
    
else

if d - rj > r
    warning('Objects do not intersect, overlap is a complex number')
end


% distance to radical line
x = abs((d^2 - rj^2 + r^2)/(2*d)); 

% if the radical line is between the two circle centers

if x <= d 

    % Proceed with the formula proposed by Protazio:
    A = do_Aij(d,rj,r);

else
    
    % Define a new circle, cX.
    
    if r > rj
    % cX has radius r and distance from c2 of 2x - d
   
    % Delete this stuff:
    % Th = atan2((Dj(2)-Di(2)),(Dj(1)-Di(1)));  % Find angle from circle 1 to circle 2
   % [X1 Y1] = pol2cart(Th,2*x);
    
    % Area is total c2 - overlap between cX & c2 + overlap between c1 & cX
    c2_area = pi*rj^2;
    cXc2_overlap = do_Aij(2*x-d,r,rj);
    c1cX_overlap = do_Aij(2*x,r,r);
    
    
    A = c2_area - cXc2_overlap + c1cX_overlap;
    
    if ~isreal(A)
        keyboard
    end
    
    else  % if r < rj, then the intersection happens behind the center of r
        
    c1_area = pi*r^2;
    cXc1_overlap = do_Aij(2*x-d,rj,r);
    c2cX_overlap = do_Aij(2*x,rj,rj);
    
    A = c1_area - cXc1_overlap + c2cX_overlap;
    
    if ~isreal(A)
        keyboard
    end
    
    end
    
end % end if x <= d
end % end if d == 0


% END Afxn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = Afxn_rect(Di,Vertices,r) 
% overlap between circle Di of radius r and rectangle defined by Vertices
% (NW, NE, SE, SW)

% Based on solution at http://stackoverflow.com/questions/622287/area-of-intersection-between-circle-and-rectangle

% Start with full circle area; subtract areas in each of the four external
% half-planes
Atotal = pi.*r.^2;

% signed (negative = outside) distances to each of the rectangle edges
d(1) = (Di(1) - Vertices(1,1)); % distance to W side
d(2) = (Vertices(1,2) - Di(2)); % distance to N side
d(3) = (Vertices(3,1) - Di(1)); % distance to E side
d(4) = (Di(2) - Vertices(3,2)); % distance to S side

% should have already tested for circle being interior
if any(d) < 0
    warning('Circle vertex is outside of study area')
end

% Calculate overlap with each edge of the rectangle:
A_edge = r .* (r .* acos( d./r ) - d.*sin(acos( d./r)) );
A_edge(d > r) = 0;

% Special case: overlaps a vertex
d2 = d; d2(5) = d(1); % for use in the following loop
A_c = zeros(4,1);
Corners = {'NW','NE','SE','SW'};
for i = 1:4 % rotate through each pair of corners
   if sqrt( (Di(1)-Vertices(i,1))^2 + (Di(2)-Vertices(i,2))^2) < r && d2(i) < r && d2(i+1) < r  % it intersects both sides & crosses the vertex
       A_c(i) = A_rect_corner(Di,Vertices(i,:),r,Corners{i});
   end
end
    
% Area is total area, minus overlaps, plus corner overlap (to avoid
% double-counting)

A = Atotal - sum(A_edge) + sum(A_c);
% END Afxn_rect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = A_rect_corner(Di,V,r,Corner)
% Overlap between circle (center Di, radius r) and quadrant defined by
% rectangle vertex V

% 1) get intersection points with each of the two lines:
IntX(1) = (r.^2 - (V(2) - Di(2)).^2).^0.5 + Di(1);
IntX(2) = -(r.^2 - (V(2) - Di(2)).^2).^0.5 + Di(1);

IntY(1) = (r.^2 - (V(1) - Di(1)).^2).^0.5 + Di(2);
IntY(2) = -(r.^2 - (V(1) - Di(1)).^2).^0.5 + Di(2);

% 2) choose intersection points that are in the correct quadrant
switch Corner
    case 'NW'
        Int(1,:) = [V(1), max(IntY)]; % northmost along Y line
        Int(2,:) = [min(IntX), V(2)]; % westmost along X line
    case 'NE'
        Int(1,:) = [V(1), max(IntY)]; % northmost along Y line
        Int(2,:) = [max(IntX), V(2)]; % eastmost along X line
    case 'SE'
        Int(1,:) = [V(1), min(IntY)]; % southmost along Y line
        Int(2,:) = [max(IntX), V(2)]; % eastmost along X line
    case 'SW'
        Int(1,:) = [V(1), min(IntY)]; % southmost along Y line
        Int(2,:) = [min(IntX), V(2)]; % westmost along X line
end

% 3) Find the length of the chord connecting the two intersection points
Z = sqrt( (Int(1,1) - Int(2,1))^2 + (Int(1,2) - Int(2,2))^2);

% 4) Find the angle of the segment
Th = asin( 0.5*Z/r);

% 5) Area of the segment
A_seg = r.^2 .* Th;

% 6) Area of the large isoceles triangle
A_tri1 = r^2 * cos(Th) * sin(Th);

% 7) Area of the smaller triangle
A_tri2 = 0.5 * ( abs(diff(Int(:,1))) * abs(diff(Int(:,2))) );

A = A_seg - A_tri1 + A_tri2;
% END A_rect_corner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Aij = do_Aij(d,rj,r)

if d == 0
   if r > rj
       Aij = pi*rj^2;
   else
       Aij = pi*r^2;
   end
   
else

% Formula proposed by Protazio (and found correctly on
    % Wolfram Mathworld) for an intersection forming a lens
Aij = (rj.^2.*acos( (d.^2 + rj.^2 - r.^2) ./ (2.*d.*rj) ) + ...
    r.^2.*acos( (d.^2 + r.^2 - rj.^2) ./ (2.*d.*r) ) ) - ...
    (0.5 .* ( (-d + rj + r).*(d + rj - r).*(d - rj + r).*(d + rj + r) ).^0.5);
% Note typo in Protazio's (2007) Eq. 4.3.  This is the correct version
end
% END do_Aij
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OL = is_overlap(Di,Dj,r)

% Is circle at i with radius r completely inside circle j?
d = ( (Di(1)-Dj(1)).^2 + (Di(2)-Dj(2)).^2).^0.5;
rj = Dj(3);

OL = d + r < rj;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OL = is_overlap_rect(Di,Vertices,r)

% Is circle at i with radius r completely inside rectangle with Vertices?
% Vertices listed as NW, NE, SE, SW

OL = abs(Di(1)) + r <= abs(Vertices(1,1)) && ...
     abs(Di(2)) + r <= abs(Vertices(1,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


