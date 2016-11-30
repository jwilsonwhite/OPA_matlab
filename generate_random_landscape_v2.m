function [XX, Rad] = generate_random_landscape_v2(maxR,n,sill_points,sill_radii,radius_mean,radius_sd,Doplot)

% generate landscape with various types of weighted random settlement
% for Deignan project

% maxR = maximum dimension of the (square) study area
% n = # points
% sill_points gives the target autocorrelation scale of the point pattern data
%       - a value 0 < x <= 1 will produce a random distribution (the actual
%       value does not matter if it falls within that range)
%       - a value > 1 produces an autocorrelated pattern with that
%       autocorrelation scale
%       - a value < 0 produces an underdispersed pattern with that value
%       giving the minimum distance between points
%       - a value of 0 produces an error
% sill_radii acts in the same way as sill_points, but acts on the object
%       radii.  So values 0 < x <= 1 assign radii at random; values > 1
%       produce objects with radii that are spatially autocorrelated, and
%       values < 0 produce a pattern in which objects repel each other at a
%       scale equal to sill_radii/2 beyond the boundary of their object
%       edge. In other words there will be a minimum distance of sill_radii
%       between the edges of any object.
% radius_mean and radius_sd give the mean & sd of the objects.  Take care
%       with this; if values are too large then very few points will be
%       able to fit on the landscape.  If the SD is very large relative to
%       the mean then many objects will be constrained to have radius 0 and
%       the radii will have a truncated size distribution
% Doplot - logical; plot the object pattern or not?

% OUTPUTS:
% XX: matrix with columns X & Y giving point coordinates (assuming NW
% corner is 0,0)
% Rad: vector of radii of each point
%  


if sill_points == 0
    error('sill_points cannot equal 0')
end
if sill_radii == 0
    error('sill_radii cannot equal 0')
end
if radius_sd <= 0
    error('radius_sd cannot be ? 0')
end

if ~exist('Doplot','var')
    Doplot = false;
end


% 1) Simulate points.  Choose whether they will be random, autocorrelated,
% or underdispersed with respect to the point-pattern

% this will be a square of dimension maxR, with values distributed ~N(0,1)
% with spatial scale sill_points

% If X should be underdispersed, apply that correction:
if sill_points < 0
    X = simulate_autocorr_v2(round(maxR),0,1,1,false);
    X = Overdisperse(X,abs(sill_points),n);
else
    % otherwise, random or autocorrelated
    X = simulate_autocorr_v2(round(maxR),0,1,sill_points,false);
    X = Discretize(X,n); % convert into discrete points
end

% Simulated radii.  Eliminate overlapping points.  Choose whether they
% should be random, autocorrelated, or dispersed with respect to radii.
if sill_radii < 0
    R = simulate_autocorr_v2(round(maxR),radius_mean,radius_sd,1,false); % random distribution of radii, randomly distributed
    R = max(0,R); % constrain to be >= 0
    % Eliminate points that overlap given radius + sill_radii
    X = No_overlap(X,R + abs(sill_radii)/2); % use same code as below, but add in 1/2 of sill to each point
else % 
    R = simulate_autocorr_v2(round(maxR),radius_mean,radius_sd,sill_radii,false);
    R = max(0,R); % constrain to be >= 0
    % Eliminate points that overlap
    X = No_overlap(X,R);
end


% Convert into data table of X, Y points + Radii
[x,y] = find(X);
XX = [x,y];
Rad = R(X);

% Optional: plot points + radii
if Doplot
        
    figure
    hold on
    plot(x,y,'ko','markerfacecolor','k');
    
    % Plot circles
for i = 1:length(Rad)
    v(1:2) = XX(i,:);
    theta = linspace(0,2*pi,360);
    [x1,y1] = pol2cart(theta,Rad(i));
    x2 = v(1)+x1;
    y2 = v(2)+y1;
    plot(x2,y2,'k')
end
axis equal

    xlim([min(x),max(x)]);
    ylim([min(x),max(x)]);
end % end if Doplot



%----------------------------------------------------------------------
function Out = No_overlap(In,R)
% Eliminate overlapping points, with respect to radii

Oa = -0.05; % Proportional overlap allowed between two objects (in proportion of distance between point centers)
% Distance matrix for each point:
X = repmat(1:size(In,2),[size(In,1),1]);
Y = X'; % assumes square matrix
Q = dist([X(:)';Y(:)']); % distance between each
Q = Q + diag(nan(length(Q),1)); % eliminate 'self' distances

% Add in radii
R2 = R(:);
Rmat = repmat(R2,[1,size(In,2)^2]) + repmat(R2',[size(In,2)^2,1]); % all pairwise radius sums
QR = (Q-Rmat)./Q; % proportional difference between distance & summed radii.  If negative, points overlap


% Only take points that are in X
isX = logical(In(:));
QRtmp = QR(isX,isX);
OK = ~any(QRtmp(:)<Oa); % criterion: all QR differences are positive

while ~OK
    
    % Find the point with the most failures 
    Tmp = sum(QRtmp<Oa);
    Tmp = Tmp == max(Tmp); 
    isXI = find(isX);
    isXI = isXI(~Tmp); % remove that point
    isX = false(length(isX),1); % reset isX
    isX(isXI) = true; % restore isX, now without the offending point
    
    QRtmp = QR(isX,isX);
    OK = ~any(QRtmp(:)<Oa);
    
    if sum(isX) == 1 % in case we are down to one point
        OK = true;
    end
    
end % end while
% taking all the points, not just some of the offending points?
InC = false(size(In(:)));
InC(isX) = true;
Out = reshape(InC,size(In));
%----------------------------------------------------------------------

%----------------------------------------------------------------------
function Out = Overdisperse(In,Sill,n)
% Distance matrix for each point:
X = repmat(1:size(In,2),[size(In,1),1]);
Y = X'; % assumes square matrix
Q = dist([X(:)';Y(:)']); % distance between each
Q = Q + diag(nan(length(Q),1)); % eliminate 'self' distances

[H,I] = sort(In(:),'descend');
OKrows = 1:n;
Htmp = H(OKrows);
Itmp = I(OKrows);
ntmp = n+1;

Qtmp = Q(Itmp,Itmp); 

OK = nanmin(Qtmp(:))>Sill; % if the minimum distance of that set is < Sill
while ~OK
   
    % find column of Qtmp with smallest mean distances
    Smallrows = find(nanmin(Qtmp)<Sill,1);
    OKtmp = true(size(OKrows));
    OKtmp(Smallrows)=false;
    OKrows = OKrows(OKtmp);
    OKrows(end+1) = ntmp; % add additional new point to fix the one removed
    
    Htmp = H(OKrows);
    Itmp = I(OKrows);
    Qtmp = Q(Itmp,Itmp);
    
    OK = nanmin(Qtmp(:))>Sill;
    ntmp = ntmp+1;
    if ntmp > length(I)
        OK = true;
    end
end

Out = zeros(size(H));
Out(Itmp) = 1;
Out = reshape(Out,size(In));
%----------------------------------------------------------------------
   
%----------------------------------------------------------------------
function Out = Discretize(In,n)
% Converts a random normal landscape into a discrete set of n points
H = flipud(sortrows(In(:))); % sort in descending order
Hmin = H(n); % pick the first n values
Out = In>=Hmin; % Out is a logical with n TRUEs
Out = reshape(Out,size(In));
%----------------------------------------------------------------------

