function [K, R, Kint] = Ripleys_Kc_intervals(Data,Type,Eps,R)

% Calculate Ripleys Kc (object pattern analysis)
% Along with confidence intervals

n = size(Data,1); % number of data points
N = 1e4; % number of randomizations (should be at least 1e3)

% Center the data
Data(:,1:2) = Data(:,1:2) - repmat(mean(Data(:,1:2)),[length(Data(:,1)),1]);
maxR = max((Data(:,1).^2 + Data(:,2).^2).^0.5)+max(Data(:,3)); % max radius of center points + largest radius
maxR = round(maxR); % has to be an integer for the random landscape scenario

if ~exist('Type','var')
    error('Type must be circle or rect')
end

if ~exist('Eps','var')
    Eps = 0.01;
end

if ~exist('R','var')
% Scales of influence
R = linspace(quantile(Data(:,3),0.75),(maxR/2),10); 
end

% The actual values for this dataset
    K = Ripleys_Kc(Data,Type,Eps,R);
 
% Intervals around them:  
XX = nan(n,2,N);
Rad = nan(n,N);
Kt = nan(N,length(R));
    
   % Generate random landscapes:
   for i = 1:N
   XX(:,1,i) = rand(n,1); % random X coordinates
   XX(:,2,i) = rand(n,1); % random Y coordinates
   Rad(:,i) = Data(randi(length(Data),1),3);
   end
   
   XX = XX.*maxR; % scale to size of original plot
   
   for i = 1:N
       Kt(i,:) = Ripleys_Kc([XX(:,:,i),Rad(:,i)],Type,Eps,R);
   end
   
   Kint = quantile(Kt,[0.025, 0.975],1);
