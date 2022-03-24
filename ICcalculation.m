%% calculate IC - simple VS partial correlation 
function [IC IC_lag] = ICcalculation(allsignal)

% INPUT: allsignal contain fMRI timecourses in a 1*N cell, 
%        where N is # of subjects.
%        Each cell element contains a i*j matrix in double,
%        where i is fMRI time points/volumes and j is # of ROIs, which has to be an
%        even number.
%        The arrangement of ROIs in timecourses has to follow the order of LEFT
%        ROIs first and then RIGHT ROIs.
% OUTPUT: IC is a 1*N vector with interhemispheric connectivity values
%         IC_lag is a 1*N vector with IC values that fMRI timecourses have
%         been shifted by 1 time point. (this measure is optional)
%
% ref: "Both Stationary and Dynamic Functional Interhemispheric Connectivity Are 
% Strongly Associated with Performance on Cognitive Tests in Multiple Sclerosis." Frontiers in Neurology, 2020
% Sue-Jin Lin 2020 April UBC/McGill University


% chack data 
nROIs = size(allsignal{1,1},2);
d= mod(nROIs,2);
if d~=0
    sprintf('error: the number of ROIs has to be even')
end 
n2 = nROIs/2;

% calculate IC - no lag 
for ij = 1:n2
    LH = setxor(1:n2,ij); RH = setxor(n2+1:nROIs,n2+ij); OtherROIs = setxor(1:nROIs,[LH,RH]);
    % for lag-0
    for k = 1:length(allsignal)
        signals = allsignal{k}';
        X = [signals(LH,:)']; % Left timecourse
        Y = [signals(RH,:)']; % Right timecourse
        OL = [ signals(OtherROIs,:)']; % timecourse in other ROIs
        r_homo = corrcoef([X,Y]); % full correlation matrix
        R_homo = r_homo(1:size(X,2),(size(X,2)+1):end); % matrix for IC 
        ICind = atanh(partialcorr(X,Y,OL)) - atanh(R_homo); % fisherZ add in 2020April
        IC(k, ij) = sum(abs(diag(ICind)));
    end
end

% calculate IC - lag 1 (optional)
for ij = 1:n2
    LH = setxor(1:n2,ij); RH = setxor(n2+1:nROIs,n2+ij); OtherROIs = setxor(1:nROIs,[LH,RH]);
    % for lag-1
    for k = 1:length(allsignal)
        signals = allsignal{k}';
        X = signals(LH,:)';
        Y = signals(RH,:)';
        OL = [circshift(signals(OtherROIs,:)',1)]; % shift timecourse by 1 time point
        r_homo_lag = corrcoef([X,Y]); % full correlation matrix
        R_homo_lag = r_homo_lag(1:size(X,2),(size(X,2)+1):end); % matrix for IC 
        ICind_lag = atanh(partialcorr(X,Y,OL)) - atanh(R_homo_lag); % fisherZ add in 2020April
        IC_lag(k, ij) = sum(abs(diag(ICind_lag)));
    end
end

