function [spike_inds, quality_flags] = BGC_spiketest(floatID,cycle,profiledata,param,savedir,varargin)

% ************************************************************************
% BGC_spiketest.m
% ************************************************************************
%
% Function to identify spikes in profile data for BGC parameters.
% Identification is based on thresholds placed on a calculated "test value" which
% is the difference between the point in question within the profile and
% the median value of a running 5-pt filter.  Due to the nature of the
% filtering, the first and last two points within a profile are inherently
% not subject to the spiketest.
%
% Whenever spikes are identified, they are written to a textfile for record
% keeping to track performance of the algorithm.
%
% USE AS: [spike_inds, quality_flag] = BGC_spiketest('9666SOOCN',12,profiledata,'NO3')
%
% INPUTS:
%    floatID      :  float for which the test is being performed (to track identified spikes over time)
%    cycle        :  cycle for which the test is being performed (to track identified spikes over time)
%    profiledata  :  array of profile data for parameter of interst.  Must
%                    include 2 columns: [pressure; paramter data]
%    param        :  string identifier for parameter of interest, either 'NO3','PH',or 'O2'.
%    (fillvale)      : (optional input argument)  if you'd like to screen for
%                       fill values (ie if using 99999 instead of nan for
%                       missing values), enter fill value here.  They will be
%                       screened out prior to imposing the spike algorithm.
%                       Original spike indices are reassigned to function output.
%    (QFscreen)     : (optional input argument) logical array of
%                       size(profiledata, 1).  1 = data to use; 0 = data to omit (ie data that
%                       has already been identified as bad, say, from range checks.)
%    
%
% OUTPUTS: 
%   spike_inds    : row indices of values within 'profiledata' that have been
%                   identified as spikes
%   quality_flags : Argo quality flag that should be assigned to identified
%                   spike (vector the same size as 'spike_inds')
%
%
% USE AS:
%   [spike_inds, quality_flags] = BGC_spiketest('9099SOOCN',25,[PRESvec; PHvec],'PH','C:\Users\tmaurer\Documents\MATLAB\ARGO_PROCESSING\MFILES\FLOATS\',99999)
%   [spike_inds, quality_flags] = BGC_spiketest('9099SOOCN',25,[PRESvec; PHvec],'PH','C:\Users\tmaurer\Documents\MATLAB\ARGO_PROCESSING\MFILES\FLOATS\',nan,[QFvec])
%
%
% AUTHOR: Tanya Maurer
%         MBARI
%         tmaurer@mbari.org
%
% DATE: 05/10/2018
% UPDATES:
% NOTES: 
% ************************************************************************
%
% ************************************************************************
%--------------------------------------------------------------------------

% disp(['RUNNING ',param,' SPIKE TEST ALGORITHM FOR FLOAT ',floatID,',
% CYCLE ',num2str(cycle),'...']) %too much text output!

% CHECK INPUTS
if nargin < 5 || nargin > 7
    disp('WRONG NUMBER OF INPUT ARGUMENTS!')
    help BGC_spiketest
    return
end

if nargin == 6
    if length(varargin{1}) ~=1
        disp('Invalid input arguments.')
        help BGC_spiketest
        return
    end
    fillvalue = varargin{1};
end

if nargin == 7
    if (length(varargin{1}) ~=1) || (length(varargin{2}) ~= size(profiledata,1))
        disp('Invalid input arguments.')
        help BGC_spiketest
        return
    end
    QFscreen = varargin{2};
    fillvalue = varargin{1};
end

%--------------------------------------------------------------------------
% FIRST, IDENTIFY WHICH PARAMETER THE TEST IS ACTING ON 
% (Can add CHL later)
if strcmp(param,'NO3') == 1
    thresh1 = 5;
    thresh2 = 3;
elseif strcmp(param,'PH') == 1
    thresh1 = 0.04;
    thresh2 = 0.03;
elseif strcmp(param,'O2') == 1
    thresh1 = 40;
    thresh2 = 30;
end

%--------------------------------------------------------------------------
% REMOVE PRE-SCREENED BAD DATA, FILL VALUE AND/OR NANS FROM RECORD (I do it this way, as opposed
% to incrementally removing flagged data, fill, then nan, in order to maintain all
% original row indices in one index vector (I)

% IF USER INPUTS INCLUDE FILL VALUE, CONVERT TO NAN BEFORE PROCEEDING
if exist('fillvalue','var')==1
    [K,L] = find(profiledata==fillvalue);
    profiledata(K,L) = nan;
end

% IF USER INPUTS INCLUDE LOGICAL ARRAY TO IDENTIFY PRE-SCREENED BAD DATA,
% CONVERT TO NAN BEFORE PROCEEDING
if exist('QFscreen','var')==1
    profiledata(QFscreen,:) = nan;
end

% REMOVE NANS FROM RECORD PRIOR TO IMPOSING SPIKE ALGORITHM
if ~isempty(profiledata)
    [I,J] = find(~isnan(profiledata)); %find nans in either pressure or parameter data
    myindices = unique(I);
    nonan_data = profiledata(myindices,:); %remove both column indices for either column index with identified nan (although columns should always be similar)
    depth_data = nonan_data(:,1);
    param_data = nonan_data(:,2);

    %--------------------------------------------------------------------------
    % LOOP THROUGH PROFILE AND CALCULATE "TEST VALUE" FOR EACH POINT ALONG
    % PROFILE
    K = 1;
    TV = [];
    TV_P = [];
    if ~isempty(nonan_data) && size(nonan_data,1) >=5 %can only perform spiketest on profile with 5 or more points!
        for j = 3:size(nonan_data,1)-2  %because of 5pt filter, first and last 2 datapoints are not subject to the test
                    V0 = nonan_data(j-2,2);
                    V1 = nonan_data(j-1,2);
                    V2 = nonan_data(j,2);
                    V3 = nonan_data(j+1,2);
                    V4 = nonan_data(j+2,2);
                    tv = abs(V2-nanmedian([V0 V1 V2 V3 V4])); %test value algorithm
                    TV(K) = tv; %holds test values for entire profile 
                    TV_P(K) = nonan_data(j,1); %corresponding pressure values
                    K=K+1;
         end
    else
        disp(['PROFILE IS TOO SPARSE (OR EMPTY) FOR FLOAT ',floatID,', CYCLE ',num2str(cycle),'.  NO SPIKETEST WAS PERFORMED.'])
        spike_inds = [];
        quality_flags = [];
        return
    end

    %--------------------------------------------------------------------------
    % IMPOSE DEFINED THRESHOLD AND ASSIGN FLAGGING ACCORDINGLY 
    nTV = [nan nan TV nan nan]; %fill out the array (no test values for first and last two row indices)
    nTV_P = [nan nan TV_P nan nan]; %fill out the array
    TVarray = [myindices,nTV',nTV_P']; % [original row index; original column index; test value; corresponding pressure value]
    tmp_inds = find(TVarray(:,2) > thresh1); %get spike indices within no-nan array
    if ~isempty(tmp_inds)
        spike_inds1 = TVarray(tmp_inds,1); %back out to find index of initial array
        spike_inds = spike_inds1;
        quality_flags = repmat(4,length(spike_inds1),1); %assign Argo quality flag 4 (bad data) for identified spikes
        disp([num2str(length(spike_inds)),' ',param,' SPIKES WERE DETECTED FOR FLOAT ',floatID,' CYCLE ',num2str(cycle),'.'])
        % *-*-Simple plots for testing purposes.  Can comment out operationally.
    %     figure
    %     plot(profiledata(:,2),profiledata(:,1),'b.-')
    %     hold on
    %     plot(profiledata(spike_inds,2),profiledata(spike_inds,1),'ro')
    %     title([floatID,' cycle',num2str(cycle)])
    %     set(gca,'Ydir', 'Reverse')
        %*-*-
    else
        spike_inds = [];
        quality_flags = [];
    %     disp(['NO ',param,' SPIKES WERE DETECTED FOR FLOAT ',floatID,' CYCLE ',num2str(cycle),'.'])
    end
    %---
    %can add code here to implement a second threshold test and set quality
    %flags accordingly
    %---

    %--------------------------------------------------------------------------
    % NOW SOME RECORD KEEPING TO HELP TRACK AND ASSESS THE ALGORITHM
    % PERFORMANCE OVER TIME.  LOG SPIKES DETECTED AND FLAGGED FOR LATER REVIEW.
    spike_fname = [savedir,'MBARI_BGCArgo_float_spikes.txt'];
    if ~isempty(spike_inds)
        thedtg = datestr(now,31); %grab current time
        mycycs = repmat(cycle,length(spike_inds),1);
        MYD = [mycycs, profiledata(spike_inds,1), profiledata(spike_inds,2)];
        C = cell(length(spike_inds1),1);
        for jj = 1:size(MYD,1)
            C(jj,4) = {MYD(jj,1)};
            C(jj,5) = {MYD(jj,2)};
            C(jj,6) = {MYD(jj,3)};
        end
        C(:,1) = {thedtg};
        C(:,2) = {floatID};
        C(:,3) = {param};
        fid = fopen(spike_fname,'a');
        formatspec = '\n%s %s %s %d %d %6.2f %6.2f\n';
        [nrows, ncols] = size(C);
        for row = 1:nrows
            fprintf(fid,formatspec,C{row,:});
        end
        fclose(fid);
    end
else
    spike_inds = [];
    quality_flags = [];  
    disp(['HR or LR ',param,' PROFILE DATA IS EMPTY FOR',floatID,' CYCLE ',num2str(cycle),'.'])
end

% disp('SPIKE TEST COMPLETE.')
