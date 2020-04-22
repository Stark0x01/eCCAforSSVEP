% -------------------------------------------------------------------------
% Main for Extended Canonical Correlation Analysis[1]
%
% Dataset (Sx.mat):
%   A 40-target SSVEP dataset recorded from a single subject. The stimuli
%   were generated by the j oint frequency-phase modulation (JFPM)
%     - Stimulus frequencies    : 8.0 - 15.8 Hz with an interval of 0.2 Hz
%     - Stimulus phases         : 0pi, 0.5pi, 1.0pi, and 1.5pi
%     - # of channels           : Oz
%     - # of recording blocks   : 6
%     - Data length of epochs   : 1.5 [seconds]
%     - Sampling rate           : 250 [Hz]
%     - Data format             : # channels, # points, # targets, # blocks
% 
% See also:
%   CCA.m
%
% Reference:
%   [1] Chen X, Wang Y, Nakanishi M, Gao X, Jung TP, Gao S (2015)
%       High-speed spelling with a noninvasive brain-computer interface. 
%       PNAS 112:E6058-6067.
% -------------------------------------------------------------------------

clear all
close all

load ('Freq_Phase.mat')
load('subject2.mat')
eeg = subject2;
[N_channel, N_point, N_target, N_block] = size(eeg);
% sample rate
fs = 250;
t=1/fs:1/fs:N_point/fs;   
%% ------------classification-------------
tic
% LOO cross-validation
for loocv_i = 1:N_block
     Testdata = eeg(:, :, :, loocv_i);
     Traindata = eeg;
     Traindata(:, :, :, loocv_i) = [];
     % number of harmonics
     N_harmonic = 2;
    for targ_i = 1:N_target
        % Template
        Template(:, :, targ_i) = mean(squeeze(Traindata(:,:,targ_i,:)),3);
        % Reference
        Y=[];
        for har_i=1:N_harmonic
             Y=cat(1,Y,cat(1, sin(2*pi*freqs(targ_i)*har_i*t), ...
                 cos(2*pi*freqs(targ_i)*har_i*t)));        
        end
        Reference(:, :, targ_i) = Y;
    end
    
    % labels assignment according to testdata
    truelabels=freqs;
    
    N_testTrial=size(Testdata, 3);
    for trial_i=1:N_testTrial
        Allcoefficience = [];
        for targ_j=1:length(freqs)             
            % 老1 (Filter: Test data & Reference)
            [wn1, wn2]= CCA(Testdata(:,:,trial_i), Reference(:, :, targ_j));
            weighted_train = wn2'*Reference(:,:,targ_j);
            weighted_test = wn1'*Testdata(:,:,trial_i);
            coefficienceMatrix = corrcoef(weighted_test,weighted_train);
            coefficience(1) = abs(coefficienceMatrix(1,2));
            % 老2 (Filter: Test data & Template)
            [wn, ~] = CCA(Testdata(:,:,trial_i), Template(:, :, targ_j));
            weighted_train = wn'*Template(:,:,targ_j);
            weighted_test = wn'*Testdata(:,:,trial_i);
            coefficienceMatrix = corrcoef(weighted_test,weighted_train);
            coefficience(2) = coefficienceMatrix(1,2);
             % 老3 (Filter: Test data & Reference)
            [wn, ~] = CCA(Testdata(:,:,trial_i), Reference(:, :, targ_j));
            weighted_train = wn'*Template(:,:,targ_j);
            weighted_test = wn'*Testdata(:,:,trial_i);
            coefficienceMatrix = corrcoef(weighted_test,weighted_train);
            coefficience(3) = coefficienceMatrix(1,2);
            % 老4 (Filter: Template & Reference)
            [wn, ~] = CCA(Template(:, :, targ_j), Reference(:, :, targ_j));
            weighted_train = wn'*Template(:,:,targ_j);
            weighted_test = wn'*Testdata(:,:,trial_i);
            coefficienceMatrix = corrcoef(weighted_test,weighted_train);
            coefficience(4) = coefficienceMatrix(1,2);
%              % 老5 (Filter: Test data & Template)
%             [wn1, wn2] = CCA(Testdata(:,:,trial_i), Template(:, :, targ_j));
%             weighted_train = wn2'*Template(:,:,targ_j);
%             weighted_test = wn1'*Template(:,:,targ_j);
%             coefficienceMatrix = corrcoef(weighted_test,weighted_train);
%             coefficience(5) = coefficienceMatrix(1,2);
            % compute correlation values
            Allcoefficience(targ_j) = sum(sign(coefficience).*coefficience.^2);
        end % end targ_i
            % target detection
            [~, index] = max(Allcoefficience);
            outputlabels(trial_i) = freqs(index);
            
    end % end trial_i
    trueNum = sum((outputlabels-truelabels)==0);
    acc(loocv_i) = trueNum/length(truelabels);
    fprintf('The %d-th CV accuracy is: %.4f, samples: %d/%d\n',loocv_i,...
        acc(loocv_i),trueNum, N_testTrial)
end % end looCv_i
t=toc;
% data visualization
fprintf('\n-----------------------------------------\n')
disp(['total time: ',num2str(t),' s']);
fprintf('6-fold CV average accuracy is: %.4f\n',mean(acc))
