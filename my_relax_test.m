%%  ENF harmonic estimation algorithms with harmonic enhancement
%   This is the program comparing multi-tone harmonic model based ENF 
%   estimation algorithms, including a RELAX-based algorithm. The RELAX 
%   algorithm was proposed by Liu and Li in Ap. 1998.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close; clc;
%% parameter setting and initialization
FS                       = 800; % constant sampling frequency
HARMONIC_INDEX           = [2,3,4,5,6,7]; % constant value for ENF harmonic processing
fc                       = 50*HARMONIC_INDEX; % nominal frequencies at each harmonic
bound                    = 0.1*HARMONIC_INDEX; % tolerable IF deviations at each harmonic
filter_length            = 256;
[BPF_coeffs, coeffs_2nd] = func_BPF(filter_length);
index                    = dir('/media/blue/ckorgial/ENF-WHU-Dataset-master/ENF-WHU-Dataset/H1/*.wav');
index_ref                = dir('/media/blue/ckorgial/ENF-WHU-Dataset-master/ENF-WHU-Dataset/H1_ref/*.wav');

load('selected_harmonic_index.mat')
load('selected_harmonic_index0.mat')

bound_MLE                = selected_harmonic_index; % cell array - output of the GHSA
bound_MLE                = cellfun(@(x) x*0.1, bound_MLE, 'UniformOutput', false); % tolerable IF deviations at each harmonic

bound_WMLE               = selected_harmonic_index0; % cell array - output of the GHSA
bound_WMLE               = cellfun(@(x) x*0.1, bound_WMLE, 'UniformOutput', false); % tolerable IF deviations at each harmonic

f_ref                    = cell(1,length(index));

f_single_relax          = cell(1,length(index));
f_E_single_relax        = cell(1,length(index));
f_MLE_relax             = cell(1,length(index));
f_WMLE_relax            = cell(1,length(index));
f_E_MLE_relax           = cell(1,length(index));
f_E_WMLE_relax          = cell(1,length(index));
f_S_MLE_relax           = cell(1,length(index));
f_S_WMLE_relax          = cell(1,length(index));
f_P_MLE_relax           = cell(1,length(index));
f_P_WMLE_relax          = cell(1,length(index));

MSE_single_relax        = zeros(1,length(index));
MSE_E_single_relax      = zeros(1,length(index));
MSE_MLE_relax           = zeros(1,length(index));
MSE_WMLE_relax          = zeros(1,length(index));
MSE_E_MLE_relax         = zeros(1,length(index));
MSE_E_WMLE_relax        = zeros(1,length(index));
MSE_S_MLE_relax         = zeros(1,length(index));
MSE_S_WMLE_relax        = zeros(1,length(index));
MSE_P_MLE_relax         = zeros(1,length(index));
MSE_P_WMLE_relax        = zeros(1,length(index));
 
selected_harmonic_index0 = cell(1,length(index));
selected_harmonic_index  = cell(1,length(index));

tic;
for i =1:length(index)
    disp(['i=',num2str(i)]); 
    %% import audio recording and reference
    [audio, fs_audio] = audioread(strcat('/media/blue/ckorgial/ENF-WHU-Dataset-master/ENF-WHU-Dataset/H1/',index(i).name));
    [ref, fs_ref]     = audioread(strcat('/media/blue/ckorgial/ENF-WHU-Dataset-master/ENF-WHU-Dataset/H1_ref/',index_ref(i).name));
    ref               = ref';
    audio             = audio(:,1);
    audio             = audio';
    raw_wave          = resample(audio, FS, fs_audio);
    N                 = length(raw_wave);
    %% import the enhanced data
    input_E_single    = audioread(strcat('/media/blue/ckorgial/ENF-WHU-Dataset-master/ENF-WHU-Dataset/H1_enhanced_single/',index(i).name));
    input_E_single    = input_E_single(:,1);
    input_E_single    = input_E_single';

    input_E_multi     = audioread(strcat('/media/blue/ckorgial/ENF-WHU-Dataset-master/ENF-WHU-Dataset/H1_enhanced_multi/',index(i).name));
    input_E_multi     = input_E_multi(:,1);
    input_E_multi     = input_E_multi';
    %% bandpass filtering
    input             = filtfilt(BPF_coeffs,1,raw_wave);
    %% harmonic enhancement
    N_ite             = 2; % enhancement iterations
    h_rfa             = 3000; % window length of RFA
    initial_guess     = fc(1)*ones(1,N); % initial IF of RFA is fixed to 100 Hz
    TS                = 1; % constant time-step for RFA, 1 second
    window_dur_rfa    = 8; % enhancement window size 8 seconds
    FFT_res_rfa       = 200; % FFT resolution for RFA 1/FFT_res_rfa Hz
    refined_guess     = initial_guess;
    % % step 1: single-tone enhancement
    for k = 1:N_ite
        [input_denoised_single,~,refined_guess] = func_RFA(input,h_rfa,FS,TS,...
            refined_guess,fc(1),bound(1),window_dur_rfa,FFT_res_rfa);
    end
    % % step 2: use refined_guess to construct initial guesses for harmonics
    initial_guesses   = kron(HARMONIC_INDEX(2:end)'/2,refined_guess);
    refined_guesses   = initial_guesses;
    % % step 3: harmonic enhancement
    for k=1:N_ite
        [input_denoised_multi,~,refined_guesses] = func_RFA_multi(input,h_rfa,FS,TS,...
            refined_guesses,fc(2:end),bound(2:end),window_dur_rfa,FFT_res_rfa);
    end
    input_E_single    = input_denoised_single;
    input_E_multi     = sum([input_denoised_single;input_denoised_multi],1);
    %% ENF estimators 
    % set up parameters for frame-based processing
    window_dur        = 16; % duration of overlapping frame in second
    step_size_dur     = 1; % frame step-size usually 1 second
    FFT_res_factor    = 2000; % FFT resolution = 1/FFT_res_factor Hz
    f_ref{i}        = func_STFT_single_tone_relax(ref,fs_ref,window_dur,step_size_dur,50,0.1,FFT_res_factor)*2; % reference
    %% 1. single-tone estimation - RELAX
    f_single_relax{i}  = func_STFT_single_tone_relax(input,FS,window_dur,step_size_dur,fc(1),bound(1),FFT_res_factor);
    %% 2. single-tone enhancement - RELAX 
    f_E_single_relax{i}= func_STFT_single_tone_relax(input_E_single,FS,window_dur,step_size_dur,fc(1),bound(1),FFT_res_factor);
    %% 3. search within sum of harmonic components - RELAX
    f_MLE_relax{i}     = func_STFT_multi_tone_search_relax(input,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);
    %% 4. search within weighted sum of harmonic components - RELAX
    f_WMLE_relax{i}         = func_STFT_multi_tone_search_weighted_relax(input,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);
    %% 5. only enhancement MLE - RELAX
    f_E_MLE_relax{i}   = func_STFT_multi_tone_search_relax(input_E_multi,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);
    %% 6. only enhancement WMLE - RELAX
    f_E_WMLE_relax{i}  = func_STFT_multi_tone_search_weighted_relax(input_E_multi,FS,window_dur,step_size_dur,fc,bound,2*FFT_res_factor);
    %% 7. only harmonic selection - RELAX
    f_S_MLE_relax{i}   = func_STFT_multi_tone_search_relax(input,FS,window_dur,step_size_dur,fc,bound_MLE{i},2*FFT_res_factor);
    f_S_WMLE_relax{i}  = func_STFT_multi_tone_search_weighted_relax(input,FS,window_dur,step_size_dur,fc,bound_WMLE{i},2*FFT_res_factor);
    %% 8. both harmonic selection and enhancement - RELAX
    f_P_MLE_relax{i}    = func_STFT_multi_tone_search_relax(input_E_multi,FS,window_dur,step_size_dur,fc,bound_MLE{i},2*FFT_res_factor);
    f_P_WMLE_relax{i}   = func_STFT_multi_tone_search_weighted_relax(input_E_multi,FS,window_dur,step_size_dur,fc,bound_WMLE{i},2*FFT_res_factor);

    MSE_single_relax(i)   = 1/length(f_ref{i})*norm(f_single_relax{i}-f_ref{i}).^2;
    MSE_E_single_relax(i) = 1/length(f_ref{i})*norm(f_E_single_relax{i}-f_ref{i}).^2;
    MSE_MLE_relax(i)      = 1/length(f_ref{i})*norm(f_MLE_relax{i}-f_ref{i}).^2;
    MSE_WMLE_relax(i)     = 1/length(f_ref{i})*norm(f_WMLE_relax{i}-f_ref{i}).^2;
    MSE_E_MLE_relax(i)    = 1/length(f_ref{i})*norm(f_E_MLE_relax{i}-f_ref{i}).^2;
    MSE_E_WMLE_relax(i)   = 1/length(f_ref{i})*norm(f_E_WMLE_relax{i}-f_ref{i}).^2;
    MSE_S_MLE_relax(i)    = 1/length(f_ref{i})*norm(f_S_MLE_relax{i}-f_ref{i}).^2;
    MSE_S_WMLE_relax(i)   = 1/length(f_ref{i})*norm(f_S_WMLE_relax{i}-f_ref{i}).^2;
    MSE_P_MLE_relax(i)     = 1/length(f_ref{i})*norm(f_P_MLE_relax{i}-f_ref{i}).^2;
    MSE_P_WMLE_relax(i)    = 1/length(f_ref{i})*norm(f_P_WMLE_relax{i}-f_ref{i}).^2;
end
toc;