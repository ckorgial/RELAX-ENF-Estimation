function IF = func_STFT_multi_tone_search_relax(signal,fs,window_dur,step_size_dur,fc,bound,FFT_res_factor)
window_length   = window_dur*fs;
window_func     = rectwin(window_length)';
step_size       = step_size_dur*fs;
NFFT            = FFT_res_factor*fs;
% signal          = [signal, zeros(1, window_length-1)];
window_pos      = 1:step_size:(length(signal)-window_length+1);
IF              = zeros(1,length(window_pos)); 
%% set harmonic search region
search_region_1st = round((50-bound(1)/2)*FFT_res_factor):round((50+bound(1)/2)*FFT_res_factor); % fundamental IF search region
length_per_band   = length(search_region_1st);
search_region     = kron(search_region_1st,(fc/50))';

% RELAX Algorithm Initialization
max_K = 10;
beta_hat = zeros(max_K, 1);
omega_hat = zeros(max_K, 1);
initial_omega = zeros(1, length(window_pos));

for i = 1:length(window_pos)
    current_window  = signal(window_pos(i):window_pos(i)+window_length-1) .* window_func;
    temp            = fft(signal(window_pos(i):window_pos(i)+window_length-1).*window_func,NFFT);
    HalfTempFFT     = temp(1:end/2);
    absHalfTempFFT  = abs(HalfTempFFT).';
    fbin_candidate  = absHalfTempFFT(search_region);
    fbin_candidate  = reshape(fbin_candidate,[length(fc),length_per_band]);
    weighted_fbin   = sum(fbin_candidate.^2,1);
    ValueMax        = max(weighted_fbin);
    PeakLoc         = search_region_1st((weighted_fbin==ValueMax(1)));
    initial_omega(i)= 2 * pi * (PeakLoc - 1) / NFFT;

    % RELAX Algorithm
    for K = 1:max_K
        for iter = 1:10
            for k = 1:K
                y_k = current_window;
                for j = 1:k-1
                    y_k = y_k - beta_hat(j) * exp(1j * omega_hat(j) * (0:length(current_window)-1));
                end
                % Grid search
                freq_range = linspace(initial_omega(i) - 0.0002, initial_omega(i) + 0.0002, 100);
                g_k_min = inf;
                for f_idx = 1:length(freq_range)
                    test_freq = freq_range(f_idx);
                    test_beta = y_k * exp(-1j * test_freq * (0:length(y_k)-1))'/length(y_k);
                    g_k_test = sum(abs(y_k - test_beta * exp(1j * test_freq * (0:length(y_k)-1))).^2);
                    if g_k_test < g_k_min
                        g_k_min = g_k_test;
                        omega_hat(k) = test_freq;
                        beta_hat(k) = test_beta;
                    end
                end
            end
            if g_k_min < 1e-6
                break;
            end
        end
    end
    [~, idx] = min(abs(omega_hat - initial_omega(i)));
    IF(i) = omega_hat(idx) * fs / (2 * pi)*2; 
end
norm_fc            = 100;
norm_bound         = 100*bound(1)/fc(1);
IF(IF<norm_fc-norm_bound)=norm_fc-norm_bound;IF(IF>norm_fc+norm_bound)=norm_fc+norm_bound;
end