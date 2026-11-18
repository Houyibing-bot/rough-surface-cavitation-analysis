%% 多参数可控高斯粗糙面生成 + 主流粗糙面算法定量对比
% 参考文献算法：
% 1. 本文：迭代校正多参数高斯随机过程(改进Jawaid)
% 2. 原始Jawaid高斯(无迭代校正)
% 3. 傅里叶变换法(Fourier transform)
% 4. 分形函数法(Fractal)
% 5. 蒙特卡洛随机高斯(Monte‑Carlo)
% 评价指标：RMS相对误差eRMS、峰高相对误差ePeak、谷深相对误差eValley、综合误差Etotal、CPU耗时
clear; clc; close all;
restoredefaultpath;
rehash;
rng(0);  %% 固定随机种子，保证结果可复现

%% ====================== 全局设置 ======================
N = 256;                % 网格分辨率，和原文保持一致
L = 1;                  % 表面几何尺寸
betax = 20; betay = 20; % x/y方向相关长度
C = 1;                  % 功率谱常数
max_iterations = 15;    % 本文算法最大迭代次数
rms_tolerance = 0.01;   % RMS收敛容差
pv_tolerance  = 0.01;   % 峰谷差收敛容差

%% 论文3组测试工况 [RMS目标，目标最低谷，目标最高峰]
case_set = {
    {'Slight rough',    0.10,  -0.50,   0.50};
    {'Middle rough',    0.40,  -2.00,   2.00};
    {'Heavy rough',     0.80,  -4.00,   4.00};
};

method_list = {
    'Present_iter',...
    'Jawaid_original',...
    'Fourier',...
    'Fractal',...
    'MonteCarlo'
};
n_case = length(case_set);
n_method = length(method_list);

%% 预分配结果存储：每列：eRMS(%),ePeak(%),eValley(%),Etotal(%),time(s)
result_table = zeros(n_case, n_method*5);

%% ====================== 循环遍历所有工况 ======================
for ic = 1:n_case
    casename  = case_set{ic}{1};
    delta_tar = case_set{ic}{2};
    valley_tar= case_set{ic}{3};
    peak_tar  = case_set{ic}{4};
    PV_tar    = peak_tar - valley_tar;
    fprintf('\n========== 工况 %d / %d : %s ==========\n',ic,n_case,casename);
    
    %% ------------------- 1.本文迭代校正算法 Present_iter -------------------
    t0 = tic;
    [z_present, stat_present] = func_present_iter_surface(N,L,delta_tar,valley_tar,peak_tar,betax,betay,C,max_iterations,rms_tolerance,pv_tolerance);
    t1 = toc(t0);
    [eRMS_p, ePeak_p, eValley_p, Etotal_p] = calc_error(stat_present, delta_tar, peak_tar, valley_tar, PV_tar);
    result_table(ic,1:5) = [eRMS_p, ePeak_p, eValley_p, Etotal_p, t1];
    
    %% -------------------2.原始Jawaid(无迭代，只做FFT，不做峰谷RMS迭代校正)----------------
    t0 = tic;
    [z_jawaid, stat_jawaid] = func_jawaid_original(N,L,delta_tar,betax,betay,C);
    t1 = toc(t0);
    [eRMS_j, ePeak_j, eValley_j, Etotal_j] = calc_error(stat_jawaid, delta_tar, peak_tar, valley_tar, PV_tar);
    result_table(ic,6:10) = [eRMS_j, ePeak_j, eValley_j, Etotal_j, t1];
    
    %% -------------------3.傅里叶变换粗糙面 Fourier(仅约束RMS，峰谷自然随机输出)----------------
    t0 = tic;
    [z_ft, stat_ft] = func_fourier_surface(N,L,delta_tar,betax,betay);
    t1 = toc(t0);
    [eRMS_ft, ePeak_ft, eValley_ft, Etotal_ft] = calc_error(stat_ft, delta_tar, peak_tar, valley_tar, PV_tar);
    result_table(ic,11:15) = [eRMS_ft, ePeak_ft, eValley_ft, Etotal_ft, t1];
    
    %% -------------------4.分形粗糙面 Fractal ----------------
    t0 = tic;
    Hurst = 0.8; % 分形Hurst指数，粗糙面常用0.7‑0.9
    [z_frac, stat_frac] = func_fractal_surface(N,L,delta_tar,Hurst);
    t1 = toc(t0);
    [eRMS_fr, ePeak_fr, eValley_fr, Etotal_fr] = calc_error(stat_frac, delta_tar, peak_tar, valley_tar, PV_tar);
    result_table(ic,16:20) = [eRMS_fr, ePeak_fr, eValley_fr, Etotal_fr, t1];
    
    %% -------------------5.蒙特卡洛高斯粗糙面 Monte‑Carlo ----------------
    t0 = tic;
    [z_mc, stat_mc] = func_montecarlo_surface(N,L,delta_tar,betax,betay);
    t1 = toc(t0);
    [eRMS_mc, ePeak_mc, eValley_mc, Etotal_mc] = calc_error(stat_mc, delta_tar, peak_tar, valley_tar, PV_tar);
    result_table(ic,21:25) = [eRMS_mc, ePeak_mc, eValley_mc, Etotal_mc, t1];
    
    %% ========= 绘图：当前工况5种方法高度分布直方图 =========
    figure('Position',[50,50,1200,700],'Name',[casename,' Height distribution']);
    subplot(2,3,1); histogram(z_present(:),40); title('Present‑iter'); xlabel('Height'); grid on;
    subplot(2,3,2); histogram(z_jawaid(:),40); title('Jawaid original'); xlabel('Height'); grid on;
    subplot(2,3,3); histogram(z_ft(:),40); title('Fourier'); xlabel('Height'); grid on;
    subplot(2,3,4); histogram(z_frac(:),40); title('Fractal'); xlabel('Height'); grid on;
    subplot(2,3,5); histogram(z_mc(:),40); title('Monte‑Carlo'); xlabel('Height'); grid on;
    sgtitle([casename,' Height histogram comparison']);
end

%% ===================== 统计：3个工况取平均，用于论文表格Table1‑A =====================
avg_result = mean(result_table,1);
fprintf('\n\n==================== 多工况平均结果(用于论文表格) ====================\n');
fprintf('%-22s %-10s %-10s %-10s %-12s %-10s\n',...
    'Method','eRMS(%)','ePeak(%)','eValley(%)','Etotal(%)','Time(s)');

fprintf('%-22s %-10.2f %-10.2f %-10.2f %-12.2f %-10.2f\n',...
    'Present iterative‑corrected',avg_result(1),avg_result(2),avg_result(3),avg_result(4),avg_result(5));
fprintf('%-22s %-10.2f %-10.2f %-10.2f %-12.2f %-10.2f\n',...
    'Original Jawaid Gaussian',avg_result(6),avg_result(7),avg_result(8),avg_result(9),avg_result(10));
fprintf('%-22s %-10.2f %-10.2f %-10.2f %-12.2f %-10.2f\n',...
    'Fourier transform',avg_result(11),avg_result(12),avg_result(13),avg_result(14),avg_result(15));
fprintf('%-22s %-10.2f %-10.2f %-10.2f %-12.2f %-10.2f\n',...
    'Fractal function',avg_result(16),avg_result(17),avg_result(18),avg_result(19),avg_result(20));
fprintf('%-22s %-10.2f %-10.2f %-10.2f %-12.2f %-10.2f\n',...
    'Monte‑Carlo Gaussian',avg_result(21),avg_result(22),avg_result(23),avg_result(24),avg_result(25));

%% 导出csv，可直接复制到Excel生成论文Table1‑A
out_table = table(...
    avg_result(1),avg_result(2),avg_result(3),avg_result(4),avg_result(5),...
    avg_result(6),avg_result(7),avg_result(8),avg_result(9),avg_result(10),...
    avg_result(11),avg_result(12),avg_result(13),avg_result(14),avg_result(15),...
    avg_result(16),avg_result(17),avg_result(18),avg_result(19),avg_result(20),...
    avg_result(21),avg_result(22),avg_result(23),avg_result(24),avg_result(25),...
    'VariableNames',...
    {'Pre_eRMS','Pre_ePeak','Pre_eValley','Pre_Etotal','Pre_time',...
    'Jaw_eRMS','Jaw_ePeak','Jaw_eValley','Jaw_Etotal','Jaw_time',...
    'FT_eRMS','FT_ePeak','FT_eValley','FT_Etotal','FT_time',...
    'Frac_eRMS','Frac_ePeak','Frac_eValley','Frac_Etotal','Frac_time',...
    'MC_eRMS','MC_ePeak','MC_eValley','MC_Etotal','MC_time'});
writetable(out_table,'rough_compare_table.csv');
fprintf('\nCSV文件 rough_compare_table.csv 已经保存，可以导入Excel制作论文表格\n');

%% ===================== 子函数定义 =====================
%% 1.本文迭代校正多参数高斯粗糙面
function [z_final, stat_out] = func_present_iter_surface(N,L,delta_tar,lowest_valley,highest_peak,betax,betay,C,max_iter,rms_tol,pv_tol)
    [x,y] = meshgrid(linspace(0,L,N));
    target_PV = highest_peak - lowest_valley;
    eta = randn(N,N);
    A = fft2(eta);
    R = delta_tar^2 * exp(-2.3*sqrt((x/(L*betax/N)).^2 + (y/(L*betay/N)).^2));
    Gz = fft2(R);
    H = sqrt(abs(Gz)./C);
    Z = H .* A;
    z = real(ifft2(Z));
    
    z_final = z;
    iter = 1; converged = false;
    while ~converged && iter <= max_iter
        curr_mean = mean(z_final(:));
        curr_min  = min(z_final(:));
        curr_max  = max(z_final(:));
        curr_PV   = curr_max - curr_min;
        curr_rms  = std(z_final(:));
        
        % 步骤1：缩放峰谷差
        pv_scale = target_PV / curr_PV;
        z_scaled = (z_final - curr_mean)*pv_scale + curr_mean;
        % 步骤2：平移谷底至目标最低谷
        z_shifted = z_scaled - min(z_scaled(:)) + lowest_valley;
        % 步骤3：缩放RMS
        rms_scale = delta_tar / std(z_shifted(:));
        z_final = (z_shifted - mean(z_shifted(:)))*rms_scale + mean(z_shifted(:));
        
        % 更新状态
        f_min = min(z_final(:));
        f_max = max(z_final(:));
        f_PV = f_max - f_min;
        f_rms = std(z_final(:));
        err_rms = abs(f_rms - delta_tar);
        err_pv  = abs(f_PV - target_PV);
        if err_rms < rms_tol*delta_tar && err_pv < pv_tol*target_PV
            converged = true;
        end
        iter = iter + 1;
    end
    stat_out.rms = std(z_final(:));
    stat_out.peak = max(z_final(:));
    stat_out.valley = min(z_final(:));
    stat_out.PV = max(z_final(:))-min(z_final(:));
end

%% 2.原始Jawaid 无迭代校正，只做FFT，不做峰谷RMS迭代
function [z, stat_out] = func_jawaid_original(N,L,delta_tar,betax,betay,C)
    [x,y] = meshgrid(linspace(0,L,N));
    eta = randn(N,N);
    A = fft2(eta);
    R = delta_tar^2 * exp(-2.3*sqrt((x/(L*betax/N)).^2 + (y/(L*betay/N)).^2));
    Gz = fft2(R);
    H = sqrt(abs(Gz)./C);
    Z = H .* A;
    z = real(ifft2(Z));
    stat_out.rms = std(z(:));
    stat_out.peak = max(z(:));
    stat_out.valley = min(z(:));
    stat_out.PV = max(z(:))-min(z(:));
end

%% 3.传统傅里叶变换粗糙面（仅约束RMS，峰谷随机输出）
function [z, stat_out] = func_fourier_surface(N,L,delta_tar,betax,betay)
    [x,y] = meshgrid(linspace(0,L,N));
    noise = randn(N,N);
    fft_noise = fft2(noise);
    kx = 2*pi/L * ((0:N-1)-N/2);
    ky = 2*pi/L * ((0:N-1)-N/2);
    [KX,KY] = meshgrid(kx,ky);
    k_mod = sqrt(KX.^2+KY.^2);
    k_mod(k_mod==0)=eps;
    % 功率谱指数型
    S = delta_tar^2 * betax*betay ./ (1 + (k_mod.*betax).^2).^(1.5);
    Amp = sqrt(S);
    z = real(ifft2( fft_noise .* Amp ));
    z = z/std(z(:)) * delta_tar;
    stat_out.rms = std(z(:));
    stat_out.peak = max(z(:));
    stat_out.valley = min(z(:));
    stat_out.PV = max(z(:))-min(z(:));
end

%% 4.分形粗糙面（fractal，Weierstrass‑Mandelbrot分形）
function [z, stat_out] = func_fractal_surface(N,L,delta_tar,Hurst)
    [x,y] = meshgrid(linspace(0,L,N));
    H = Hurst;
    nf = 8; % 分形谐波阶数
    z = zeros(N,N);
    for n = 1:nf
        freq = 2^n;
        phase1 = 2*pi*rand();
        phase2 = 2*pi*rand();
        amp = freq^(-H);
        z = z + amp * (cos(freq*x*2*pi + phase1).*cos(freq*y*2*pi+phase2));
    end
    z = z/std(z(:)) * delta_tar;
    stat_out.rms = std(z(:));
    stat_out.peak = max(z(:));
    stat_out.valley = min(z(:));
    stat_out.PV = max(z(:))-min(z(:));
end

%%5.蒙特卡洛法粗糙面：高斯随机场，相关长度，直接生成，仅归一RMS
function [z, stat_out] = func_montecarlo_surface(N,L,delta_tar,betax,betay)
    [x,y] = meshgrid(linspace(0,L,N));
    eta = randn(N,N);
    R = exp(-2.3*sqrt((x/(L*betax/N)).^2 + (y/(L*betay/N)).^2));
    fftR = fft2(R);
    fftEta = fft2(eta);
    z = real(ifft2( sqrt(abs(fftR)).*fftEta ));
    z = z/std(z(:)) * delta_tar;
    stat_out.rms = std(z(:));
    stat_out.peak = max(z(:));
    stat_out.valley = min(z(:));
    stat_out.PV = max(z(:))-min(z(:));
end

%% 误差计算子函数：计算相对误差、综合几何误差
function [eRMS, ePeak, eValley, Etotal] = calc_error(stat_actual, tar_rms, tar_peak, tar_valley, tar_PV)
    eRMS = abs(stat_actual.rms - tar_rms)/tar_rms *100;
    ePeak = abs(stat_actual.peak - tar_peak)/abs(tar_peak)*100;
    eValley = abs(stat_actual.valley - tar_valley)/abs(tar_valley)*100;
    Etotal = (eRMS + ePeak + eValley)/3.0;
end
