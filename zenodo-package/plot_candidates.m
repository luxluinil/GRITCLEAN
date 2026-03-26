%Script to generate plots

%Define thresholds (change accordingly)
detection_threshold = 9; %set as 9 or 8
complexmass_threshold = 0.9;
deltatoa_threshold = 0.15;
deltasnr_threshold = 0.1;
flag = 1;
%-----------DO NOT CHANGE BELOW------------------------
if detection_threshold == 8
    deltasnr_threshold = 0.07;
end
thresholds = [detection_threshold, complexmass_threshold, deltatoa_threshold, deltasnr_threshold];

%Get All Candidates with real snr injections
posfile = 'Output/positiverun_candidates_w_siginjections.txt'; 
negfile = 'Output/negativerun_candidates_w_siginjections.txt';

[glitch_status,cidxs, toadiffs, snrdiffs, snrs,noise_pos, noise_neg] = run_GRITCLEAN(posfile, negfile, thresholds);

%Get all low snr signal injection candidates

lowposfile = 'Output/lowsnr_positiverun_candidates.txt';
lownegfile = 'Output/lowsnr_negativerun_candidates.txt';

[lowsnr_glitch_status,lowsnr_cidxs, lowsnr_toadiffs, lowsnr_snrdiffs, low_snrs,~,~] = run_GRITCLEAN(lowposfile, lownegfile, thresholds);


%Get all high snr signal injection candidates
highposfile = 'Output/highsnr_positiverun_candidates.txt';
highnegfile = 'Output/highsnr_negativerun_candidates.txt';

[highsnr_glitch_status,highsnr_cidxs, highsnr_toadiffs, highsnr_snrdiffs,high_snrs,~,~] = run_GRITCLEAN(highposfile, highnegfile, thresholds);


%Get mass gap injections

massgapposfile = 'Output/GVS_massgap_allsnrs_siginj_positive_sorted.txt';
massgapnegfile = 'Output/GVS_massgap_allsnrs_siginj_negative_sorted.txt';

[massgap_glitch_status,massgap_cidxs, massgap_toadiffs, massgap_snrdiffs,massgap_snrs,~,~] = run_GRITCLEAN(massgapposfile, massgapnegfile, thresholds);


%Get [25,40] M_sun injections

Msunposfile = 'Output/GVS_veryhighmass_25to40Msun_siginj_positive_sorted.txt';
Msunnegfile = 'Output/GVS_veryhighmass_25to40Msun_siginj_negative_sorted.txt';

[Msun_glitch_status,Msun_cidxs, Msun_toadiffs, Msun_snrdiffs,Msun_snrs,~,~] = run_GRITCLEAN(Msunposfile, Msunnegfile, thresholds);

%[25,40] M_sun injections Low SNR = [10,13]
Msunposfile_low = 'Output/GVS_veryhighmass_25to40Msun_lowsnr_siginj_positive_sorted.txt';
Msunnegfile_low = 'Output/GVS_veryhighmass_25to40Msun_lowsnr_siginj_negative_sorted.txt';

[Msun_glitch_status_low,Msun_cidxs_low, Msun_toadiffs_low, Msun_snrdiffs_low,Msun_snrs_low,~,~] = run_GRITCLEAN(Msunposfile_low, Msunnegfile_low, thresholds);

%%[25,40] M_sun injections High SNR = [100,500]
Msunposfile_high = 'Output/GVS_veryhighmass_25to40Msun_highsnr_siginj_positive_sorted.txt';
Msunnegfile_high = 'Output/GVS_veryhighmass_25to40Msun_highsnr_siginj_negative_sorted.txt';

[Msun_glitch_status_high,Msun_cidxs_high, Msun_toadiffs_high, Msun_snrdiffs_high,Msun_snrs_high,~,~] = run_GRITCLEAN(Msunposfile_high, Msunnegfile_high, thresholds);

%Get Very high mass injections: [40,80] and [80,110] Msun
vhighmass_posfile = 'Output/GVS_veryhighmass_allsnrs_siginj_positive_sorted.txt';
vhighmass_negfile = 'Output/GVS_veryhighmass_allsnrs_siginj_negative_sorted.txt';

[vhm_glitch_status, vhm_cidxs, vhm_toadiffs, vhm_snrdiffs, vhm_snrs, ~, ~] = run_GRITCLEAN(vhighmass_posfile, vhighmass_negfile, thresholds);

%Get PyCBC IMRPHenomXHM injections: [25,40] Msun and [40,80], [80,110] Msun
pycbc_posfile = 'Output/GVS_pycbc_siginj_positive_sorted.txt';
pycbc_negfile = 'Output/GVS_pycbc_siginj_negative_sorted.txt';

[pycbc_glitch_status, pycbc_cidxs, pycbc_toadiffs, pycbc_snrdiffs, pycbc_snrs, ~, ~] = run_GRITCLEAN(pycbc_posfile, pycbc_negfile, thresholds);

%Get PyCBC Lower Mass (Original mass ranges) IMRPHenomXHM injections: 
pycbc_lower_posfile = 'Output/GVS_pycbc_lowermass_siginj_positive_sorted.txt';
pycbc_lower_negfile = 'Output/GVS_pycbc_lowermass_siginj_negative_sorted.txt';

[pycbc_lower_glitch_status, pycbc_lower_cidxs, pycbc_lower_toadiffs, pycbc_lower_snrdiffs, pycbc_lower_snrs, ~, ~] = run_GRITCLEAN(pycbc_lower_posfile, pycbc_lower_negfile, thresholds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get Absolute differences in SNR and TOA for noise 
noisetoadiffs = abs(noise_pos(:,end-1) - noise_neg(:,end-1));
noisesnrdiffs = abs((noise_pos(:,end) - noise_neg(:,end))./noise_pos(:,end));

%Load and get simulated noise events
%Noise Files
noise_pos_file = 'Output/GVS_noise_positiverun_sorted.txt';
noise_neg_file = 'Output/GVS_noise_negativerun_sorted.txt';
noise_pos_file_2 = 'Output/GVS_noiseonly_v2_sorted.txt';
noise_neg_file_2 = 'Output/GVS_noiseonly_negative_v2_sorted.txt';
%Add the extra simulated noise estimates
posnoise = textread(noise_pos_file, '%s', 'delimiter', '\n');
negnoise = textread(noise_neg_file, '%s', 'delimiter', '\n');
numnoise = length(posnoise);

simulatednoisesnrdiff = [];
simulatednoisetoadiff = [];
simulatednoisesnrs = [];
for i = 1:numnoise
    negnoisevals = str2num(negnoise{i});
    posnoisevals = str2num(posnoise{i});
    negnoiseTOA = negnoisevals(end - 1);
    negnoiseSNR = negnoisevals(end);
    posnoiseTOA = posnoisevals(end - 1);
    posnoiseSNR = posnoisevals(end);

    relsnrnoisediff = (posnoiseSNR - negnoiseSNR)./posnoiseSNR;
    reltoanoisediff = posnoiseTOA - negnoiseTOA;
    
%     noisesnrs = [noisesnrs, posnoiseSNR];
    simulatednoisesnrs = [simulatednoisesnrs, posnoiseSNR];
    simulatednoisesnrdiff = [simulatednoisesnrdiff, relsnrnoisediff];
    simulatednoisetoadiff =[simulatednoisetoadiff, reltoanoisediff];
end

abssimnoiseSNRdiff = abs(simulatednoisesnrdiff);
abssimnoiseTOAdiff = abs(simulatednoisetoadiff);


%Add the second run of noise estimates

posnoise_2 = textread(noise_pos_file_2, '%s', 'delimiter', '\n');
negnoise_2 = textread(noise_neg_file_2, '%s', 'delimiter', '\n');
numnoise = length(posnoise_2);

simulatednoisesnrdiff_2 = [];
simulatednoisetoadiff_2 = [];
simulatednoisesnrs_2 = [];
for i = 1:numnoise
    negnoisevals_2 = str2num(negnoise_2{i});
    posnoisevals_2 = str2num(posnoise_2{i});
    negnoiseTOA_2 = negnoisevals_2(end - 1);
    negnoiseSNR_2 = negnoisevals_2(end);
    posnoiseTOA_2 = posnoisevals_2(end - 1);
    posnoiseSNR_2 = posnoisevals_2(end);

    relsnrnoisediff_2 = (posnoiseSNR_2 - negnoiseSNR_2)./posnoiseSNR_2;
    reltoanoisediff_2 = posnoiseTOA_2 - negnoiseTOA_2;
    
%     noisesnrs = [noisesnrs, posnoiseSNR];
    simulatednoisesnrs_2 = [simulatednoisesnrs_2, posnoiseSNR_2];
    simulatednoisesnrdiff_2 = [simulatednoisesnrdiff_2, relsnrnoisediff_2];
    simulatednoisetoadiff_2 =[simulatednoisetoadiff_2, reltoanoisediff_2];
end

abssimnoiseSNRdiff_2 = abs(simulatednoisesnrdiff_2);
abssimnoiseTOAdiff_2 = abs(simulatednoisetoadiff_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Segregate all glitches based on veto steps
chirplength_vetoed_glitches_idxs = [];
chirplength_vetoed_glitches_snrs = [];
chirplength_vetoed_glitches_toadiffs = [];
chirplength_vetoed_glitches_snrdiffs = [];
complexmass_vetoed_glitches_idxs = [];
complexmass_vetoed_glitches_snrs = [];
complexmass_vetoed_glitches_toadiffs = [];
complexmass_vetoed_glitches_snrdiffs = [];
SP_bypassed_glitches_idxs = [];
SP_bypassed_glitches_snrs = [];
SP_bypassed_glitches_toadiffs = [];
SP_bypassed_glitches_snrdiffs = [];
SN_vetoed_glitches_idxs = [];
SN_vetoed_glitches_snrs = [];
SN_vetoed_glitches_toadiffs = [];
SN_vetoed_glitches_snrdiffs = [];

%Glitch params
glitchidxs = [];
glitchsnrs = [];
glitch_toadiffs = [];
glitch_snrdiffs = [];

%Injected Signal (Real SNRs) params
realsnr_cidxs = [];
realsnr_glitch_status = [];
real_snrs = [];
realsnr_toadiffs = [];
realsnr_snrdiffs = [];

%Load Signal injection segments
S = load('signalinjectionsegments.mat');

realsnr_sigsegments = S.realsnrsigsegments;
lowsnrsigsegments = S.lowsnrsigsegments;
highsnrsigsegments = S.highsnrsigsegments;
veryhighmasssigsegments = S.veryhighmasssigsegments;
pycbcsigsegments = S.pycbcsigsegments;
pycbclowersigsegments = S.pycbclowersigsegments;

for i = 1:length(glitch_status)
    glitch = glitch_status(i,:);

    if ismember(cidxs(i),realsnr_sigsegments)
        realsnr_cidxs = [realsnr_cidxs,cidxs(i)];
        realsnr_glitch_status = [realsnr_glitch_status; glitch_status(i)];
        real_snrs = [real_snrs,snrs(i)];
        realsnr_toadiffs = [realsnr_toadiffs,toadiffs(i)];
        realsnr_snrdiffs = [realsnr_snrdiffs,snrdiffs(i)];
    else
        glitchidxs = [glitchidxs,cidxs(i)];
        glitchsnrs = [glitchsnrs,snrs(i)];
        glitch_toadiffs = [glitch_toadiffs,toadiffs(i)];
        glitch_snrdiffs = [glitch_snrdiffs,snrdiffs(i)];
        if glitch(2) == 1
            chirplength_vetoed_glitches_idxs = [chirplength_vetoed_glitches_idxs,cidxs(i)];
            chirplength_vetoed_glitches_snrs = [chirplength_vetoed_glitches_snrs,snrs(i)];
            chirplength_vetoed_glitches_toadiffs = [chirplength_vetoed_glitches_toadiffs,toadiffs(i)];
            chirplength_vetoed_glitches_snrdiffs = [chirplength_vetoed_glitches_snrdiffs,snrdiffs(i)];
        elseif glitch(2) == 2
            complexmass_vetoed_glitches_idxs = [complexmass_vetoed_glitches_idxs,cidxs(i)];
            complexmass_vetoed_glitches_snrs = [complexmass_vetoed_glitches_snrs,snrs(i)];
            complexmass_vetoed_glitches_toadiffs = [complexmass_vetoed_glitches_toadiffs,toadiffs(i)];
            complexmass_vetoed_glitches_snrdiffs = [complexmass_vetoed_glitches_snrdiffs,snrdiffs(i)];
        else
            SP_bypassed_glitches_idxs = [SP_bypassed_glitches_idxs,cidxs(i)];
            SP_bypassed_glitches_snrs = [SP_bypassed_glitches_snrs,snrs(i)];
            SP_bypassed_glitches_toadiffs = [SP_bypassed_glitches_toadiffs,toadiffs(i)];
            SP_bypassed_glitches_snrdiffs = [SP_bypassed_glitches_snrdiffs,snrdiffs(i)];
            if glitch(2) == 3
                SN_vetoed_glitches_idxs = [SN_vetoed_glitches_idxs,cidxs(i)];
                SN_vetoed_glitches_snrs = [SN_vetoed_glitches_snrs,snrs(i)];
                SN_vetoed_glitches_toadiffs = [SN_vetoed_glitches_toadiffs,toadiffs(i)];
                SN_vetoed_glitches_snrdiffs = [SN_vetoed_glitches_snrdiffs,snrdiffs(i)];
            end
        end
    end


end

if detection_threshold == 9
    type1_anomaly_glitch_idxs = [322,336,662,943];
    type2_anomaly_glitch_idxs = [620,621];
    type3_anomaly_glitch_idxs = [630,772];
    type4_anomaly_glitch_idxs = 517;
end

%Segregate signal injections into Low, Intermediate and High Mass
%injections
lowmass_toadiffs = [];
lowmass_snrdiffs = [];
lowmass_snrs = [];
intermmass_toadiffs = [];
intermmass_snrdiffs = [];
intermmass_snrs = [];
highmass_toadiffs = [];
highmass_snrdiffs = [];
highmass_snrs = [];

for i = 1:length(realsnr_cidxs)
    if ismember(realsnr_cidxs(i), realsnr_sigsegments(1:66))
        lowmass_toadiffs = [lowmass_toadiffs,realsnr_toadiffs(i)];
        lowmass_snrdiffs = [lowmass_snrdiffs,realsnr_snrdiffs(i)];
        lowmass_snrs = [lowmass_snrs,real_snrs(i)];
    elseif ismember(realsnr_cidxs(i), realsnr_sigsegments(67:132))
        intermmass_toadiffs = [intermmass_toadiffs,realsnr_toadiffs(i)];
        intermmass_snrdiffs = [intermmass_snrdiffs,realsnr_snrdiffs(i)];
        intermmass_snrs = [intermmass_snrs,real_snrs(i)];
    elseif ismember(realsnr_cidxs(i), realsnr_sigsegments(133:end))
        highmass_toadiffs = [highmass_toadiffs,realsnr_toadiffs(i)];
        highmass_snrdiffs = [highmass_snrdiffs,realsnr_snrdiffs(i)];
        highmass_snrs = [highmass_snrs,real_snrs(i)];
    end
end

for i = 1:length(lowsnr_cidxs)
    if ismember(lowsnr_cidxs(i), lowsnrsigsegments(1:66))
        lowmass_toadiffs = [lowmass_toadiffs,lowsnr_toadiffs(i)];
        lowmass_snrdiffs = [lowmass_snrdiffs,lowsnr_snrdiffs(i)];
        lowmass_snrs = [lowmass_snrs,low_snrs(i)];
    elseif ismember(lowsnr_cidxs(i), lowsnrsigsegments(67:132))
        intermmass_toadiffs = [intermmass_toadiffs,lowsnr_toadiffs(i)];
        intermmass_snrdiffs = [intermmass_snrdiffs,lowsnr_snrdiffs(i)];
        intermmass_snrs = [intermmass_snrs,low_snrs(i)];
    elseif ismember(lowsnr_cidxs(i), lowsnrsigsegments(133:end))
        highmass_toadiffs = [highmass_toadiffs,lowsnr_toadiffs(i)];
        highmass_snrdiffs = [highmass_snrdiffs,lowsnr_snrdiffs(i)];
        highmass_snrs = [highmass_snrs,low_snrs(i)];
    end
end

for i = 1:length(highsnr_cidxs)
    if ismember(highsnr_cidxs(i), highsnrsigsegments(1:20))
        lowmass_toadiffs = [lowmass_toadiffs,highsnr_toadiffs(i)];
        lowmass_snrdiffs = [lowmass_snrdiffs,highsnr_snrdiffs(i)];
        lowmass_snrs = [lowmass_snrs,high_snrs(i)];
    elseif ismember(highsnr_cidxs(i), highsnrsigsegments(21:40))
        intermmass_toadiffs = [intermmass_toadiffs,highsnr_toadiffs(i)];
        intermmass_snrdiffs = [intermmass_snrdiffs,highsnr_snrdiffs(i)];
        intermmass_snrs = [intermmass_snrs,high_snrs(i)];
    elseif ismember(highsnr_cidxs(i), highsnrsigsegments(41:end))
        highmass_toadiffs = [highmass_toadiffs,highsnr_toadiffs(i)];
        highmass_snrdiffs = [highmass_snrdiffs,highsnr_snrdiffs(i)];
        highmass_snrs = [highmass_snrs,high_snrs(i)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Print Statistics
disp(['Detection threshold = ',num2str(detection_threshold)]);
disp(['Complex mass ratio threshold = ',num2str(complexmass_threshold)]);
disp(['Total number of glitch events = ', num2str(length(glitchidxs))]);
str1 = ['Number of chirplength vetoed glitch events = ', num2str(length(chirplength_vetoed_glitches_idxs)),' = ',num2str(100*length(chirplength_vetoed_glitches_idxs)/length(glitchidxs)),' % of total'];
disp(str1);

str2 = ['Number of complex mass vetoed glitch events = ', num2str(length(complexmass_vetoed_glitches_idxs)),' = ',num2str(100*length(complexmass_vetoed_glitches_idxs)/length(glitchidxs)),' % of total'];
disp(str2);

SP_vetoed_glitchnum = length(chirplength_vetoed_glitches_idxs) + length(complexmass_vetoed_glitches_idxs);
str3 = ['Number of SP vetoed glitch events = ', num2str(SP_vetoed_glitchnum),' = ',num2str(100*SP_vetoed_glitchnum/length(glitchidxs)),' % of total'];
disp(str3);
str4 = ['Number of SP bypass glitch events = ', num2str(length(SP_bypassed_glitches_snrdiffs)),' = ',num2str(100*length(SP_bypassed_glitches_snrdiffs)/length(glitchidxs)),' % of total'];
disp(str4);
str5 = ['Number of SN vetoed glitch events = ', num2str(length(SN_vetoed_glitches_idxs)),' = ',num2str(100*length(SN_vetoed_glitches_idxs)/length(glitchidxs)),' % of total'];
disp(str5);
total_vetoed_glitchnum = length(SN_vetoed_glitches_idxs) + SP_vetoed_glitchnum;
total_missed_glitchnum = length(glitchidxs) - total_vetoed_glitchnum;
str6 = ['Total number of VETOED glitch events = ', num2str(total_vetoed_glitchnum),' out of ',num2str(length(glitchidxs)),' = ',num2str(100*total_vetoed_glitchnum/length(glitchidxs)),' % of total'];
disp(str6);
str6 = ['Total number of MISSED glitch events = ', num2str(total_missed_glitchnum),' out of ',num2str(length(glitchidxs)),' = ',num2str(100*total_missed_glitchnum/length(glitchidxs)),' % of total'];
disp(str6);
%Injected signal statistics
numinj = 454 + 140 + 18;
%total_detected_injections = length(real_snrs) + length(low_snrs) + length(high_snrs) + length(Msun_snrs_high) + length(Msun_snrs_low) + length(Msun_snrs) + length(massgap_snrs);
total_detected_injections = length(real_snrs) + length(low_snrs) + length(high_snrs);
detection_threshold_missed_siginjs = numinj - total_detected_injections;

%Number of signal injections classified as glitches
%missed_signal_injections = sum(lowsnr_glitch_status(:,1)) + sum(highsnr_glitch_status(:,1)) + sum(realsnr_glitch_status(:,1)) + sum(Msun_glitch_status_high(:,1)) + sum(Msun_glitch_status(:,1)) + sum(Msun_glitch_status_low(:,1)) + sum(massgap_glitch_status(:,1));
missed_signal_injections = sum(lowsnr_glitch_status(:,1)) + sum(highsnr_glitch_status(:,1)) + sum(realsnr_glitch_status(:,1));
str5 = ['Number of signal injections missed due to detection threshold = ',num2str(detection_threshold_missed_siginjs), ' out of ',num2str(numinj)];
disp(str5);
str6 = ['Number of signal injections classified as glitches = ',num2str(missed_signal_injections),' = ',num2str(100*missed_signal_injections/total_detected_injections),'% of total detected signals'];
disp(str6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make Plots
if detection_threshold == 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot Histogram of Vetoed and Non-Vetoed SNRs
    sorted_nonvetoedglitchsnrs = sort(SP_bypassed_glitches_snrs);
    sorted_glitchsnrs = sort(glitchsnrs);
    counts_glitchs = linspace(1,length(glitchsnrs),length(glitchsnrs));
    counts_nonvetoedglitchs = linspace(1,length(SP_bypassed_glitches_snrs),length(SP_bypassed_glitches_snrs));
    chirplengthbypassedglitchsnrs = [complexmass_vetoed_glitches_snrs, SP_bypassed_glitches_snrs];
    sorted_glitches_afterchirplengthveto = sort(chirplengthbypassedglitchsnrs);
    counts_afterchirplengthvetoed = linspace(1, length(chirplengthbypassedglitchsnrs ),length(chirplengthbypassedglitchsnrs ));
    figure; 
    plot(log10(sorted_glitchsnrs),length(glitchsnrs) - counts_glitchs,'LineWidth',4,'DisplayName','All Glitch Events');
    hold on; 
    plot(log10(sorted_glitches_afterchirplengthveto),length(chirplengthbypassedglitchsnrs) - counts_afterchirplengthvetoed,'LineWidth',4,'DisplayName','Glitches that bypass chirp length vetoes');
    hold on; 
    plot(log10(sorted_nonvetoedglitchsnrs), length(SP_bypassed_glitches_snrs) - counts_nonvetoedglitchs,'LineWidth',4,'DisplayName', 'Glitch Events that bypass $S_p$ vetoes'); 
    hold off;
    xlabel('$\log{(\hat{\rho}_P)}$','Interpreter','latex');
    ylabel('Counts');
    grid on;
    ax = gca; ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
    legend('FontSize',55,'Interpreter','latex'); axis tight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Make Plot of glitches colored according to log SNR
    sz = 100;
    numBins = 5;
    binEdges = linspace(min(log10(glitchsnrs)), max(log10(glitchsnrs)), numBins);
    [~, binIdx] = histc(log10(glitchsnrs), binEdges);
    colors = jet(numBins);
    figure;
    hold on;
    scatter(glitch_snrdiffs, glitch_toadiffs,sz,colors(binIdx, :),'filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold off;
    colormap(colors);
    clim([min(log10(glitchsnrs)), max(log10(glitchsnrs))]);
    c = colorbar('Ticks', binEdges, 'TickLabels', cellstr(num2str(binEdges', '%.2f')));
    c.FontSize = 20;
    grid on;
    xlabel('$|\Delta \rho|$','Interpreter','latex');
    ylabel('$|\Delta t_a|$','Interpreter','latex');
    ylim([0.0001,1000]);
    ax = gca;
    ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %3D Distribution Plot
    figure;
    sz = 10;
    sz3 = 20;
    hold on;
    plot3(glitch_snrdiffs, glitch_toadiffs,log10(glitchsnrs),'o','MarkerFaceColor','red','MarkerSize',sz,'MarkerEdgeColor','black','DisplayName','GLITCHES');
    hold on;
    plot3(realsnr_snrdiffs, realsnr_toadiffs,log10(real_snrs),'o','MarkerFaceColor','yellow','MarkerSize',sz,'MarkerEdgeColor','black','DisplayName','GLITCHES');
    hold on;
    plot3(lowsnr_snrdiffs, lowsnr_toadiffs,log10(low_snrs(:,1)),'o','MarkerFaceColor','cyan','MarkerSize',sz,'MarkerEdgeColor','black','DisplayName','GLITCHES');
    hold on;
    plot3(highsnr_snrdiffs, highsnr_toadiffs,log10(high_snrs(:,1)),'o','MarkerFaceColor','blue','MarkerSize',sz,'MarkerEdgeColor','black','DisplayName','GLITCHES');
    hold on;
%     plot3(massgap_snrdiffs, massgap_toadiffs,log10(massgap_snrs(:,1)),'o','MarkerFaceColor','black','MarkerSize',sz,'MarkerEdgeColor','black','DisplayName','GLITCHES');
%     hold on;
%     plot3(Msun_snrdiffs, Msun_toadiffs,log10(Msun_snrs(:,1)),'o','MarkerFaceColor','yellow','MarkerSize',sz,'MarkerEdgeColor','black','DisplayName','GLITCHES');
%     hold on;
%     plot3(Msun_snrdiffs_low, Msun_toadiffs_low,log10(Msun_snrs_low(:,1)),'o','MarkerFaceColor','cyan','MarkerSize',sz,'MarkerEdgeColor','black','DisplayName','GLITCHES');
%     hold on;
%     plot3(Msun_snrdiffs_high, Msun_toadiffs_high,log10(Msun_snrs_high(:,1)),'o','MarkerFaceColor','blue','MarkerSize',sz,'MarkerEdgeColor','black','DisplayName','GLITCHES');
    hold off;
    grid on;
    xlabel('$|\Delta \rho|$','Interpreter','latex');
    ylabel('$|\Delta t_a|$','Interpreter','latex');
    zlabel('$\log(\hat{\rho}_P)$','Interpreter','latex');
    ax = gca;
    ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;ax.ZAxis.FontSize = 40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %Delta TOA vs LogSNR Plot
    x = [deltatoa_threshold, deltatoa_threshold, 1000, 1000];
    y = [log10(detection_threshold), 4.5, 4.5, log10(detection_threshold)];
    sz = 100;
    figure;
    scatter(glitch_toadiffs, log10(glitchsnrs),sz,'red','filled','MarkerEdgeColor','black','DisplayName','GLITCHES');
    hold on;
    scatter(realsnr_toadiffs, log10(real_snrs),sz,'yellow','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(lowsnr_toadiffs, log10(low_snrs(:,1)),sz,'cyan','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(highsnr_toadiffs, log10(high_snrs(:,1)),sz,'blue','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
%     scatter(massgap_toadiffs, log10(massgap_snrs(:,1)),sz,'black','filled','MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
    scatter(Msun_toadiffs, log10(Msun_snrs(:,1)),sz,'yellow','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(Msun_toadiffs_low, log10(Msun_snrs_low(:,1)),sz,'cyan','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(Msun_toadiffs_high, log10(Msun_snrs_high(:,1)),sz,'blue','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    patch(x,y,'green','FaceAlpha',0.04,'HandleVisibility','off');
    hold on;
%     scatter(nan,nan,sz3,"black",'filled','MarkerEdgeColor','black','DisplayName','Mass-gap Injections'); 
%     hold on;
    scatter(nan,nan,sz3,"yellow",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,40]$'); 
    hold on;
    scatter(nan,nan,sz3,"cyan",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,13]$');
    hold on; 
    scatter(nan,nan,sz3,"blue",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [100,500]$');
    hold off;
    legend;
    grid on;
    xlabel('$|\Delta t_a|$','Interpreter','latex');
    ylabel('$\log(\hat{\rho}_P)$','Interpreter','latex');
    ax = gca;
    ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
    legend('FontSize',45,'Interpreter','latex');
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %Delta SNR vs LogSNR Plot
    x = [deltasnr_threshold, deltasnr_threshold, 1, 1];
    y = [log10(detection_threshold), 4.5, 4.5, log10(detection_threshold)];
    sz = 100;
    figure;
    scatter(glitch_snrdiffs, log10(glitchsnrs),sz,'red','filled','MarkerEdgeColor','black','DisplayName','GLITCHES');
    hold on;
    scatter(realsnr_snrdiffs, log10(real_snrs),sz,'yellow','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(lowsnr_snrdiffs, log10(low_snrs(:,1)),sz,'cyan','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(highsnr_snrdiffs, log10(high_snrs(:,1)),sz,'blue','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
%     scatter(massgap_snrdiffs, log10(massgap_snrs(:,1)),sz,'black','filled','MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(Msun_snrdiffs, log10(Msun_snrs(:,1)),sz,'yellow','filled','MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(Msun_snrdiffs_low, log10(Msun_snrs_low(:,1)),sz,'cyan','filled','MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(Msun_snrdiffs_high, log10(Msun_snrs_high(:,1)),sz,'blue','filled','MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
    patch(x,y,'green','FaceAlpha',0.04,'HandleVisibility','off');
    hold on;
%     scatter(nan,nan,sz3,"black",'filled','MarkerEdgeColor','black','DisplayName','Mass-gap Injections'); 
%     hold on;
    scatter(nan,nan,sz3,"yellow",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,40]$'); 
    hold on;
    scatter(nan,nan,sz3,"cyan",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,13]$');
    hold on; 
    scatter(nan,nan,sz3,"blue",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [100,500]$');
    hold off;
    legend;
    grid on;
    xlabel('$|\Delta \rho|$','Interpreter','latex');
    ylabel('$\log(\hat{\rho}_P)$','Interpreter','latex');
    ax = gca;
    ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
    [~, icons] = legend('FontSize',45,'Interpreter','latex');
    h = findobj(icons,'type','patch');
    set(h,'MarkerSize',15);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Vetoed glitch plot for detection threshold = 9

    sz = 100;
    sz2 = 200;
    % temp_snr_thresh = 0.0161;
    temp_snr_thresh = 0.0955;
    % temp_toa_thresh = 0.104;
    x = [deltasnr_threshold,1,1,deltasnr_threshold];
    y = [deltatoa_threshold, deltatoa_threshold, 1000, 1000];
    % x = [temp_snr_thresh,1,1,temp_snr_thresh];
    % y = [temp_toa_thresh, temp_toa_thresh, 1000, 1000];
    if flag
    
        %%%%%PyCBC Injections Plot%%%%%%%%%%%%%%%%%%%%%
        figure;
        rectangle_x = [0.0161, 1];
        rectangle_y = [0.104, 1000];
        rectangle('Position', [rectangle_x(1), rectangle_y(1), rectangle_x(2)-rectangle_x(1), rectangle_y(2)-rectangle_y(1)], ...
          'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
        hold on;
        signal_mrkrfacecolor = [224, 224, 224]./255;
        scatter(chirplength_vetoed_glitches_snrdiffs, chirplength_vetoed_glitches_toadiffs,sz,'MarkerFaceColor', [255, 204, 204]./255,'MarkerEdgeColor','black','DisplayName','CHIRPLENGTH VETOED GLITCHES');
        hold on;
        scatter(complexmass_vetoed_glitches_snrdiffs, complexmass_vetoed_glitches_toadiffs,sz,'MarkerFaceColor', [255,102,255]./255,'MarkerEdgeColor','black','DisplayName','COMPLEX MASS VETOED GLITCHES');
        hold on;
        scatter(SP_bypassed_glitches_snrdiffs, SP_bypassed_glitches_toadiffs,sz,'red','filled','MarkerEdgeColor','black','DisplayName','$S_P$ VETO BYPASSED GLITCHES');
        hold on;
        scatter(lowmass_snrdiffs(1:66), lowmass_toadiffs(1:66),sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        scatter(lowmass_snrdiffs(67:end-20), lowmass_toadiffs(67:end-20),sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        scatter(lowmass_snrdiffs(end-19:end), lowmass_toadiffs(end-19:end),sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        scatter(intermmass_snrdiffs(1:66), intermmass_toadiffs(1:66),sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        scatter(intermmass_snrdiffs(67:end-20), intermmass_toadiffs(67:end-20),sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        scatter(intermmass_snrdiffs(end-19:end), intermmass_toadiffs(end-19:end),sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        scatter(highmass_snrdiffs(1:66), highmass_toadiffs(1:66),sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        scatter(highmass_snrdiffs(67:end-20), highmass_toadiffs(67:end-20),sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        scatter(highmass_snrdiffs(end-19:end), highmass_toadiffs(end-19:end),sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        scatter(Msun_snrdiffs, Msun_toadiffs,sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        scatter(Msun_snrdiffs_low, Msun_toadiffs_low,sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        scatter(Msun_snrdiffs_high, Msun_toadiffs_high,sz,'MarkerFaceColor', signal_mrkrfacecolor,'MarkerEdgeColor','black','HandleVisibility','off');
        hold on;
        %Added PyCBC Lower Mass injections
        for i = 1:length(pycbc_lower_cidxs)
            if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(1:60)) 
                if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(1:20))
                    scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'cyan','filled','square','MarkerEdgeColor','black','HandleVisibility','off');
                hold on;
                end
                if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(21:40))
                    scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'yellow','filled','square','MarkerEdgeColor','black','HandleVisibility','off');
                hold on;
                end
                if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(41:60))
                    scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'blue','filled','square','MarkerEdgeColor','black','HandleVisibility','off');
                hold on;
                end
            end
            if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(61:120)) 
                if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(61:80))
                    scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'cyan','filled','^','MarkerEdgeColor','black','HandleVisibility','off');
                hold on;
                end
                if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(81:100))
                    scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'yellow','filled','^','MarkerEdgeColor','black','HandleVisibility','off');
                hold on;
                end
                if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(101:120))
                    scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'blue','filled','^','MarkerEdgeColor','black','HandleVisibility','off');
                hold on;
                end
            end
            if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(121:180)) 
                if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(121:140))
                    scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'cyan','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
                hold on;
                end
                if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(141:160))
                    scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'yellow','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
                hold on;
                end
                if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(161:180))
                    scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'blue','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
                hold on;
                end
            end
        end
        %Added PyCBC High Mass injections
        % for i = 1:length(pycbc_cidxs)
        %     if ismember(pycbc_cidxs(i),pycbcsigsegments(21:40)) || ismember(pycbc_cidxs(i),pycbcsigsegments(81:100))
        %         scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz,'yellow','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
        %         hold on;
        %     end
        %     if ismember(pycbc_cidxs(i),pycbcsigsegments(1:20)) || ismember(pycbc_cidxs(i),pycbcsigsegments(61:80))
        %         scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz,'cyan','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
        %         hold on;
        %     end
        %     if ismember(pycbc_cidxs(i),pycbcsigsegments(41:60)) || ismember(pycbc_cidxs(i),pycbcsigsegments(101:120))
        %         scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz,'blue','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
        %         hold on;
        %     end
        % end
        
        for i = 1:length(pycbc_cidxs)
            if ismember(pycbc_cidxs(i),pycbcsigsegments(1:60))
                if ismember(pycbc_cidxs(i),pycbcsigsegments(1:20))
                    scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'cyan','filled','pentagram','MarkerEdgeColor','black','HandleVisibility','off');
                    hold on;
                end
                if ismember(pycbc_cidxs(i),pycbcsigsegments(21:40))
                    scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'yellow','filled','pentagram','MarkerEdgeColor','black','HandleVisibility','off');
                    hold on;
                end
                if ismember(pycbc_cidxs(i),pycbcsigsegments(41:60))
                    scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'blue','filled','pentagram','MarkerEdgeColor','black','HandleVisibility','off');
                    hold on;
                end
            end
            if ismember(pycbc_cidxs(i),pycbcsigsegments(61:120))
                if ismember(pycbc_cidxs(i),pycbcsigsegments(61:80))
                    scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'cyan','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
                    hold on;
                end
                if ismember(pycbc_cidxs(i),pycbcsigsegments(81:100))
                    scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'yellow','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
                    hold on;
                end
                if ismember(pycbc_cidxs(i),pycbcsigsegments(101:120))
                    scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'blue','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
                    hold on;
                end
            end
        end
    
        %Add [25,40] M_sun injections in realistic SNRs
    %     scatter(massgap_snrdiffs, massgap_toadiffs,sz,'MarkerFaceColor', [224, 224, 224]./255','MarkerEdgeColor','black','HandleVisibility','off');
    %     hold on;
        
        patch(x,y,'green','FaceAlpha',0.04,'HandleVisibility','off');
        hold on;
    %     scatter(nan,nan,sz2,'^','MarkerEdgeColor','black','DisplayName','INTERMEDIATE MASS INJECTIONS'); 
        scatter(nan,nan,sz2,'square','MarkerEdgeColor','black','DisplayName','$m_1,m_2 \in [1.4,3] M_{\odot}$');
         hold on; 
         scatter(nan,nan,sz2,'^','MarkerEdgeColor','black','DisplayName','$m_1 \in [1.4,3] M_{\odot}$ and $m_2 \in [5,10] M_{\odot}$');
         hold on;
        scatter(nan,nan,sz2,'diamond','MarkerEdgeColor','black','DisplayName','$m_1,m_2 \in [10,25] M_{\odot}$');
         hold on;
         scatter(nan,nan,sz2,'pentagram','MarkerEdgeColor','black','DisplayName','$m_1,m_2 \in [25,40] M_{\odot}$');
         hold on;
         scatter(nan,nan,sz2,'hexagram','MarkerEdgeColor','black','DisplayName','$m_1 \in [40,80] M_{\odot}$ and $m_2 \in [80,110] M_{\odot}$');
        hold on;
         
    %     scatter(nan,nan,sz2,'black','filled','MarkerEdgeColor','black','DisplayName','Mass-gap Injections');
    %     hold on;
    %     scatter(nan,nan,sz2,'pentagram','MarkerEdgeColor','black','DisplayName','TYPE 1'); 
    %     hold on;
    %     scatter(nan,nan,sz2,'v','MarkerEdgeColor','black','DisplayName','TYPE 2');
    %     hold on; 
    %     scatter(nan,nan,sz2,'<','MarkerEdgeColor','black','DisplayName','TYPE 3');
    %     hold on;
    %     scatter(nan,nan,sz2,'*','MarkerEdgeColor','black','DisplayName','TYPE 4');
    %     hold on;
        scatter(nan,nan,sz3,"cyan",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,13]$');
        hold on; 
        scatter(nan,nan,sz3,"yellow",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,40]$'); 
        hold on;
        scatter(nan,nan,sz3,"blue",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [100,500]$');
        hold off;
        legend;
        grid on;
        xlabel('$|\Delta \rho|$','Interpreter','latex');
        ylabel('$|\Delta t_a|$','Interpreter','latex');
        xlim([1e-7, 1]);
        ax = gca;
        ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
        [leg, icons] = legend('FontSize',30,'Interpreter','latex');
        h = findobj(icons,'type','patch');
        set(h,'MarkerSize',15);
end

    figure;
    %Plot Chirplength Vetoed Glitches
    %Plot different types of outliers first
    for i = 1:length(type1_anomaly_glitch_idxs)
        if ismember(type1_anomaly_glitch_idxs(i),chirplength_vetoed_glitches_idxs)
            index = find(chirplength_vetoed_glitches_idxs == type1_anomaly_glitch_idxs(i));
            scatter(chirplength_vetoed_glitches_snrdiffs(index), chirplength_vetoed_glitches_toadiffs(index),sz2,'p','MarkerFaceColor', [255, 204, 204]./255,'MarkerEdgeColor','black','HandleVisibility', 'off');
            hold on;
            chirplength_vetoed_glitches_idxs(index) = [];
            chirplength_vetoed_glitches_snrs(index) = [];
            chirplength_vetoed_glitches_snrdiffs(index) = [];
            chirplength_vetoed_glitches_toadiffs(index) = [];
        end
    end
    for i = 1:length(type2_anomaly_glitch_idxs)
        if ismember(type2_anomaly_glitch_idxs(i),chirplength_vetoed_glitches_idxs)
            index = find(chirplength_vetoed_glitches_idxs == type2_anomaly_glitch_idxs(i));
            scatter(chirplength_vetoed_glitches_snrdiffs(index), chirplength_vetoed_glitches_toadiffs(index),sz2,'v','MarkerFaceColor', [255, 204, 204]./255,'MarkerEdgeColor','black','HandleVisibility', 'off');
            hold on;
            chirplength_vetoed_glitches_idxs(index) = [];
            chirplength_vetoed_glitches_snrs(index) = [];
            chirplength_vetoed_glitches_snrdiffs(index) = [];
            chirplength_vetoed_glitches_toadiffs(index) = [];
        end
    end
    for i = 1:length(type3_anomaly_glitch_idxs)
        if ismember(type3_anomaly_glitch_idxs(i),chirplength_vetoed_glitches_idxs)
            index = find(chirplength_vetoed_glitches_idxs == type3_anomaly_glitch_idxs(i));
            scatter(chirplength_vetoed_glitches_snrdiffs(index), chirplength_vetoed_glitches_toadiffs(index),sz2,'<','MarkerFaceColor', [255, 204, 204]./255,'MarkerEdgeColor','black','HandleVisibility', 'off');
            hold on;
            chirplength_vetoed_glitches_idxs(index) = [];
            chirplength_vetoed_glitches_snrs(index) = [];
            chirplength_vetoed_glitches_snrdiffs(index) = [];
            chirplength_vetoed_glitches_toadiffs(index) = [];
        end
    end
    for i = 1:length(type4_anomaly_glitch_idxs)
        if ismember(type4_anomaly_glitch_idxs(i),chirplength_vetoed_glitches_idxs)
            index = find(chirplength_vetoed_glitches_idxs == type4_anomaly_glitch_idxs(i));
            scatter(chirplength_vetoed_glitches_snrdiffs(index), chirplength_vetoed_glitches_toadiffs(index),sz2,'*','MarkerFaceColor', [255, 204, 204]./255,'MarkerEdgeColor', [255, 204, 204]./255, 'LineWidth',1.5,'HandleVisibility', 'off');
            hold on;
            chirplength_vetoed_glitches_idxs(index) = [];
            chirplength_vetoed_glitches_snrs(index) = [];
            chirplength_vetoed_glitches_snrdiffs(index) = [];
            chirplength_vetoed_glitches_toadiffs(index) = [];
        end
    end
    scatter(chirplength_vetoed_glitches_snrdiffs, chirplength_vetoed_glitches_toadiffs,sz,'MarkerFaceColor', [255, 204, 204]./255,'MarkerEdgeColor','black','DisplayName','CHIRPLENGTH VETOED GLITCHES');
    hold on;
    
    %Plot Complex Mass Vetoed Glitches
    %Plot Outliers first
    for i = 1:length(type1_anomaly_glitch_idxs)
        if ismember(type1_anomaly_glitch_idxs(i),complexmass_vetoed_glitches_idxs)
            index = find(complexmass_vetoed_glitches_idxs == type1_anomaly_glitch_idxs(i));
            scatter(complexmass_vetoed_glitches_snrdiffs(index), complexmass_vetoed_glitches_toadiffs(index),sz2,'p','MarkerFaceColor',[255,102,255]./255,'MarkerEdgeColor','black','HandleVisibility', 'off');
            hold on;
            complexmass_vetoed_glitches_idxs(index) = [];
            complexmass_vetoed_glitches_snrs(index) = [];
            complexmass_vetoed_glitches_snrdiffs(index) = [];
            complexmass_vetoed_glitches_toadiffs(index) = [];
        end
    end
    for i = 1:length(type2_anomaly_glitch_idxs)
        if ismember(type2_anomaly_glitch_idxs(i),complexmass_vetoed_glitches_idxs)
            index = find(complexmass_vetoed_glitches_idxs == type2_anomaly_glitch_idxs(i));
            scatter(complexmass_vetoed_glitches_snrdiffs(index), complexmass_vetoed_glitches_toadiffs(index),sz2,'v','MarkerFaceColor',[255,102,255]./255,'MarkerEdgeColor','black','HandleVisibility', 'off');
            hold on;
            complexmass_vetoed_glitches_idxs(index) = [];
            complexmass_vetoed_glitches_snrs(index) = [];
            complexmass_vetoed_glitches_snrdiffs(index) = [];
            complexmass_vetoed_glitches_toadiffs(index) = [];
        end
    end
    for i = 1:length(type3_anomaly_glitch_idxs)
        if ismember(type3_anomaly_glitch_idxs(i),complexmass_vetoed_glitches_idxs)
            index = find(complexmass_vetoed_glitches_idxs == type3_anomaly_glitch_idxs(i));
            scatter(complexmass_vetoed_glitches_snrdiffs(index), complexmass_vetoed_glitches_toadiffs(index),sz2,'<','MarkerFaceColor',[255,102,255]./255,'MarkerEdgeColor','black','HandleVisibility', 'off');
            hold on;
            complexmass_vetoed_glitches_idxs(index) = [];
            complexmass_vetoed_glitches_snrs(index) = [];
            complexmass_vetoed_glitches_snrdiffs(index) = [];
            complexmass_vetoed_glitches_toadiffs(index) = [];
        end
    end
    for i = 1:length(type4_anomaly_glitch_idxs)
        if ismember(type4_anomaly_glitch_idxs(i),complexmass_vetoed_glitches_idxs)
            index = find(complexmass_vetoed_glitches_idxs == type4_anomaly_glitch_idxs(i));
            scatter(complexmass_vetoed_glitches_idxs(index), chirplength_vetoed_glitches_toadiffs(index),sz2,'*','MarkerFaceColor',[255,102,255]./255,'MarkerEdgeColor', [255,102,255]./255, 'LineWidth',1.5,'HandleVisibility', 'off');
            hold on;
            complexmass_vetoed_glitches_idxs(index) = [];
            complexmass_vetoed_glitches_snrs(index) = [];
            complexmass_vetoed_glitches_snrdiffs(index) = [];
            complexmass_vetoed_glitches_toadiffs(index) = [];
        end
    end
    scatter(complexmass_vetoed_glitches_snrdiffs, complexmass_vetoed_glitches_toadiffs,sz,'MarkerFaceColor', [255,102,255]./255,'MarkerEdgeColor','black','DisplayName','COMPLEX MASS VETOED GLITCHES');
    hold on;
    %Plot S_P bypassed glitches 
    %Plot Outliers first
    for i = 1:length(type1_anomaly_glitch_idxs)
        if ismember(type1_anomaly_glitch_idxs(i),SP_bypassed_glitches_idxs)
            index = find(SP_bypassed_glitches_idxs == type1_anomaly_glitch_idxs(i));
            scatter(SP_bypassed_glitches_snrdiffs(index), SP_bypassed_glitches_toadiffs(index),sz2,'p','MarkerFaceColor', 'red','MarkerEdgeColor','black','HandleVisibility', 'off');
            hold on;
            SP_bypassed_glitches_idxs(index) = [];
            SP_bypassed_glitches_snrs(index) = [];
            SP_bypassed_glitches_snrdiffs(index) = [];
            SP_bypassed_glitches_toadiffs(index) = [];
        end
    end
    for i = 1:length(type2_anomaly_glitch_idxs)
        if ismember(type2_anomaly_glitch_idxs(i),SP_bypassed_glitches_idxs)
            index = find(SP_bypassed_glitches_idxs == type2_anomaly_glitch_idxs(i));
            scatter(SP_bypassed_glitches_snrdiffs(index), SP_bypassed_glitches_toadiffs(index),sz2,'v','MarkerFaceColor', 'red','MarkerEdgeColor','black','HandleVisibility', 'off');
            hold on;
            SP_bypassed_glitches_idxs(index) = [];
            SP_bypassed_glitches_snrs(index) = [];
            SP_bypassed_glitches_snrdiffs(index) = [];
            SP_bypassed_glitches_toadiffs(index) = [];
        end
    end
    for i = 1:length(type3_anomaly_glitch_idxs)
        if ismember(type3_anomaly_glitch_idxs(i),SP_bypassed_glitches_idxs)
            index = find(SP_bypassed_glitches_idxs == type3_anomaly_glitch_idxs(i));
            scatter(SP_bypassed_glitches_snrdiffs(index), SP_bypassed_glitches_toadiffs(index),sz2,'<','MarkerFaceColor', 'red','MarkerEdgeColor','black','HandleVisibility', 'off');
            hold on;
            SP_bypassed_glitches_idxs(index) = [];
            SP_bypassed_glitches_snrs(index) = [];
            SP_bypassed_glitches_snrdiffs(index) = [];
            SP_bypassed_glitches_toadiffs(index) = [];
        end
    end
    for i = 1:length(type4_anomaly_glitch_idxs)
        if ismember(type4_anomaly_glitch_idxs(i),SP_bypassed_glitches_idxs)
            index = find(SP_bypassed_glitches_idxs == type4_anomaly_glitch_idxs(i));
            scatter(SP_bypassed_glitches_snrdiffs(index), SP_bypassed_glitches_toadiffs(index),sz2,'*','MarkerFaceColor', 'red','MarkerEdgeColor', 'red', 'LineWidth',1.5,'HandleVisibility', 'off');
            hold on;
            SP_bypassed_glitches_idxs(index) = [];
            SP_bypassed_glitches_snrs(index) = [];
            SP_bypassed_glitches_snrdiffs(index) = [];
            SP_bypassed_glitches_toadiffs(index) = [];
        end
    end
    scatter(SP_bypassed_glitches_snrdiffs, SP_bypassed_glitches_toadiffs,sz,'red','filled','MarkerEdgeColor','black','DisplayName','$S_P$ VETO BYPASSED GLITCHES');
    hold on;
    scatter(lowmass_snrdiffs(1:66), lowmass_toadiffs(1:66),sz,'yellow','filled','square','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(lowmass_snrdiffs(67:end-20), lowmass_toadiffs(67:end-20),sz,'cyan','filled','square','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(lowmass_snrdiffs(end-19:end), lowmass_toadiffs(end-19:end),sz,'blue','filled','square','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(intermmass_snrdiffs(1:66), intermmass_toadiffs(1:66),sz,'yellow','filled','^','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(intermmass_snrdiffs(67:end-20), intermmass_toadiffs(67:end-20),sz,'cyan','filled','^','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(intermmass_snrdiffs(end-19:end), intermmass_toadiffs(end-19:end),sz,'blue','filled','^','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(highmass_snrdiffs(1:66), highmass_toadiffs(1:66),sz,'yellow','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(highmass_snrdiffs(67:end-20), highmass_toadiffs(67:end-20),sz,'cyan','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(highmass_snrdiffs(end-19:end), highmass_toadiffs(end-19:end),sz,'blue','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    %Added Very High Mass injections
%     for i = 1:length(vhm_cidxs)
%         if ismember(vhm_cidxs(i),veryhighmasssigsegments(1:66))
%             scatter(vhm_snrdiffs(i), vhm_toadiffs(i),sz,'yellow','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%         end
%         if ismember(vhm_cidxs(i),veryhighmasssigsegments(67:132))
%             scatter(vhm_snrdiffs(i), vhm_toadiffs(i),sz,'cyan','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%         end
%         if ismember(vhm_cidxs(i),veryhighmasssigsegments(133:end))
%             scatter(vhm_snrdiffs(i), vhm_toadiffs(i),sz,'blue','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%         end
%     end
    %Add [25,40] M_sun injections in realistic SNRs
    scatter(massgap_snrdiffs, massgap_toadiffs,sz,'black','+','MarkerEdgeColor','black','LineWidth',2,'HandleVisibility','off');
    hold on;
    scatter(Msun_snrdiffs, Msun_toadiffs,sz2,'yellow','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(Msun_snrdiffs_low, Msun_toadiffs_low,sz2,'cyan','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(Msun_snrdiffs_high, Msun_toadiffs_high,sz2,'blue','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    patch(x,y,'green','FaceAlpha',0.04,'HandleVisibility','off');
    hold on;
%     scatter(nan,nan,sz2,'^','MarkerEdgeColor','black','DisplayName','INTERMEDIATE MASS INJECTIONS'); 
    scatter(nan,nan,sz2,'square','MarkerEdgeColor','black','DisplayName','$m_1,m_2 \in [1.4,3] M_{\odot}$');
    hold on; 
    scatter(nan,nan,sz2,'^','MarkerEdgeColor','black','DisplayName','$m_1 \in [1.4,3] M_{\odot}$ and $m_2 \in [5,10] M_{\odot}$');
    hold on;
    scatter(nan,nan,sz2,'diamond','MarkerEdgeColor','black','DisplayName','$m_1,m_2 \in [10,25] M_{\odot}$');
    hold on;
    scatter(nan,nan,sz2,'hexagram','MarkerEdgeColor','black','DisplayName','$m_1,m_2 \in [25,40] M_{\odot}$');
    hold on;
    scatter(nan,nan,sz2,'+','MarkerEdgeColor','black','LineWidth',2,'DisplayName','Mass-gap injections');
    hold on;
    scatter(nan,nan,sz2,'pentagram','MarkerEdgeColor','black','DisplayName','TYPE 1'); 
    hold on;
    scatter(nan,nan,sz2,'v','MarkerEdgeColor','black','DisplayName','TYPE 2');
    hold on; 
    scatter(nan,nan,sz2,'<','MarkerEdgeColor','black','DisplayName','TYPE 3');
    hold on;
    scatter(nan,nan,sz2,'*','MarkerEdgeColor','black','DisplayName','TYPE 4');
    hold on;
    scatter(nan,nan,sz3,"cyan",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,13]$');
    hold on; 
    scatter(nan,nan,sz3,"yellow",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,40]$'); 
    hold on;
%     scatter(nan,nan,sz3,"black",'filled','MarkerEdgeColor','black','DisplayName','[40,80] and [80,110] $M_{sun}$ injections');
%     hold on;
    scatter(nan,nan,sz3,"blue",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [100,500]$');
    hold off;
    legend;
    grid on;
    xlabel('$|\Delta \rho|$','Interpreter','latex');
    ylabel('$|\Delta t_a|$','Interpreter','latex');
    ax = gca;
    ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
    [leg, icons] = legend('FontSize',25,'Interpreter','latex');
    h = findobj(icons,'type','patch');
    set(h,'MarkerSize',15);

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vetoed glitch plot for detection threshold = 9 with very high mass
%injections
% if flag
%     sz = 100;
%     sz2 = 200;
%     sz3 = sz2;
% %     x = [deltasnr_threshold,1,1,deltasnr_threshold];
% %     y = [deltatoa_threshold, deltatoa_threshold, 1000, 1000];
%     figure;
%    scatter(glitch_snrdiffs, glitch_toadiffs,sz,'MarkerFaceColor', [224, 224, 224]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(lowmass_snrdiffs(1:66), lowmass_toadiffs(1:66),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(lowmass_snrdiffs(67:end-20), lowmass_toadiffs(67:end-20),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(lowmass_snrdiffs(end-19:end), lowmass_toadiffs(end-19:end),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(intermmass_snrdiffs(1:66), intermmass_toadiffs(1:66),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(intermmass_snrdiffs(67:end-20), intermmass_toadiffs(67:end-20),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(intermmass_snrdiffs(end-19:end), intermmass_toadiffs(end-19:end),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(highmass_snrdiffs(1:66), highmass_toadiffs(1:66),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(highmass_snrdiffs(67:end-20), highmass_toadiffs(67:end-20),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(highmass_snrdiffs(end-19:end), highmass_toadiffs(end-19:end),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     %Added Very High Mass injections
%     for i = 1:length(vhm_cidxs)
%         if ismember(vhm_cidxs(i),veryhighmasssigsegments(1:66))
%             scatter(vhm_snrdiffs(i), vhm_toadiffs(i),sz,'yellow','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%         end
%         if ismember(vhm_cidxs(i),veryhighmasssigsegments(67:132))
%             scatter(vhm_snrdiffs(i), vhm_toadiffs(i),sz,'cyan','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%         end
%         if ismember(vhm_cidxs(i),veryhighmasssigsegments(133:end))
%             scatter(vhm_snrdiffs(i), vhm_toadiffs(i),sz,'blue','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%         end
%     end
%     %Add [25,40] M_sun injections in realistic SNRs
% %     scatter(massgap_snrdiffs, massgap_toadiffs,sz,'MarkerFaceColor', [224, 224, 224]./255','MarkerEdgeColor','black','HandleVisibility','off');
% %     hold on;
% %     scatter(Msun_snrdiffs, Msun_toadiffs,sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
% %     hold on;
% %     scatter(Msun_snrdiffs_low, Msun_toadiffs_low,sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
% %     hold on;
% %     scatter(Msun_snrdiffs_high, Msun_toadiffs_high,sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
% %     hold on;
% %     patch(x,y,'green','FaceAlpha',0.04,'HandleVisibility','off');
% %     hold on;
% %     scatter(nan,nan,sz2,'^','MarkerEdgeColor','black','DisplayName','INTERMEDIATE MASS INJECTIONS'); 
% %     scatter(nan,nan,sz2,'square','MarkerEdgeColor','black','DisplayName','$m_1,m_2 \in [1.4,3] M_{\odot}$');
% %     hold on; 
% %     scatter(nan,nan,sz2,'^','MarkerEdgeColor','black','DisplayName','$m_1 \in [1.4,3] M_{\odot}$ and $m_2 \in [5,10] M_{\odot}$');
% %     hold on;
%     scatter(nan,nan,sz2,'diamond','MarkerEdgeColor','black','DisplayName','$m_1 \in [40,80] M_{\odot}$ and $m_2 \in [80,110] M_{\odot}$');
%     hold on;
% %     scatter(nan,nan,sz2,'black','filled','MarkerEdgeColor','black','DisplayName','Mass-gap Injections');
% %     hold on;
% %     scatter(nan,nan,sz2,'pentagram','MarkerEdgeColor','black','DisplayName','TYPE 1'); 
% %     hold on;
% %     scatter(nan,nan,sz2,'v','MarkerEdgeColor','black','DisplayName','TYPE 2');
% %     hold on; 
% %     scatter(nan,nan,sz2,'<','MarkerEdgeColor','black','DisplayName','TYPE 3');
% %     hold on;
% %     scatter(nan,nan,sz2,'*','MarkerEdgeColor','black','DisplayName','TYPE 4');
% %     hold on;
%     scatter(nan,nan,sz3,"cyan",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,13]$');
%     hold on; 
%     scatter(nan,nan,sz3,"yellow",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,40]$'); 
%     hold on;
%     scatter(nan,nan,sz3,"blue",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [100,500]$');
%     hold off;
%     legend;
%     grid on;
%     xlabel('$|\Delta \rho|$','Interpreter','latex');
%     ylabel('$|\Delta t_a|$','Interpreter','latex');
%     ax = gca;
%     ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
%     [leg, icons] = legend('FontSize',45,'Interpreter','latex');
%     h = findobj(icons,'type','patch');
%     set(h,'MarkerSize',15);
% 
% 
% 
% 
%     %%%%%PyCBC Injections Plot%%%%%%%%%%%%%%%%%%%%%
%     figure;
%     scatter(chirplength_vetoed_glitches_snrdiffs, chirplength_vetoed_glitches_toadiffs,sz,'MarkerFaceColor', [255, 204, 204]./255,'MarkerEdgeColor','black','DisplayName','CHIRPLENGTH VETOED GLITCHES');
%     hold on;
%     scatter(complexmass_vetoed_glitches_snrdiffs, complexmass_vetoed_glitches_toadiffs,sz,'MarkerFaceColor', [255,102,255]./255,'MarkerEdgeColor','black','DisplayName','COMPLEX MASS VETOED GLITCHES');
%     hold on;
%     scatter(SP_bypassed_glitches_snrdiffs, SP_bypassed_glitches_toadiffs,sz,'red','filled','MarkerEdgeColor','black','DisplayName','$S_P$ VETO BYPASSED GLITCHES');
%     hold on;
%     scatter(lowmass_snrdiffs(1:66), lowmass_toadiffs(1:66),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(lowmass_snrdiffs(67:end-20), lowmass_toadiffs(67:end-20),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(lowmass_snrdiffs(end-19:end), lowmass_toadiffs(end-19:end),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(intermmass_snrdiffs(1:66), intermmass_toadiffs(1:66),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(intermmass_snrdiffs(67:end-20), intermmass_toadiffs(67:end-20),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(intermmass_snrdiffs(end-19:end), intermmass_toadiffs(end-19:end),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(highmass_snrdiffs(1:66), highmass_toadiffs(1:66),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(highmass_snrdiffs(67:end-20), highmass_toadiffs(67:end-20),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(highmass_snrdiffs(end-19:end), highmass_toadiffs(end-19:end),sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     %Added PyCBC Lower Mass injections
%     for i = 1:length(pycbc_lower_cidxs)
%         if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(1:60)) 
%             if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(1:20))
%                 scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'cyan','filled','square','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%             end
%             if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(21:40))
%                 scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'yellow','filled','square','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%             end
%             if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(41:60))
%                 scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'blue','filled','square','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%             end
%         end
%         if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(61:120)) 
%             if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(61:80))
%                 scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'cyan','filled','^','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%             end
%             if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(81:100))
%                 scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'yellow','filled','^','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%             end
%             if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(101:120))
%                 scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'blue','filled','^','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%             end
%         end
%         if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(121:180)) 
%             if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(121:140))
%                 scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'cyan','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%             end
%             if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(141:160))
%                 scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'yellow','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%             end
%             if ismember(pycbc_lower_cidxs(i),pycbclowersigsegments(161:180))
%                 scatter(pycbc_lower_snrdiffs(i), pycbc_lower_toadiffs(i),sz,'blue','filled','diamond','MarkerEdgeColor','black','HandleVisibility','off');
%             hold on;
%             end
%         end
%     end
%     %Added PyCBC High Mass injections
%     % for i = 1:length(pycbc_cidxs)
%     %     if ismember(pycbc_cidxs(i),pycbcsigsegments(21:40)) || ismember(pycbc_cidxs(i),pycbcsigsegments(81:100))
%     %         scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz,'yellow','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
%     %         hold on;
%     %     end
%     %     if ismember(pycbc_cidxs(i),pycbcsigsegments(1:20)) || ismember(pycbc_cidxs(i),pycbcsigsegments(61:80))
%     %         scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz,'cyan','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
%     %         hold on;
%     %     end
%     %     if ismember(pycbc_cidxs(i),pycbcsigsegments(41:60)) || ismember(pycbc_cidxs(i),pycbcsigsegments(101:120))
%     %         scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz,'blue','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
%     %         hold on;
%     %     end
%     % end
% 
%     for i = 1:length(pycbc_cidxs)
%         if ismember(pycbc_cidxs(i),pycbcsigsegments(1:60))
%             if ismember(pycbc_cidxs(i),pycbcsigsegments(1:20))
%                 scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'cyan','filled','pentagram','MarkerEdgeColor','black','HandleVisibility','off');
%                 hold on;
%             end
%             if ismember(pycbc_cidxs(i),pycbcsigsegments(21:40))
%                 scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'yellow','filled','pentagram','MarkerEdgeColor','black','HandleVisibility','off');
%                 hold on;
%             end
%             if ismember(pycbc_cidxs(i),pycbcsigsegments(41:60))
%                 scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'blue','filled','pentagram','MarkerEdgeColor','black','HandleVisibility','off');
%                 hold on;
%             end
%         end
%         if ismember(pycbc_cidxs(i),pycbcsigsegments(61:120))
%             if ismember(pycbc_cidxs(i),pycbcsigsegments(61:80))
%                 scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'cyan','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
%                 hold on;
%             end
%             if ismember(pycbc_cidxs(i),pycbcsigsegments(81:100))
%                 scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'yellow','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
%                 hold on;
%             end
%             if ismember(pycbc_cidxs(i),pycbcsigsegments(101:120))
%                 scatter(pycbc_snrdiffs(i), pycbc_toadiffs(i),sz2,'blue','filled','hexagram','MarkerEdgeColor','black','HandleVisibility','off');
%                 hold on;
%             end
%         end
%     end
% 
%     %Add [25,40] M_sun injections in realistic SNRs
% %     scatter(massgap_snrdiffs, massgap_toadiffs,sz,'MarkerFaceColor', [224, 224, 224]./255','MarkerEdgeColor','black','HandleVisibility','off');
% %     hold on;
% %     scatter(Msun_snrdiffs, Msun_toadiffs,sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
% %     hold on;
% %     scatter(Msun_snrdiffs_low, Msun_toadiffs_low,sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
% %     hold on;
% %     scatter(Msun_snrdiffs_high, Msun_toadiffs_high,sz,'MarkerFaceColor', [160,160,160]./255,'MarkerEdgeColor','black','HandleVisibility','off');
% %     hold on;
%     patch(x,y,'green','FaceAlpha',0.04,'HandleVisibility','off');
%     hold on;
% %     scatter(nan,nan,sz2,'^','MarkerEdgeColor','black','DisplayName','INTERMEDIATE MASS INJECTIONS'); 
%     scatter(nan,nan,sz2,'square','MarkerEdgeColor','black','DisplayName','$m_1,m_2 \in [1.4,3] M_{\odot}$');
%      hold on; 
%      scatter(nan,nan,sz2,'^','MarkerEdgeColor','black','DisplayName','$m_1 \in [1.4,3] M_{\odot}$ and $m_2 \in [5,10] M_{\odot}$');
%      hold on;
%     scatter(nan,nan,sz2,'diamond','MarkerEdgeColor','black','DisplayName','$m_1,m_2 \in [10,25] M_{\odot}$');
%      hold on;
%      scatter(nan,nan,sz2,'pentagram','MarkerEdgeColor','black','DisplayName','$m_1,m_2 \in [25,40] M_{\odot}$');
%      hold on;
%      scatter(nan,nan,sz2,'hexagram','MarkerEdgeColor','black','DisplayName','$m_1 \in [40,80] M_{\odot}$ and $m_2 \in [80,110] M_{\odot}$');
%     hold on;
% 
% %     scatter(nan,nan,sz2,'black','filled','MarkerEdgeColor','black','DisplayName','Mass-gap Injections');
% %     hold on;
% %     scatter(nan,nan,sz2,'pentagram','MarkerEdgeColor','black','DisplayName','TYPE 1'); 
% %     hold on;
% %     scatter(nan,nan,sz2,'v','MarkerEdgeColor','black','DisplayName','TYPE 2');
% %     hold on; 
% %     scatter(nan,nan,sz2,'<','MarkerEdgeColor','black','DisplayName','TYPE 3');
% %     hold on;
% %     scatter(nan,nan,sz2,'*','MarkerEdgeColor','black','DisplayName','TYPE 4');
% %     hold on;
%     scatter(nan,nan,sz3,"cyan",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,13]$');
%     hold on; 
%     scatter(nan,nan,sz3,"yellow",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,40]$'); 
%     hold on;
%     scatter(nan,nan,sz3,"blue",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [100,500]$');
%     hold off;
%     legend;
%     grid on;
%     xlabel('$|\Delta \rho|$','Interpreter','latex');
%     ylabel('$|\Delta t_a|$','Interpreter','latex');
%     xlim([1e-7, 1]);
%     ax = gca;
%     ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
%     [leg, icons] = legend('FontSize',40,'Interpreter','latex');
%     h = findobj(icons,'type','patch');
%     set(h,'MarkerSize',15);
% end
% end



%Vetoed Glitch Plot for detection threshold = 8.0
if detection_threshold == 8
    sz = 100;
    sz2 = 200;
    sz3 = 200;
    x = [deltasnr_threshold,1,1,deltasnr_threshold];
    y = [deltatoa_threshold, deltatoa_threshold, 1000, 1000];
    figure;
    hold on;
    scatter(noisesnrdiffs, noisetoadiffs,sz, 'MarkerFaceColor', [224, 224, 224]./255,'MarkerEdgeColor','black','DisplayName','$\widehat{\rho}_P < 8.0$ EVENTS');
    hold on;
    scatter(abssimnoiseSNRdiff, abssimnoiseTOAdiff,sz, 'MarkerFaceColor', [224, 224, 224]./255,'MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(abssimnoiseSNRdiff_2, abssimnoiseTOAdiff_2,sz, 'MarkerFaceColor', [224, 224, 224]./255,'MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(chirplength_vetoed_glitches_snrdiffs, chirplength_vetoed_glitches_toadiffs,sz,'MarkerFaceColor', [255, 204, 204]./255,'MarkerEdgeColor','black','DisplayName','CHIRPLENGTH VETOED GLITCHES');
    hold on;
    scatter(complexmass_vetoed_glitches_snrdiffs, complexmass_vetoed_glitches_toadiffs,sz,'MarkerFaceColor', [255,102,255]./255,'MarkerEdgeColor','black','DisplayName','COMPLEX MASS VETOED GLITCHES');
    hold on;
    scatter(SP_bypassed_glitches_snrdiffs, SP_bypassed_glitches_toadiffs,sz,'red','filled','MarkerEdgeColor','black','DisplayName','$S_P$ VETO BYPASSED GLITCHES');
    hold on;
    hold on;
    scatter(realsnr_snrdiffs, realsnr_toadiffs,sz,'yellow','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(lowsnr_snrdiffs, lowsnr_toadiffs,sz,'cyan','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
    scatter(highsnr_snrdiffs, highsnr_toadiffs,sz,'blue','filled','MarkerEdgeColor','black','HandleVisibility','off');
    hold on;
%     scatter(massgap_snrdiffs, massgap_toadiffs,sz,'black','filled','MarkerEdgeColor','black','DisplayName','Mass-gap Injections');
%     hold on;
%     scatter(Msun_snrdiffs, Msun_toadiffs,sz,'yellow','filled','MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(Msun_snrdiffs_low, Msun_toadiffs_low,sz,'cyan','filled','MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
%     scatter(Msun_snrdiffs_high, Msun_toadiffs_high,sz,'blue','filled','MarkerEdgeColor','black','HandleVisibility','off');
%     hold on;
    patch(x,y,'green','FaceAlpha',0.04,'HandleVisibility','off');
    hold on;
    scatter(nan,nan,sz3,"cyan",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,13]$');
    hold on; 
    scatter(nan,nan,sz3,"yellow",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [10,40]$'); 
    hold on;
    scatter(nan,nan,sz3,"blue",'filled','MarkerEdgeColor','black','DisplayName','SNR $\in [100,500]$');
    hold off;
    legend;
    grid on;
    xlabel('$|\Delta \rho|$','Interpreter','latex');
    ylabel('$|\Delta t_a|$','Interpreter','latex');
    ax = gca;
    ax.XAxis.FontSize = 40; ax.YAxis.FontSize = 40;
    [leg, icons] = legend('FontSize',25,'Interpreter','latex');
    h = findobj(icons,'type','patch');
    set(h,'MarkerSize',15);
end













