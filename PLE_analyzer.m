%Code to fit NV PLE data

%Date: 07/05/2020  (this is up-to-date)
%Srivatsa Chakravarthi

%This code can process multiple PLE data files generated either by manual
%PLE scans or automated PLE and extract statistics for the NV linewidth

%When fitting single scans, pressing enter will accept the fit and
%pressing backspace will reject the fit (other buttons do nothing). 
%Only accepted fits will be used for data analysis and final plot

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
close all
clear all

% Annoying matlab
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')

dir = 'D:\Dropbox\Projects\Angled Implant\PLE_data_mar_2021\round-2\';

% User parameters (%1=Yes; 0=No;)
autople = 0;    % Are the files saved by autople code?
make_fig = 1;   % Show raw data with frequency axis
show_stats = 1; % Show the scan parameters

processLines = 0; % Process data (1=Yes; 0=No)

    show_fits = 1;    % Show individual fits for verification (1= Manual; 0=Automated)
    gauss_fit = 0;    % Use Gaussian curve fitting
    lorentz_fit = 1;  % Use Lorentzian curve fitting
    min_height = 6;   % Threshold for fiting (in photons/bin)
    pix2fit = 20;     % No. of pixels to fit
    
    first_scan = 1;  % First scan # (to ignore initial laser variation)
        % Data range for processing [Truncate scan range; useful for picking specific SB]
        %last = upPixels-100;
        begin = 30;
       % last = 90;
            
    %Final plot    
    pretty_plot = 1;        % Make final plot with only fitted data
        center_pixel = 193; % Where is the PLE trace?
        plot_range = 13;     % (in GHz);
        
make_hist=0; % Make histogram for autople data (1=Yes; 0=No)
start_wl = 637; % wavelength in nm, not critical

% Piezo volt to frequency Calibration
%lambda_per_volt = 0.02;   % Old calibration used by Emma (2019)
%lambda_per_volt = 0.02348; % Updated with wavemeter measurement (March 2020)

lambda_per_volt = 0.014273;


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

% Get the files
[fle, dirt] = uigetfile([dir '*.mat'], 'MultiSelect','on'); 

if iscell(fle)
    numFile = size(fle,2); 
else
    fle = {fle}; 
    numFile = 1; 
end

% Initialize variables
min_width = zeros(1,numFile);
avg_width = zeros(1,numFile);
minwidths=[]; avgwidths=[]; scanstot = []; specdifstds = [];

%fprintf('\n Min width \t Avg width \t Max width \t Std. diff \t Scans \n')

% Iterate through all files
for n = 1:numFile
    
    fname = fle{n}; 
    [tmp, saveFname, ext] = fileparts(fname); 
    data = load([dirt '/' fname]); 
    
    if autople
        data.data = data.data_d;
    end
    
    % Extract variables from data structure
    plotData = cell2mat(data.data.data); 
    titleStr = strrep(saveFname, '_', ' '); 
    scanAxis = data.data.scans{1}; 
    inputs = data.data.inputs{1};    
    voltageAxis = inputs.xaxis; 
    upPixels = inputs.upPixels; 
    dnPixels = length(voltageAxis)-upPixels; 
    yaxis = voltageAxis(1:upPixels);
    
    % Convert y-axis to wl(nm) or GHz
    lambdaAxis = [linspace(start_wl+lambda_per_volt*voltageAxis(1), start_wl+lambda_per_volt*voltageAxis(upPixels), upPixels)]; 
    freqAxis = (2.998e8./(lambdaAxis.*1e-9))./1e9;
    %freqAxis = freqAxis - freqAxis(round(upPixels/2)); 
    freqAxis = abs(freqAxis - freqAxis(1));
    
    last = upPixels;
    
%     bkSknLmbda = linspace(637.202+0.02*voltageAxis(numPixels), 637.202+0.02*voltageAxis(1), szDown); 
%     bkSknFreq = (2.998e8./(bkSknLmbda.*1e-9))./1e9; 

    % Make raw figure (Pl vs Freq)
    if make_fig==1
        h=figure;
        surf(scanAxis, freqAxis, plotData(:,1:upPixels)'); 
        view([0 90]); shading flat; axis tight; colormap bone; colorbar; 
        title(titleStr); 
        xlabel('Scan #'); 
        ylabel('Frequency (GHz)'); 
        caxis([mean(plotData(1,upPixels-10:upPixels)) max(max(plotData(2:end,1:upPixels)))]);
        %savefig(h1, [dirt '\' saveFname '.fig']); 
        %close(h);
        
%         hh=figure;
%         surf(scanAxis, 1:length(plotData(1,:)), plotData(:,:)'); 
%         view([0 90]); shading flat; axis tight; colormap bone; colorbar;
    end
    
    % Get stats
    if show_stats == 1
        max_counts = max(max(plotData(2:end,1:upPixels)))/(data.data.inputs{1,1}.upTime/data.data.inputs{1,1}.upPixels)...
            - mean(plotData(1,upPixels-10:upPixels));
        speed = abs(freqAxis(end)-freqAxis(1))/(data.data.inputs{1,1}.upTime);
        bin_size = abs(freqAxis(2)-freqAxis(1))*1e3; 
        fom = max_counts/bin_size; % Arb. figure of merit

        fprintf('\nFile = %s \n\n Max count rate = %0.2f cts/s \n Scan rate = %0.2f GHz/s \n Bin = %0.2f Mhz\n FOM = %0.2f \n', saveFname, max_counts, speed, bin_size, fom)
    end
    
    % Initialize variables    
    final_data = zeros(size(plotData,1),upPixels+10); %Ten extra pixels at end indicating repump
    final_count=1;

    if processLines  
%         fprintf('\n Min width \t Avg width \t Max width \t Std. diff \t Scans \n')
        
        %intialize variables       
        widths = NaN*ones(size(plotData,1), 1);
        peak = NaN*ones(size(plotData,1), 1);
        peak_in_scans = 0;
        bkscan_old = 0;
        
        % Iterate through each scans
        for m=first_scan:size(plotData,1)
        % for m=first_scan:30  
            %----------------------------------------

            %-----------------------------------------
            
            lineData = plotData(m, begin:last); % Truncated upscan data

            upData = plotData(m,1:upPixels); % Upscan data
            bkscan = plotData(m, upPixels+1:end); % Backscan data
            
            bkscan_new = mean(bkscan(~isnan(bkscan)));
            
            [val, ind] = max(lineData); %This assumes peak intensity is NV PLE line

            min_range = ind-pix2fit;
            max_range = ind+pix2fit;
            
            %Check data range
            if min_range <1
                min_range=1;
                max_range = 2*pix2fit;
            elseif max_range>last-(begin-1);
                max_range = last-(begin-1);
                min_range = max_range - 2*pix2fit;
            end

            %Extract data for fit
            x2fit = freqAxis(min_range:max_range);
            y2fit = lineData(min_range:max_range);
            
            %Try to fit the data
            try
                %Define curve expression and find fit
                if gauss_fit==1
                    gauss_a=fittype('d1+a1*exp(-((x-b1)/c1)^2)','dependent',{'y'},'independent',{'x'},...
                                        'coefficients',{'a1','b1','c1','d1'});
                    f1 = fit(x2fit', y2fit', gauss_a, 'StartPoint',[10,  freqAxis(ind), bin_size/1000, 0]);

                elseif lorentz_fit==1
                    lorentz_a = fittype('d1 + a1*c1./(2*pi*( ((x-b1).^2) + ((c1/2)^2) ))','dependent',{'y'},'independent',{'x'},...
                                        'coefficients',{'a1','b1','c1','d1'});
                    f1 = fit(x2fit', y2fit', lorentz_a, 'StartPoint',[10,  freqAxis(ind), bin_size/1000, 0]);
                end

                %Check quality of fit
                if (max(y2fit)>min_height)
                    
                    if show_fits==1
                        q= figure;
                        hold all;
                        plot(x2fit, y2fit,'-o');
                        plot(f1);
                        
                        if gauss_fit==1
                            title(['Fit of scan ' num2str(m) ' width= ' num2str(2.355*f1.c1/sqrt(2))]);
                        elseif lorentz_fit==1
                            title(['Fit of scan ' num2str(m) ' width= ' num2str(f1.c1)]);
                        end

                        %Verify the peak
                        currkey=0;
                        % do not move on until enter key is pressed
                        while currkey==0
                            waitforbuttonpress % wait for a keypress
                            currkey=get(gcf,'CurrentKey');

                            if strcmp(currkey, 'return') 
                                currkey=1; % Accept the fit
                                
                                if gauss_fit==1
                                    widths(m) = 2.355*f1.c1/(sqrt(2)); %for Gaussian fit
                                elseif lorentz_fit==1
                                    widths(m) = f1.c1; %for Lorentzian fit
                                end
                                
                                peak_in_scans = peak_in_scans +1;
                                peak(m) = f1.b1;

                                %save the scan data for final plot
                                final_data(final_count,1:end-10)=upData; %upscan data
                                final_data(final_count,end-9:end)=bkscan_old; %backscan data
                                final_count = final_count + 1;

                            elseif strcmp(currkey, 'backspace')
                                currkey=2; % Reject the fit
                                widths(m) = NaN; 
                                peak(m) = NaN;

                            else
                                currkey=0; % Wait
                            end

                        end

                            %pause(0.5)
                            close(q); 

                    else
                            if gauss_fit==1
                                if 2.355*f1.c1/(sqrt(2))>bin_size/1000 && f1.c1<50*bin_size/1000 && f1.a1>min_height
                                    widths(m) = 2.355*f1.c1/(sqrt(2)); %for Gaussian fit
                                    peak_in_scans = peak_in_scans +1;
                                    peak(m) = f1.b1;

                                    %save the scan data for final plot
                                    final_data(final_count,1:end-10)=upData; %upscan data
                                    final_data(final_count,end-9:end)=bkscan_old; %backscan data
                                    final_count = final_count + 1;
                                else
                                    widths(m) = NaN; 
                                    peak(m) = NaN;
                                end
                                
                            elseif lorentz_fit==1
                                l_amp = f1.d1 + f1.a1*2/(pi*f1.c1); 
                                if f1.c1>bin_size/1000 && f1.c1<10*bin_size/1000 && l_amp>0 && l_amp<2*max(y2fit)
                                    widths(m) = f1.c1; %for Lorentzian fit
                                    peak_in_scans = peak_in_scans +1;
                                    peak(m) = f1.b1;

                                    %save the scan data for final plot
                                    final_data(final_count,1:end-10)=upData; %upscan data
                                    final_data(final_count,end-9:end)=bkscan_old; %backscan data
                                    final_count = final_count + 1;
                                else
                                    widths(m) = NaN; 
                                    peak(m) = NaN;
                                end
                            end
                    end
                else
                    widths(m) = NaN; 
                    peak(m) = NaN;
                end
            end   
            bkscan_old = bkscan_new;
        end
    end
   
    %Close PLE fig
%     if make_fig*displayLines==1 
%        close(h);
%     end
    
    try
        %Average Scan Metrics
        min_width = 0;
        avg_width = 0;
        max_width = 0;
        std_spec_diff = 0;

        widths = widths(~isnan(widths));
        min_width = min(widths);
        avg_width = mean(widths);
        max_width = max(widths);
        %max_spec_diff = max(peak(~isnan(peak)))-min(peak(~isnan(peak)));
        std_spec_diff = std(peak(~isnan(peak))).*sqrt(peak_in_scans-1);
        
        fprintf('\n Min width \t Avg width \t Max width \t Std. diff \t Scans \n')
        fprintf('% 1.3f \t\t %1.3f \t\t %1.3f \t\t %1.3f \t\t %d\n', min_width, avg_width, max_width, std_spec_diff, peak_in_scans)
        %close all; 

        minwidths(end+1) = min_width;
        avgwidths(end+1) = avg_width;
        scanstot(end+1) = peak_in_scans;
        specdifstds(end+1) = std_spec_diff;
    end
    
    %Make pretty_plot
    if pretty_plot && processLines
       new_scan_axis = 2:final_count-1;
       plot_pixels =plot_range/(bin_size/1000);
       plot_ind = center_pixel-round(plot_pixels/2):center_pixel+round(plot_pixels/2);
        
       figure
      %plot with repump
      surf(new_scan_axis, [freqAxis(plot_ind) freqAxis(plot_ind(end))+(1:10)*bin_size/1000]-freqAxis(plot_ind(1)), (final_data(new_scan_axis,[plot_ind length(freqAxis)+(1:10)])./(data.data.inputs{1,1}.upTime/data.data.inputs{1,1}.upPixels))');  
       
      %plot with no repump
      % surf(new_scan_axis, freqAxis(plot_ind)-freqAxis(plot_ind(1)), (final_data(new_scan_axis,plot_ind)./(1e3*data.data.inputs{1,1}.upTime/data.data.inputs{1,1}.upPixels))');  
      
       view([-90 90]); shading flat; axis tight; colormap cool; colorbar; 
       xlabel('Scan #'); 
       ylabel('Frequency (GHz)'); 
       %load('blue_ple_cmap.mat');
       %colormap(cmap);
       colormap('gray')
       %caxis([500 2500]); 
       ax = gca;
       ax.YDir = 'reverse';
        
%         %For publication data only
        set(gca, 'fontname','Palatino Linotype')
        set(gca, 'FontSize',18);
        set(gcf, 'Color', 'None');
        set(gca, 'Color', 'None');

    end
    
end

if make_hist
    %%Make Histograms for autoPLE data
    figure
    histogram(minwidths);
    xlabel('min instantaneous linewidth (GHz)');

    figure
    histogram(avgwidths);
    xlabel('average instantaneous linewidth (GHz)');
    %savepngandfig(lwfilename, savedir);

    % histogram(avgwidths);
    % xlabel('spectral diffusion std (GHz)');
    %savepngandfig(diffusfilename, savedir);
end