%Change the path below to the required folder
[files,dir] = uigetfile('D:\Dropbox\Projects\*1800*.spe', 'MultiSelect','on'); 
if ischar(files)
    files = {files}; 
end

%-------------------------------------------------------------------------
% %User options
stitched_spectra = 0;  %is this a winspec stitched spectra (i.e pixels > ccd_width ?)
wl_calibrated = 1;    %is the winspec x-axis data calibrated ?
new_spectrometer = 1; %is this the new black (HP) spectrometer or the old blue box
try_fit = 1;          %try to fit the peak to a lorentzian (for cavity Q-analysis)
bkg_subtract = 0;     %try a crude bkg subtration (mean of the first 100 pixels)
integrate_curve = 0;  %integrate to find area under the curve?
save_plot_png = 0;    %save a png file?
%-------------------------------------------------------------------------


% set spectra axis calibration parameters
if ~stitched_spectra && ~wl_calibrated
    
    center_wl = 637; %in nm; maunual configuration of center wl
    
    if ~new_spectrometer %this is for the old blue box spectrometer (now on brynn mc)
        wl = 0.706.*(1:512) + center_wl-0.706*512/2; %50g grating
        full_data = zeros(length(files),1,512); 

    else
        %wl =0.0195.*(1:1340) + center_wl-0.0195*1340/2; %1200g grating (new) Calibrated jul 08 2020
        %wl =0.0111.*(1:1340) + center_wl-0.0111*1340/2; %1800g grating (new) Calibrated jul 08 2020
        wl =0.0872.*(1:1340) + center_wl-0.0872*1340/2; %300g grating (new) Calibrated jul 15 2020
        full_data = zeros(length(files),1,1340);
    end
end

% Iterate through the files
for n=1:length(files)
    
    figure; %make a figure for the final plot 
    
    fname = [dir cell2mat(files(n))];
    [image] = loadSPE(fname);  % extract the raw data from spe
    
    % If winspec calibration exists use it, else manual wl calibration
    if wl_calibrated
        data_x = image.wavelength;
    else
        data_x = wl;
    end
    
    data_y = image.int'/image.expo_time; % Normalize to cts/s
    
    % %Bkg subtraction
    if bkg_subtract
        bkg = mean(data(1:100));
        data_y = data_y - bkg;
    end
        
    %Make plot    
    plot(data_x,data_y)
    xlabel('Wavelength (nm)'); 
    ylabel('(cts/s)'); 
        
    if try_fit
        
        %Fit the spectra to a Lorentzian
        dx = (data_x(2)-data_x(1));
        
        % %Assume peak to fit is @ the max intensity ?
        %[m,m_ind] = max(data_y);   
        %start_x = data_x(m_ind);
        
        %User selection of Peak to fit
        set(gcf,'CurrentCharacter',char(1));
        h=datacursormode;
        set(h,'DisplayStyle','datatip');
        waitfor(gcf,'CurrentCharacter',char(32));
        s = getCursorInfo(h);
        start_x = s.Position(1);
        start_val = s.Position(2);
        start_ind = find(data_x == start_x);
        
        %Define a lorentzian fit
        lorentz = fittype('e1*x + d1 + a1*c1./(2*pi*( ((x-b1).^2) + ((c1/2)^2) ))', ...
                                'dependent',{'y'},'independent',{'x'}, ...
                                'coefficients',{'a1','b1','c1','d1','e1'});
         
        %Try the fit function                                                                                        
       % try
            f1 = fit(data_x(start_ind-100:start_ind+100), data_y(start_ind-100:start_ind+100), ...
                lorentz, 'StartPoint',[100, start_x, 10*dx, 500,1], ...
                'Upper', [start_val, start_x + 10*dx, 50*dx, 2*start_val,1000], ...
                'Lower', [1, start_x - 30*dx, 2*dx, 0,-1000]);
            
            %Update plot and print results
            hold on;
            plot(f1);
            title(['Fit of scan ' num2str(n) ' width=' num2str(f1.c1) ' Q=' num2str(f1.b1/f1.c1)]);
          
            fprintf('FWHM = %2.3f \t Q = %4.1f \n',f1.c1, f1.b1/f1.c1)        
       % end
        
    else
        
        title(strrep(files{n}(1:end-4), '_',' ')); 
        set(gca, 'fontname','Palatino Linotype')
        set(gca, 'FontSize',15);
        axis tight
        
    end

    if integrate_curve
        range = start_x-5:start_x+5; % Make this smarter; not hardcoded
        integ = trapz(abs(data_y(range)))
    end

    if save_plot_png
        saveas(h,[files{n}(1:end-4) '.png']); 
    end
end


%  xlim([600 750]);
% ylim([0 max((data1/time)/1000)]);