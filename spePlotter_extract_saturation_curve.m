%user parameters
center_wl = 637; %in nm
spec_size = 1340; %spectrometer ccd size
%wl =0.0195.*(1:1340) + center_wl-0.0195*1340/2; %1200g grating (new) Calibrated jul 08 2020
wl =0.0111.*(1:1340) + center_wl-0.0111*1340/2; %1800g grating (new) Calibrated jul 08 2020
%wl =0.0872.*(1:1340) + center_wl-0.0872*1340/2; %300g grating (new) Calibrated jul 15 2020

%wl =1:spec_size; %for debug

%Select the files
[files,dir] = uigetfile('D:\Dropbox\Projects\Angled Implant\PLE_data_mar_2021\round-2\ppb1\spectra\nv_103642\charge state check\*637*.spe', 'MultiSelect','on'); 
if ischar(files)
    files = {files}; 
end

%Initialize figure
h=figure;
cts = zeros(1,length(files));

for n=1:length(files)
    
    %Get data
    fname = [dir cell2mat(files(n))];
    [image,time] = readSPE(fname); 
    data = double(image);

    %figure(h)
    plot(1:size(data,2),data);
    
    %User selection of Peak to fit
    set(gcf,'CurrentCharacter',char(1));
    h=datacursormode;
    set(h,'DisplayStyle','datatip');
    waitfor(gcf,'CurrentCharacter',char(32));
    s = getCursorInfo(h);
    start_x = s.Position(1);
    start_val = s.Position(2);
        
%     %Data weights for peak detection 
%     if center_wl==637
%         wt = [linspace(0.001,1,center_wl) linspace(1,0.001,1340-center_wl)];%for NV-
%     else
%         wt = [0*linspace(0.01,0.1,center_wl-20) linspace(0.1,1,20) linspace(1,0.001,1340-center_wl)];%For NV0
%     end
%     
% %     %Custom ROI for nosiy data
% %     if n ==1 && center_wl==637 %long exposure; cosmic rays; manual peak selection
% %         p=674;
% %     else
% %         [~, p] = max(data.*wt);
% %     end
    
    %[~, p] = max(data.*wt);
    p = start_x;
    range = p-6:p+6; %integration range
    data = data-min(data(range));
    data1 = data/time;

    %Integrate area under curve
    integ = trapz((data1(range)));
    cts(n) = integ;

%     %Bkg subtraction
%         bkg = mean(data(end-50:end));
%         data1=(data-bkg)/time;
%         
%         if min(data1)<0
%             data1 = data1 + abs(min(data1));
%         end

    %Make plot for verification
    plot(wl, (data1)); hold on
    scatter(wl(range(1)),0); scatter(wl(p),data1(p)); scatter(wl(range(end)),0);
    xlabel('Wavelength (nm)'); 
    ylabel('cts/s'); 
    %title(fname)
    
    %Verify the peak
    currkey=0;
    % do not move on until enter key is pressed
    while currkey~=1
        pause; % wait for a keypress
        currkey=get(gcf,'CurrentKey'); 
        if strcmp(currkey, 'return') 
            currkey=1; % 
        else
            currkey=0; %
        end
    end

    hold off
end

% xlim([600 750]);
% ylim([0 max((data1/time)/1000)]);

%for publication data  
% set(gca, 'fontname','Palatino Linotype')
% set(gca, 'FontSize',20)
% set(gcf, 'Color', 'None');
% set(gca, 'Color', 'None');

%plot saturation curve
ex = [25 50 75 150 250 500 1000]*0.6;  
wt = linspace(1, 0.5, length(ex))';

%fit saturation curve
xData = ex';
%yData = cts';
yData = nvn./nvp;
ft = fittype( 'a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0];
opts.StartPoint = [200 0.5];
%opts.Weights = wt;

% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft, opts);

% Plot fit with data.
figure
plot( fitresult, xData, yData );
fitresult.a
fitresult.b

ylim([0 4])
set(gca, 'fontname','Palatino Linotype')
set(gca, 'FontSize',20)