% FUll description needed
function propagateGUI()
    close all;
    clear;
    function [rgb] = complex2rgb(A, scale)
       if nargin == 1
          reduce = abs(A);
          ndim = length(size(reduce));
          for i=1:ndim
              reduce = max(reduce);
          end
          scale = reduce;
       end

       H = angle(A)+pi/2;
       H (H<0)=H(H<0)+2*pi;
       H = H/(pi/3);
       S = 1;
       V = abs(A)/scale;

       C = V * S;
       X = C .* (1 - abs(rem(H,2) - 1));
       Z = zeros(size(V));

       R = C;
       G = X;
       B = Z;

       R(H>1) = X(H>1);
       G(H>1) = C(H>1);
       B(H>1) = Z(H>1);

       R(H>2) = Z(H>2);
       G(H>2) = C(H>2);
       B(H>2) = X(H>2);

       R(H>3) = Z(H>3);
       G(H>3) = X(H>3);
       B(H>3) = C(H>3);

       R(H>4) = X(H>4);
       G(H>4) = Z(H>4);
       B(H>4) = C(H>4);

       R(H>5) = C(H>5);
       G(H>5) = Z(H>5);
       B(H>5) = X(H>5);

       m = V - C;
       R = R + m;
       G = G + m;
       B = B + m;
       rgb = reshape([R(:), G(:), B(:)], [size(R), 3]);
       rgb(rgb<0)=0;
       rgb(rgb>1)=1;
    end
        
    function dataTime = getTimeStamp
    %GETTIMESTAMP Summary of this function goes here
    %   Detailed explanation goes here
        dataVector = clock;
        dataTime = strcat(num2str(dataVector(3)),'.',num2str(dataVector(2)),'.',num2str(dataVector(1)),'_',num2str(dataVector(4)),'_',num2str(dataVector(5)));    
    end
    % Create the main window, where everything is happening
    timeStamp = getTimeStamp;
    
    mainF = figure('Name','PROBE PROPAGATOR',...
        'NumberTitle','off','units',...
        'normalized','outerposition',[0 0.2 0.8 0.8]);
    
    % Plot of initial frame at the sample position
    startIm = imread('introLogo.png');
    hIm = imagesc(startIm); axis image; zoom(1); set(gca,'xtick',[],'ytick',[])
        
    uicontrol('Parent',mainF,'Style','text',...
        'String','Dzhigaev Dmitry 2016','Units','normalized',...
        'Position',[0.01 0.05 0.2 0.03]);
    
    leftPanel = uipanel('Title','Plot control','FontSize',...
        12,'Units','Normalized','Position',[.01 0.1 .2 .85]);
    
    rightPanel = uipanel('Title','Propagation control',...
        'FontSize',12,'Units','Normalized','Position',[.78 .1 .2 .85]);
     
    % Button for openning the file 
    uicontrol('Parent',leftPanel,'Units','normalized',...
        'Position',[0.05 0.87 0.9 0.1],'String','Open file',...
        'Callback',@openFile);                
            
    % Main function of the display
    function openFile(hObj,callbackdata)
        % Open a dialog for the MAT file selecting
        filterExtension              = '*.mat';
        dialogTitle                  = 'Select a MAT-file containing probe';
        defaultFileName              = 'results.mat';
        [param.fileName, param.pathName] = uigetfile(filterExtension, dialogTitle, fullfile(pwd,defaultFileName));
        inputFile = fullfile(param.pathName, param.fileName);
        
        % Dynamic determination of the dataset name
        in = load(inputFile); 
        names = fieldnames(in);
        IN = in.(names{1});
        
        clear filterExtension dialogTitle defaultFileName;        
        % Now we have the slice IN in the memory ->
        
        % Clear figure before drawing a new propagation
        clf;
        
        % Plot of initial frame at the sample position
        startIm = complex2rgb(IN);
        hIm = imagesc(startIm); axis image; zoom(1); set(gca,'xtick',[],'ytick',[])

        % Main window creation
        uicontrol('Parent',mainF,'Style','text',...
            'String','Dzhigaev Dmitry 2016',...
            'Units','normalized','Position',[0.01 0.05 0.2 0.03]);
    
        leftPanel = uipanel('Title','Plot control','FontSize',...
            12,'Units','Normalized','Position',[.01 0.1 .2 .85]);

        rightPanel = uipanel('Title','Propagation control',...
            'FontSize',12,'Units','Normalized',...
            'Position',[.78 .1 .2 .85]);

        % Button for openning the file 
        uicontrol('Parent',leftPanel,'Units','normalized',...
            'Position',[0.05 0.87 0.9 0.1],'String','Open file',...
            'Callback',@openFile);        
        
        %##################################################################
        % Right panel with a propagation parameters #######################               
        uicontrol('Parent',rightPanel,'Style','text',...
            'String','Photon energy, [eV]',...
            'Units','normalized','Position',[0.1 0.88 0.4 0.1]);
        
        hEnergy = uicontrol('Parent',rightPanel,'Style','edit',...
            'String','15250',...
            'Units','normalized','Position',[0.5 0.93 0.4 0.05]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uicontrol('Parent',rightPanel,'Style','text',...
            'String','Sample2Detector, [m]',...
            'Units','normalized','Position',[0.1 0.78 0.4 0.1]);
        
        hSampleToDetector = uicontrol('Parent',rightPanel,'Style','edit',...
            'String','2.5',...
            'Units','normalized','Position',[0.5 0.83 0.4 0.05]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uicontrol('Parent',rightPanel,'Style','text',...
            'String','Detector pixel size, [m]',...
            'Units','normalized','Position',[0.1 0.68 0.4 0.1]);
        
        hDetectorPixel = uicontrol('Parent',rightPanel,'Style','edit',...
            'String','176e-6',...
            'Units','normalized','Position',[0.5 0.73 0.4 0.05]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uicontrol('Parent',rightPanel,'Style','text',...
            'String','Propagation Interval, [m]',...
            'Units','normalized','Position',[0.1 0.58 0.4 0.1]);
        
        hIntervalPropagateFrom = uicontrol('Parent',rightPanel,'Style','edit',...
            'String','-100e-5',...
            'Units','normalized','Position',[0.5 0.63 0.19 0.05]);
        
        hIntervalPropagateTo = uicontrol('Parent',rightPanel,'Style','edit',...
            'String','100e-5',...
            'Units','normalized','Position',[0.71 0.63 0.19 0.05]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uicontrol('Parent',rightPanel,'Style','text',...
            'String','Number of Steps',...
            'Units','normalized','Position',[0.1 0.48 0.4 0.1]);
        
        hNumberOfSteps = uicontrol('Parent',rightPanel,'Style','edit',...
            'String','200',...
            'Units','normalized','Position',[0.5 0.53 0.4 0.05]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uicontrol('Parent',rightPanel,...
            'String','PROPAGATE!',...
            'Units','normalized','Position',[0.1 0.38 0.8 0.1],'Callback',@performPropagation);
    
    
        function performPropagation(hObj,callbackdata)        
            %##############################################################
            % Physical part about propagation #############################
            % These parameters should be adjustable in GUI
            binning = 1;
            h = 4.1357e-15; % Plank's constant    
            c = 2.99792458e8; % Speed of light in vacuum

            % Input parameters needed for free-space propagation
            param.photonEnergy           = str2double(get(hEnergy,'String')); % [eV]
            param.sampleDetectorDistance = str2double(get(hSampleToDetector,'String')); % [m]
            param.detectorPixelSize      = str2double(get(hDetectorPixel,'String'))*binning; % [m]
            param.propagateInterval      = [str2double(get(hIntervalPropagateFrom,'String')), str2double(get(hIntervalPropagateTo,'String'))]; % Region of propagation [m]  
            param.nStep                  = str2double(get(hNumberOfSteps,'String')); % Number of sampling slices along propagation            

            waveLength = h*c/param.photonEnergy; 
            k = 2*pi/waveLength; % wave-vector in the direction of propagation

            inSize = size(IN); % size of the input data
            nX = inSize(2); % number of pixels on detector x
            nY = inSize(1); % number of pixels on detector y

            pitchX = waveLength*param.sampleDetectorDistance/(nX*param.detectorPixelSize); % Sampling step in real space in horizontal direction
            pitchY = waveLength*param.sampleDetectorDistance/(nY*param.detectorPixelSize); % Sampling step in real space in vertical direction
            pitchZ = (param.propagateInterval(2)-param.propagateInterval(1))/param.nStep;

            pitchKx = 2*pi/(nX*pitchX); % Sampling step in recipocal space in horizontal direction
            pitchKy = 2*pi/(nY*pitchY); % Sampling step in recipocal space in vertical direction
            [Kx,Ky] = meshgrid(pitchKx*(-nX/2:nX/2-1),pitchKy*(-nY/2:nY/2-1)); % Coordinates grid in reciprocal space

            xVector = -pitchX*floor(nX/2):pitchX:pitchX*ceil(nX/2);
            yVector = -pitchY*floor(nY/2):pitchY:pitchY*ceil(nY/2);               
            zVector = (param.propagateInterval(1):pitchZ:param.propagateInterval(2)); % Positions of the slices along propagation [m]            

            waveField = zeros(nY,nX,numel(zVector)); % Array of 2D slices through the beam

            wB = waitbar(0,'Blood is leaking...');  
            
            for ii = 1:numel(zVector)
                
                % Calculation of the near-field Fresnel field
        %         waveField(:,:,ii) = exp(1j*k*zVector(ii))*(nX*ifftshift(ifft2(ifftshift(exp(-1j*zVector(ii)*(Kx.^2+Ky.^2)/(2*k)).*((1/nX)*fftshift(fft2(fftshift(IN))))))));
                waveField(:,:,ii) = (nX*ifftshift(ifft2(ifftshift(exp(-1j*(-zVector(ii))*(Kx.^2+Ky.^2)/(2*k)).*((1/nX)*fftshift(fft2(fftshift(IN))))))));                    
        %         if fix(ii/numel(zVector)*100) ~= fix((ii-1)/numel(zVector)*100)    
        %             fprintf('Propagating: %ld%%\n',fix(ii/numel(zVector)*100));    
        %         end;
                waitbar(ii / numel(zVector))
            end;
            close(wB);

            sliceHorizontal = squeeze(waveField(round(nY/2),:,:));
            sliceVertical = squeeze(waveField(:,round(nX/2),:));

            scaleVector = 1e6; % Conversion of the spatial coordinates from [m] to [microns]
            xVector = xVector.*scaleVector;
            yVector = yVector.*scaleVector;
            zVector = zVector.*scaleVector.*1e-3;  % Conversion of the spatial coordinates from [m] to [mm]

            % Search for the maximum of the beam sharpness
            integratedAmplitude = squeeze(sum(sum(abs(waveField).^4)));
            pinholePosition = find(integratedAmplitude == max(integratedAmplitude));
            clear integratedAmplitude;
            % End of physics ##################################################  
            %##################################################################

            % Clear figure before drawing a new propagation
            clf;

            % Main window creation
            uicontrol('Parent',mainF,'Style','text',...
                'String','Dzhigaev Dmitry 2016',...
                'Units','normalized','Position',[0.01 0.05 0.2 0.03]);

            leftPanel = uipanel('Title','Plot control','FontSize',...
                12,'Units','Normalized','Position',[.01 0.1 .2 .85]);

            rightPanel = uipanel('Title','Propagation control',...
                'FontSize',12,'Units','Normalized',...
                'Position',[.78 .1 .2 .85]);

            % Button for openning the file 
            uicontrol('Parent',leftPanel,'Units','normalized',...
                'Position',[0.05 0.87 0.9 0.1],'String','Open file',...
                'Callback',@openFile);

            % Define the initial switches in the GUI that should be on              
            plotType  = 'Complex'; % Initial type of the plot
            sliceType = 'Transverse'; % Initial slice type
            crossOn = 'false';
            val = pinholePosition(1);

            valZoom = 1;

            imAx = axes;
        
            hIm = imagesc(xVector,yVector,complex2rgb(waveField(:,:,val)));...
                axis image; zoom(valZoom); xlabel('Horizontal coordinate, [microns]');...
                ylabel('Vertical coordinate, [microns]');
                hTitle = title(sprintf('Distance: %2.d [mm]',zVector(val)));

            % Definition of the axes for longitudinal propagation
            sliceAx1 = axes('Units','normalized','Position',[0.255 0.08 0.5 0.4],'Visible','off');      
            sliceAx2 = axes('Units','normalized','Position',[0.255 0.55 0.5 0.4],'Visible','off');

            %##################################################################
            % Right panel with a propagation parameters #######################

            pinholeFinderText = uicontrol('Parent',rightPanel,'Style','text',...
                'String',sprintf('Auto distance from the pinhole to the sample: %.2d [mm]',...
                zVector(pinholePosition)),'Units','normalized','Position',[0.1 0.88 0.8 0.1]);                

            pinholeGoBtn = uicontrol('Parent',rightPanel,...
                'String',sprintf('Go to the pinhole'), ...
                'Units','normalized','Position',[0.1 0.2 0.8 0.05],...
                'Callback', @goToPinhole);                

            distanceText = uicontrol('Parent',rightPanel,'Style','text',...
                'String','Distance to the sample position',...
                'Units','normalized','Position',[0.1 0.08 0.8 0.1]);  
            
            uicontrol('Parent',rightPanel,'Style','text',...
            'String','Photon energy, [eV]',...
            'Units','normalized','Position',[0.1 0.78 0.4 0.1]);
        
            hEnergy = uicontrol('Parent',rightPanel,'Style','edit',...
                'String',num2str(param.photonEnergy),...
                'Units','normalized','Position',[0.5 0.83 0.4 0.05]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            uicontrol('Parent',rightPanel,'Style','text',...
                'String','Sample2Detector, [m]',...
                'Units','normalized','Position',[0.1 0.68 0.4 0.1]);

            hSampleToDetector = uicontrol('Parent',rightPanel,'Style','edit',...
                'String',num2str(param.sampleDetectorDistance),...
                'Units','normalized','Position',[0.5 0.73 0.4 0.05]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            uicontrol('Parent',rightPanel,'Style','text',...
                'String','Detector pixel size, [m]',...
                'Units','normalized','Position',[0.1 0.58 0.4 0.1]);

            hDetectorPixel = uicontrol('Parent',rightPanel,'Style','edit',...
                'String',num2str(param.detectorPixelSize),...
                'Units','normalized','Position',[0.5 0.63 0.4 0.05]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            uicontrol('Parent',rightPanel,'Style','text',...
                'String','Propagation Interval, [m]',...
                'Units','normalized','Position',[0.1 0.48 0.4 0.1]);

            hIntervalPropagateFrom = uicontrol('Parent',rightPanel,'Style','edit',...
                'String',num2str(param.propagateInterval(1)),...
                'Units','normalized','Position',[0.5 0.53 0.19 0.05]);

            hIntervalPropagateTo = uicontrol('Parent',rightPanel,'Style','edit',...
                'String',num2str(param.propagateInterval(2)),...
                'Units','normalized','Position',[0.71 0.53 0.19 0.05]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            uicontrol('Parent',rightPanel,'Style','text',...
                'String','Number of Steps',...
                'Units','normalized','Position',[0.1 0.38 0.4 0.1]);

            hNumberOfSteps = uicontrol('Parent',rightPanel,'Style','edit',...
                'String',num2str(param.nStep),...
                'Units','normalized','Position',[0.5 0.43 0.4 0.05]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            uicontrol('Parent',rightPanel,...
                'String','PROPAGATE!',...
                'Units','normalized','Position',[0.1 0.28 0.8 0.1],'Callback',@performPropagation);
            
            % Slider
            slideH = uicontrol('Parent',rightPanel,'Style',...
                'slider','Min',1,'Max',ii,...
                'Value',val,'Units','Normalized',...
                'Position', [0.1 0.05 0.8 0.05],...
                'Callback', @slideBeam);       

            function goToPinhole(source,callbackdata)        
                if get(slideH,'Value') ~= val    
                    val = pinholePosition;
                    set(slideH,'Value',val);
                    slideBeam(slideH,0);
                end
            end
            % End of the right panel ##########################################
            %##################################################################

            %##################################################################
            % Left panel for plotting options #################################      

            % Selection of the plot type for the wavefield 
            bgType = uibuttongroup('Visible','off','Parent',leftPanel,...
                'Position',[0.05 0.77 0.9 0.1],'HandleVisibility',...
                'off','Title','Type');

                uicontrol(bgType,'Style','radiobutton',...
                    'String','Complex',...
                    'Units','Normalized',...
                    'Position',[0.01 0 0.3 1]);

                uicontrol(bgType,'Style','radiobutton',...
                    'String','Intensity',...
                    'Units','Normalized',...
                    'Position',[0.35 0 0.3 1]); 

                uicontrol(bgType,'Style','radiobutton',...
                    'String','Phase',...
                    'Units','Normalized',...
                    'Position',[0.7 0 0.3 1]);

                set(bgType,'SelectionChangeFcn',@typeSelection);              
                set(bgType,'Visible','on');

            % Selection of the slice through the wavefield 
            bgSlice = uibuttongroup('Visible','off','Parent',leftPanel,...
                'Position',[0.05 0.67 0.9 0.1],...                  
                'HandleVisibility','off','Title','Slice selection');       

                uicontrol(bgSlice,'Style','radiobutton',...
                    'String','Transverse',...
                    'Units','Normalized',...
                    'Position',[0.01 0 0.4 1]);

                uicontrol(bgSlice,'Style','radiobutton',...
                    'String','Longitudinal',...
                    'Units','Normalized',...
                    'Position',[0.5 0 0.4 1]);

                set(bgSlice,'SelectionChangeFcn',@sliceSelection);              
                set(bgSlice,'Visible','on'); 

            % Selection of the layer onto the plot
            bgOverlay = uibuttongroup('Visible','off','Parent',leftPanel,...
                'Position',[0.05 0.57 0.9 0.1],...                  
                'HandleVisibility','off','Title','Overlay');  

                chBox1 = uicontrol(bgOverlay,'Style','checkbox',...
                    'String','Cross',...
                    'Units','Normalized',...
                    'Position',[0.01 0 0.4 1]);

                chBox2 = uicontrol(bgOverlay,'Style','checkbox',...
                    'String','Profile',...
                    'Units','Normalized',...
                    'Position',[0.5 0 0.4 1]);

                set(chBox1,'Callback',@crossOverlay);
                set(chBox2,'Callback',@profileOverlay);
                set(bgOverlay,'Visible','on'); 

            % Zoom buttons    
            bgZoom = uibuttongroup('Visible','off','Parent',leftPanel,'Position',...
                [0.05 0.47 0.9 0.1],'HandleVisibility','off','Title','Zooming');
                
                % Zoom in button
                uicontrol(bgZoom,...
                              'String','+',...
                              'Units','Normalized',...
                              'Position',[0.1 0.1 0.3 .8],...
                              'Callback',@changeZoom);
                
                % Zoom out button
                uicontrol(bgZoom,...
                              'String','-',...
                              'Units','Normalized',...
                              'Position',[0.6 0.1 0.3 .8],...
                              'Callback',@changeZoom);

                set(bgZoom,'Visible','on');

            % Color wheel 
            colorWheel = imread('phaseWheel-01_small.png');
            axes('Parent',leftPanel,'OuterPosition',[0 0 1 0.35]);
            imshow(colorWheel);

            % Button for taking a snapshot 
            uicontrol('Parent',leftPanel,'Units','normalized','Position',...
                [0.05 0.41 .9 .05],'String','Snapshot',...
                'Callback',@saveSnapshot);

            % Button for saving a 3D data 
            uicontrol('Parent',leftPanel,'Units','normalized','Position',...
                [0.05 0.36 .9 .05],'String','Save 3D data',...
                'Callback', @save3D);

            % End of the left panel ###########################################
            %##################################################################     

            % Functions to be called ##########################################
            function slideBeam(hObj,event)
                val = round(get(hObj,'Value'));
                switch sliceType 
                    case 'Transverse'
                        axes(sliceAx1); colorbar('off'); cla; set(sliceAx1,'Visible','off');
                        axes(sliceAx2); colorbar('off'); cla; set(sliceAx2,'Visible','off');
                        axes(imAx); set(imAx,'Visible','on');
                        switch plotType
                            case 'Complex'
                               hIm = imagesc(xVector,yVector,complex2rgb(waveField(:,:,val))); axis image; zoom on; zoom(valZoom); colorbar('off');xlabel('Horizontal coordinate, [microns]'); ylabel('Vertical coordinate, [microns]');
                               if strcmp(crossOn,'true')
                                   line(xVector,zeros(1,length(xVector)),'Color','white','LineStyle',':');
                                   line(zeros(1,length(yVector)),yVector,'Color','white','LineStyle',':');
                               end
                            case 'Intensity'
                               hIm = imagesc(xVector,yVector,abs(waveField(:,:,val)).^2);axis image; zoom on; zoom(valZoom); colorbar;xlabel('Horizontal coordinate, [microns]'); ylabel('Vertical coordinate, [microns]');
                               if strcmp(crossOn,'true')
                                   line(xVector,zeros(1,length(xVector)),'Color','white','LineStyle',':');
                                   line(zeros(1,length(yVector)),yVector,'Color','white','LineStyle',':');
                               end
                            case 'Phase'
                               hIm = imagesc(xVector,yVector,angle(waveField(:,:,val)));axis image; zoom on; zoom(valZoom); colorbar;xlabel('Horizontal coordinate, [microns]'); ylabel('Vertical coordinate, [microns]');
                               if strcmp(crossOn,'true')
                                   line(xVector,zeros(1,length(xVector)),'Color','white','LineStyle',':');
                                   line(zeros(1,length(yVector)),yVector,'Color','white','LineStyle',':');
                               end
                        end
                        hTitle = title('');
                        set(hTitle,'Str',sprintf('Distance: %.2d [mm]',zVector(val)));                
                    case 'Longitudinal'
                        switch plotType
                            case 'Complex'
                                axes(imAx); colorbar('off'); cla; set(imAx,'Visible','off');             
                                set(sliceAx1,'Visible','on'); axes(sliceAx1); 
                                imagesc(zVector,xVector,complex2rgb(sliceHorizontal)); zoom yon; zoom(valZoom); xlabel('Propagation coordinate, [mm]'); ylabel('Horizontal coordinate, [microns]');
                                line(repmat(zVector(val),1,length(yVector)),yVector,'Color','white');
                                set(sliceAx2,'Visible','on'); axes(sliceAx2); set(gca,'XDir','reverse');
                                imagesc(zVector,yVector,complex2rgb(sliceVertical)); zoom yon; zoom(valZoom);xlabel('Propagation coordinate, [mm]'); ylabel('Vertical coordinate, [microns]');
                                line(repmat(zVector(val),1,length(yVector)),yVector,'Color','white');
                                if strcmp(crossOn,'true')
                                   axes(sliceAx1);
                                   line(zVector,zeros(1,length(zVector)),'Color','white','LineStyle',':');
                                   line(zeros(1,length(yVector)),yVector,'Color','white','LineStyle',':');
                                   axes(sliceAx2);
                                   line(zVector,zeros(1,length(zVector)),'Color','white','LineStyle',':');
                                   line(zeros(1,length(yVector)),yVector,'Color','white','LineStyle',':');
                                end

                            case 'Intensity'
                                axes(imAx); colorbar('off'); cla; set(imAx,'Visible','off');             
                                set(sliceAx1,'Visible','on'); axes(sliceAx1);
                                imagesc(zVector,xVector,abs(sliceHorizontal).^2);colorbar; zoom yon; zoom(valZoom); xlabel('Propagation coordinate, [mm]'); ylabel('Horizontal coordinate, [microns]');
                                line(repmat(zVector(val),1,length(yVector)),yVector,'Color','white');
                                set(sliceAx2,'Visible','on'); axes(sliceAx2);
                                imagesc(zVector,yVector,abs(sliceVertical).^2);colorbar; zoom yon; zoom(valZoom); xlabel('Propagation coordinate, [mm]'); ylabel('Vertical coordinate, [microns]');
                                line(repmat(zVector(val),1,length(yVector)),yVector,'Color','white');
                                if strcmp(crossOn,'true')
                                   axes(sliceAx1);
                                   line(zVector,zeros(1,length(zVector)),'Color','white','LineStyle',':');
                                   line(zeros(1,length(yVector)),yVector,'Color','white','LineStyle',':');
                                   axes(sliceAx2);
                                   line(zVector,zeros(1,length(zVector)),'Color','white','LineStyle',':');
                                   line(zeros(1,length(yVector)),yVector,'Color','white','LineStyle',':');
                                end
                            case 'Phase'
                                axes(imAx); colorbar('off'); cla; set(imAx,'Visible','off');             
                                set(sliceAx1,'Visible','on'); axes(sliceAx1);
                                imagesc(zVector,xVector,angle(sliceHorizontal));colorbar; zoom yon; zoom(valZoom); xlabel('Propagation coordinate, [mm]'); ylabel('Horizontal coordinate, [microns]');
                                line(repmat(zVector(val),1,length(yVector)),yVector,'Color','white');
                                set(sliceAx2,'Visible','on'); axes(sliceAx2);
                                imagesc(zVector,yVector,angle(sliceVertical));colorbar; zoom yon; zoom(valZoom); xlabel('Propagation coordinate, [mm]'); ylabel('Vertical coordinate, [microns]');
                                line(repmat(zVector(val),1,length(yVector)),yVector,'Color','white');
                                if strcmp(crossOn,'true')
                                   axes(sliceAx1);
                                   line(zVector,zeros(1,length(zVector)),'Color','white','LineStyle',':');
                                   line(zeros(1,length(yVector)),yVector,'Color','white','LineStyle',':');
                                   axes(sliceAx2);
                                   line(zVector,zeros(1,length(zVector)),'Color','white','LineStyle',':');
                                   line(zeros(1,length(yVector)),yVector,'Color','white','LineStyle',':');
                                end                        
                        end                
                end
                set(distanceText,'Str',sprintf('Distance: %.2d [mm]',zVector(val)))
            end                

            % Function change the current image type immediately 
            function typeSelection(source,callbackdata)
                plotType = get(callbackdata.NewValue,'String');                
                slideBeam(slideH,0);
            end                              

            % Function change the current image profile immediately 
            function sliceSelection(source,callbackdata)
                sliceType = get(callbackdata.NewValue,'String');                
                slideBeam(slideH,0);
            end

            % Function adds the layers to the current image immediately 
            function crossOverlay(source,callbackdata)
                if (get(source,'Value') == get(source,'Max'))
                     crossOn = 'true';
                elseif (get(source,'Value') == get(source,'Min'))
                     crossOn = 'false';
                end             
                slideBeam(slideH,0);
            end

            % NOT IMPLEMENTED YET
            function profileOverlay(source,callbackdata)
                if (get(source,'Value') == get(source,'Max'))
                     display('Selected');
                end
            end

            % Changing the constant zooming of the views
            function changeZoom(source,callbackdata)
                zoomInc = get(source,'String');
                switch zoomInc
                    case '+'
                        valZoom = valZoom+0.5;
                    case '-'
                        if (valZoom-0.5)>=1
                            valZoom = valZoom-0.5;
                        end
                end
                slideBeam(slideH,0);
            end                

            % Save snapshot of the current view
            function saveSnapshot(hObject, eventdata)                            
                filePath = fullfile(param.pathName,'snapshots',sprintf('%s',timeStamp));
                mkdir(filePath);
                display('Snapshots of the beam is being save...');
                set(mainF,'PaperPositionMode','auto');
                print(mainF,strcat(filePath,sprintf('/fullFrame_%s_%s_%.2d mm.png',plotType,sliceType,zVector(val))),'-dpng','-r300');
                print(mainF,strcat(filePath,sprintf('/fullFrame_%s_%s_%.2d mm.eps',plotType,sliceType,zVector(val))),'-depsc');            
                display(sprintf('Snapshots of the beam is saved in %s.',filePath));            
            end

            % Save the full 3D beam matrix
            function save3D(hObject, eventdata)                
                filePath = fullfile(param.pathName,'arrays3D',sprintf('%s',timeStamp));
                mkdir(filePath);
                display('3D array of the beam is being saved...');            
                save(strcat(filePath,'/beam3D.mat'),'waveField','-v7.3');
                save(strcat(filePath,'/realSpacePitch.mat'),'pitchX','pitchY','pitchZ')
                display(sprintf('3D array of the beam is saved in %s. Array size: [%d %d %d]',filePath,nY,nX,numel(zVector)));
                fid = fopen(strcat(filePath,'/beam3d.raw'),'w');
                   absWaveField = abs(waveField);
                   fwrite(fid,'absWaveField','double');
                fclose(fid);
            end
        end
    end
end