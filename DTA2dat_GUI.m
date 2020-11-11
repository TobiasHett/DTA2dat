function DTA2dat_master_guinew(~,~)
clc, clear, close all
% Check proper setting of the EasySpin Path
if (exist('bmagn')|| exist('planck') || exist('eprload'))==0
    msgbox('Error! Please check that your EasySpin path is set properly!','DTA2dat GUI','Error')
    disp('Error! Please check that your EasySpin path is set properly!')
    return
end

fig=figure('Name','DTA2dat GUI','resize','off','NumberTitle','off');

a = uicontrol(fig,'Style','pushbutton','String','Select File','ToolTip','Select Single .DTA file','Units','Normalized'); % Push button to load Single .DTA file.
a.Position = [0,0.95,0.1,0.05]; % Position of pushbutton (xpos, ypos, xwidth, ywidth)
b = uicontrol(fig,'Style','pushbutton','String','Select Folder','ToolTip','Select all .DTA-files in a folder','Units','Normalized'); % Push button to select folder.
b.Position = [0.1,0.95,0.2,0.05]; % Position of pushbutton (xpos, ypos, xwidth, ywidth)
About = uicontrol(fig,'Style','pushbutton','String','About','ToolTip','About DTA2dat GUI...','Units','Normalized');
About.Position = [0.9,0.95,0.1,0.05];

Infoline2D = uicontrol(fig,'Style','text','String','Options for conversion of 2-dimensional .DTA-files','Units','Normalized','BackgroundColor','w'); % Infoline regarding 2D-files.
Infoline2D.Position = [0.0,0.85,0.5,0.05];
c = uicontrol(fig,'Style','checkbox','String','Save separate .dat-files','ToolTip','For 2D .DTA: Save each slice to separate .dat-file','Units','Normalized'); % Checkbox to determine whether the each slice should be saved separately.
c.Position = [0,0.80,0.3,0.05]; % Position of checkbox (xpos, ypos, xwidth, ywidth)
d = uicontrol(fig,'Style','checkbox','String','Save concatenated .dat-file','ToolTip','For 2D .DTA: Save one .dat-file containing all slices','Units','Normalized'); % Checkbox to determine whether all slices should be saved in one .dat-file.
d.Position = [0,0.75,0.3,0.05]; % Position of checkbox (xpos, ypos, xwidth, ywidth)

Infolinefrqshift = uicontrol(fig,'Style','text','String','Options for frequency shifting and g-scale conversion','ToolTip','Frequency Shifting and g-scale-computation currently only works for 1D Field Sweeps!','Units','Normalized','BackgroundColor','w'); % Infoline regarding 2D-files.
Infolinefrqshift.Position = [0.0,0.65,0.5,0.05];
e = uicontrol(fig,'Style','checkbox','String','Save raw data as .dat-file','ToolTip','Convert raw .DTA into .dat-file','Units','Normalized'); % Checkbox to determine whether all slices should be saved in one .dat-file.
e.Position = [0,0.6,0.3,0.05]; % Position of checkbox (xpos, ypos, xwidth, ywidth)
f = uicontrol(fig,'Style','checkbox','String','Save as g-scale','ToolTip','Save spectra with g-scale axis','Units','Normalized'); % Checkbox to determine whether all slices should be saved in one .dat-file.
f.Position = [0,0.55,0.3,0.05]; % Position of checkbox (xpos, ypos, xwidth, ywidth)
g = uicontrol(fig,'Style','checkbox','String','Save at target frequency (in GHz):','ToolTip','Shift magnetic field axis','Units','Normalized'); % Checkbox to determine whether all slices should be saved in one .dat-file.
g.Position = [0,0.5,0.43,0.05]; % Position of checkbox (xpos, ypos, xwidth, ywidth)
targfrq = uicontrol(fig,'Style','edit','String','9.7','ToolTip','Target Frequency in GHz','Units','Normalized'); % Checkbox to determine whether all slices should be saved in one .dat-file.
targfrq.Position = [0.335,0.5,0.1,0.05]; % Position of checkbox (xpos, ypos, xwidth, ywidth)

Infolinemirror = uicontrol(fig,'Style','text','String','Options for mirroring PDS traces','ToolTip','Mirror symmetric PDS traces (e.g. DQC) at the maximum point','Units','Normalized','BackgroundColor','w');
Infolinemirror.Position = [0.0,0.40,0.5,0.05];
h = uicontrol(fig,'Style','checkbox','String','Save raw data as .dat-file','ToolTip','Convert raw .DTA into .dat-file','Units','Normalized');
h.Position = [0,0.35,0.43,0.05]; % Position of checkbox (xpos, ypos, xwidth, ywidth)
%Do not use i here since it is used in loops in the functions!
j = uicontrol(fig,'Style','checkbox','String','Mirror PDS trace at maximum','ToolTip','Mirror symmetric PDS traces at the maximum','Units','Normalized');
j.Position = [0,0.30,0.43,0.05]; % Position of checkbox (xpos, ypos, xwidth, ywidth)

Statusbar = uicontrol(fig,'Style','text','String','Ready','ToolTip','Statusbar','Units','Normalized','BackgroundColor','w'); % Statusbar
Statusbar.Position = [0.0,0.15,1,0.05];

a.Callback = @getfile;
b.Callback = @getdir;
About.Callback = @aboutwindow;
%% Get files and directories
function getfile(src,event)
[fileName, pathName] = uigetfile('*.DTA','Select .DTA file...');              %Select .DTA file.

    if e.Value==1 || f.Value==1 || g.Value==1
       singlefrqshift(src,event,fileName,pathName)
    elseif h.Value==1 && j.Value==1
        singleconv(src,event,fileName,pathName)
        mirrortrace(src,event,fileName,pathName)
    elseif h.Value==1
       singleconv(src,event,fileName,pathName)
    elseif j.Value==1
        mirrortrace(src,event,fileName,pathName)
    else
    singleconv(src,event,fileName,pathName)
    end
end
function getdir(src,event)
pathName = uigetdir();              %Select folder.

    if e.Value==1 || f.Value==1 || g.Value==1
       multifrqshift(src,event,pathName)
    else
        multiconv(src,event,pathName)
    end
end

%% Single file conversion
function singleconv(~,~,fileName,pathName) %Only one single file selected (either 1D or 2D .DTA).
%Script for 1D .DTA file.    
    filePath = strcat(pathName,fileName);
    [x, y, Pars] = eprload(filePath);                        %Read .DTA file.
        
if isfield(Pars, 'YPTS')==0                             
        fileName = fileName(1:(end-3));                         %Create .dat file
        fileName = strcat(fileName,'dat');
        filePath = strcat(pathName,fileName);
        File = fopen(filePath,'w');
    
        for i = 1:size(x,1)
            fprintf(File,'%f %f %f\n',x(i),real(y(i)),imag(y(i)));  %Write ASCII data to file.
        end
        fclose(File); 
        else
        
%Script for 2D .DTA file.
        x_spalte1=x{1,1};                                       %Pre-processing: Split and convert two-dimensional spectrometer data.
        x_spalte2=x{1,2};
        x_spalte1=x_spalte1';
        x_spalte2=x_spalte2';
        y_dimension=size(y);
        
%Create separate .dat files for 2D type. 
  if c.Value==1 % If checkbox is activated, generate separate files.
        for i=1:y_dimension(2)                                 
            if exist('fileName2')==1                            %Check if filename variable exists and delete it.
                clearvars ('fileName2')
            end
            fileName2 = fileName(1:(end-4));                    %Generate new filename.
            fileName3 = strcat(fileName2,'_Slice_',num2str(i));
            fileName3 = strcat(fileName3,'.dat');
            
            filePath = strcat(pathName,fileName3);
            File = fopen(filePath,'w');                         %Write ASCII data to file.
            y_wert=y(:,i);
            fprintf(File,'%f %f\n', [x_spalte1,y_wert]');
            fclose(File);
        end
  end
  
%Export all slices as one .dat file
    if d.Value==1 % If checkbox is activated, generate concatenated file.
            fileName2 = fileName(1:(end-4));                    %Generate new filename.
            filePath = strcat(pathName,fileName2,'_total.dat'); 

            File = fopen(filePath,'w');
            start=' %f ';
            string=' %f ';
            for i=1:y_dimension(2)
                string=strcat(string, start);
            end
            string=strcat(string, ' \n');
            fprintf(File,string, [x_spalte1,y]');
            fclose(File);
            
%Pre-process y-axis data for 2D files and write them into separate _Params.dat file.     
        x_spalte2flip=x_spalte2';
        x_spalte2flip=fliplr(x_spalte2flip);
        x_spalte2flip=x_spalte2flip';
        
        filePathParam = strcat(pathName,fileName2,'_Params.dat');
        File = fopen(filePathParam,'w');
            fprintf(File,'%f \n', x_spalte2flip);
        fclose(File);
    end
end
    updatestatusbar_1file
end
%% Conversion of multiple files
function multiconv(src,event,pathName)
% % Generate list of .DTA files in the given folder.
  files = dir( fullfile(pathName,'*.DTA') );
 
  pathName=[pathName, '\'];
for numoffiles=1:numel(files)
 clearvars -except numoffiles files pathName c d e f g h j Statusbar targfrq src event;
 
    daten=struct2cell(files(numoffiles));
     fileName=char(daten(1,1));
% Read .DTA file
        filePath = strcat(pathName,fileName);
        [x, y, Pars] = eprload(filePath); 
        
% Determine type of .DTA file (1D or 2D)
    if isfield(Pars, 'YPTS')==0     % Then it is a 1D .DTA file
       
% Create .DAT file for 1D type
        fileName = fileName(1:(end-3));
        fileName = strcat(fileName,'dat');
        filePath = strcat(pathName,fileName);
        File = fopen(filePath,'w');
        for i = 1:size(x,1)
            fprintf(File,'%f %f %f\n',x(i),real(y(i)),imag(y(i)));
        end
        fclose(File);
        
    else                            % Then it is a 2D .DTA file
% 
% Pre-processing: Split and convert two-dimensional spectrometer data  
        x_spalte1=x{1,1};
        x_spalte2=x{1,2};
        x_spalte1=x_spalte1';
        x_spalte2=x_spalte2';
        y_dimension=size(y);
% 
% Create .DAT file for 2D type        
 if c.Value==1        
        for i=1:y_dimension(2)
            if exist('fileName2')==1
                clearvars ('fileName2')
            end
            fileName2 = fileName(1:(end-4));
            fileName3 = strcat(fileName2,'_Slice_',num2str(i));
            fileName3 = strcat(fileName3,'.dat');
            
            filePath = strcat(pathName,fileName3);
            File = fopen(filePath,'w');
            y_wert=y(:,i);
            fprintf(File,'%f %f\n', [x_spalte1,y_wert]');
            fclose(File);
        end
end
             fileName2 = fileName(1:(end-4));
% 
%Export all slices as one .dat file
    if d.Value==1
            filePath = strcat(pathName,fileName2,'_total.dat');
            File = fopen(filePath,'w');
            
            start=' %f ';
            string=' %f ';
            for i=1:y_dimension(2)
            string=strcat(string, start);
            end
            string=strcat(string, ' \n');
            fprintf(File,string, [x_spalte1,y]');
            fclose(File);
%       
% Pre-process y-axis data for 2D files and write them into separate _Params.dat file.
        x_spalte2flip=x_spalte2';
        x_spalte2flip=fliplr(x_spalte2flip);
        x_spalte2flip=x_spalte2flip';
        
        filePathParam = strcat(pathName,fileName2,'_Params.dat');
        File = fopen(filePathParam,'w');
            fprintf(File,'%f \n', x_spalte2flip);
            fclose(File);
    end
    end  
    updatestatusbar_multifile(src,event,numoffiles, files)
end
end

%% Frequency shift of one file
function singlefrqshift(src,event,fileName,pathName)  
filePath = strcat(pathName,fileName);
[x, y, Pars] = eprload(filePath);                        %Read .DTA file.

    if strcmp(Pars.XNAM,'Field')==1 

    if isfield(Pars, 'YPTS')~=0
        Statusbar.String='Error! Only 1D field sweeps are currently supported!';
    else
       if e.Value==1                                    %Simple conversion of raw data.
           singleconv(src,event,fileName,pathName)
       end
       
       if f.Value==1                                    %Compute g-scale.
           Field=x/10000;                               %Conversion Gauss to Tesla
           Zaehler=planck*Pars.MWFQ;
           Nenner=bmagn*Field;
           gscale=Zaehler./Nenner;                      %Transform to g-scale
           dlmwrite(strcat(filePath(1:end-4),'_gscale.dat'),{gscale;real(y)}','delimiter','\t','precision',10)
       end 
       
        if g.Value==1                                   %Compute new magnetic field axis.
           Field=x/10000;                               %Conversion Gauss to Tesla.
           Zaehler=planck*Pars.MWFQ;
           Nenner=bmagn*Field;
           gscale=Zaehler./Nenner;                      %Transform to g-scale
           newField=planck*(str2num(targfrq.String)*1e+9)./(bmagn.*gscale);
           newFieldvector=newField*10000;
           dlmwrite(strcat(filePath(1:end-4),'_targfreq_',targfrq.String,'GHz.dat'),{newFieldvector;real(y)}','delimiter','\t','precision',10)       
        end 
    end
    else
    Statusbar.String='Error! This is not a Field Swept Spectrum!';
    disp('This is not a Field Swept Spectrum!')
    end
end
%% Frequency shift of multiple files
function multifrqshift(src,event,pathName)
      files = dir( fullfile(pathName,'*.DTA') );
      pathName=[pathName, '\'];
    
   for numoffiles=1:numel(files)
        clearvars -except numoffiles files pathName a b c d e f g h j Statusbar src event targfrq;
        daten=struct2cell(files(numoffiles));
        fileName=char(daten(1,1));
 
% Read .DTA file
        filePath = strcat(pathName,fileName);
        [~, ~, Pars] = eprload(filePath);
% Determine type of .DTA file (1D or 2D); process only 1D .DTA files.

    if isfield(Pars, 'YPTS')==0     % Then it is a 1D .DTA file
        singlefrqshift(src,event,fileName,pathName)        
     end
    end
end
%% Mirror PDS trace
function mirrortrace(~,~,fileName,pathName)
    filePath=strcat(pathName,fileName);
    [x, y, Pars] = eprload(filePath)  ;

if strcmp(Pars.XNAM,'Time')==1 
    
       signal_realpart=real(y);
       signal_impart=imag(y);
       [realmaxval,realmaxrow]=max(signal_realpart); % realmaxval = value; realmaxrow = rownumber
       [immaxval,immaxrow]=max(signal_impart); % immaxval = value; immaxrow = rownumber
    
    if realmaxval>immaxval %Real part has a higher signal amplitude than the imaginary part
% Extract right and left part of the trace next to the maximum.
 time=x(realmaxrow:end);
 rightpart_real=signal_realpart(realmaxrow:end); % Right part of trace
 rightpart_im=signal_impart(realmaxrow:end);
 
 leftpart_real=signal_realpart(1:realmaxrow); % Left part of trace
 leftpart_im=signal_impart(1:realmaxrow);
 
% Flip left part of time trace prior to computing of mean values.
 leftpart_real_flip=flipud(leftpart_real);
 leftpart_im_flip=flipud(leftpart_im);
 
% Bring left and right time traces to the same vector lengths.
 relevantsize=min(size(leftpart_real_flip),size(rightpart_real)); %Size of vector
 
rightpart_real=rightpart_real(1:relevantsize(1));
rightpart_im=rightpart_im(1:relevantsize(1));
leftpart_real_flip=leftpart_real_flip(1:relevantsize(1));
leftpart_im_flip=leftpart_im_flip(1:relevantsize(1));
 
% Compute the mean values.
meanreal=(rightpart_real+leftpart_real_flip)/2;
meanim=(rightpart_im+leftpart_im_flip)/2;
 
% Generate export array.
 export(:,1)=time(1:length(meanreal));
 export(1:end,2)=meanreal; %Averaged value of left and right part
 export(1:end,3)=meanim;
 
 dlmwrite(strcat(filePath(1:end-4),'_mirrored.dat'),export,'delimiter','\t','precision','%.20f');
 
 else %Imaginary part has a higher signal amplitude than the real part
time=x(immaxrow:end);
rightpart_real=signal_realpart(immaxrow:end); % Right part of trace
rightpart_im=signal_impart(immaxrow:end);

leftpart_real=signal_realpart(1:immaxrow); % Left part of trace
leftpart_im=signal_impart(1:immaxrow);

% Flip left part of time trace prior to computing of mean values.
leftpart_real_flip=flipud(leftpart_real);
leftpart_im_flip=flipud(leftpart_im);
 
% Bring left and right time traces to the same vector lengths.
relevantsize=min(size(leftpart_real_flip),size(rightpart_real)); %Size of vector

rightpart_real=rightpart_real(1:relevantsize(1));
rightpart_im=rightpart_im(1:relevantsize(1));
leftpart_real_flip=leftpart_real_flip(1:relevantsize(1));
leftpart_im_flip=leftpart_im_flip(1:relevantsize(1));

% Compute the mean values.
meanreal=(rightpart_real+leftpart_real_flip)/2;
meanim=(rightpart_im+leftpart_im_flip)/2;
    
% Generate export array.
export(:,1)=time(1:length(meanreal));
export(1:end,2)=meanreal; %Averaged value of left and right part
export(1:end,3)=meanim;
 
dlmwrite(strcat(pathName,fileName(1:end-4),'_mirrored.dat'),export,'delimiter','\t','precision','%.20f');
  
    end
    else
    Statusbar.String='Error! This is not a PDS trace!';
    disp('This is not a PDS trace!')
end

end

%% Update statusbar if one file is processed
function updatestatusbar_1file
    Statusbar.String='Done!';
    Statusbar.BackgroundColor='b';
    Statusbar.ForegroundColor='w';
    pause(2)
    Statusbar.String='Ready';
    Statusbar.BackgroundColor='w';
    Statusbar.ForegroundColor='k';
end

%% Update statusbar if multiple files are processed
function updatestatusbar_multifile(~,~,numoffiles, files)
    
    Statusbar.String=strcat(num2str(numoffiles/length(files)*100),' %');
    Statusbar.Position = [0.0,0.15,numoffiles/length(files),0.05];
    Statusbar.BackgroundColor='b';
    Statusbar.ForegroundColor='w';
    drawnow
    if numoffiles==length(files)
        pause(1)
       Statusbar.String='Done! All files converted!';
       Statusbar.BackgroundColor='w';
       Statusbar.ForegroundColor='k';
       pause(0.75)
       Statusbar.String='Ready';
    end
end
%% Infobox about-window
function aboutwindow(~,~)
msgbox({'DTA2dat GUI';'Version 1.0';'September 2020'},'About','Help')
end
end
