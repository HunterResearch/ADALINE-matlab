
%%% THIS CODE CREATES A QUANTILE PLOT OF THE ADALINE RUNS FROM DATA FILES
%%% THE DATA FILES ARE GENERATED BY servo_parfor.m


clear; clc;
close all;


% set up parameters
%myOracleName = 'BusScheduling';
myOracleName = 'DiscreteQuadratic';

strreps = '10'; %string with number of replications, for loading data files

% This function saves a specifically-sized pdf file of a figure with an 
% array of rowsubs x colsubs sub-plots.

%***********************************************************************
%   Define output 
%***********************************************************************
% How wide should the output figure be (inches)?
fig_wd = 6.5; %4.7;%<-actual size 
% How tall should the output figure be (inches)?
fig_ht = fig_wd/1.8; %3;
% What is the primary font size for each plot?
fs=10;
% What is the secondary font size for each plot?
fs2=8;
% How many sub-plots per row?
rowsubs = 2;
% How many sub-plots per column?
colsubs = 2;
% How much of a gap should there be between subplots (inches)?
space = 0.05;
%***********************************************************************

% Create figure
figure('Name','Plot','Units','inches','Position',[0 0 fig_wd fig_ht]);
set(0,'DefaultTextInterpreter','latex')

% Figure out normalized width/height of sub-plots
sp_wd1 = (fig_wd-space*(rowsubs+1))/colsubs;
sp_wd = sp_wd1/fig_wd;
sp_ht1 = (fig_ht-space*(colsubs+1))/rowsubs;
sp_ht = sp_ht1/fig_ht;
hspace = space/fig_wd;
vspace = space/fig_ht;


% Create Sub-plots
for j=1:2%rowsubs  
    for i=1:2%colsubs
        
        %% LOAD DATA FILES FOR THE CURRENT PLOT
        if j==1 %FIRST ROW
            if i==1 %FIRST PLOT
                dim = '25';
            else
                dim = '50';
            end
        else %SECOND ROW
            if i==1 %FIRST PLOT
                dim = '100';
            else
                dim = '200';
            end
        end
        
        %LOAD ADALINE
        myfilename=strcat('Data',myOracleName,'ADALINE',dim,'reps',strreps,'.mat');
        load(myfilename,'quantilecurve','qtpts');
        ADALINEquantilecurve = quantilecurve;
        
        %LOAD RSPLINE
        myfilename=strcat('Data',myOracleName,'RSPLINE',dim,'reps',strreps,'.mat');
        load(myfilename,'quantilecurve','qtpts');
        RSPLINEquantilecurve = quantilecurve;
            
        
        % Create sub-plot at custom position in the figure (based on spacing determined above)
        pos = [hspace*i+sp_wd*(i-1) 1-(vspace*j+sp_ht*j) sp_wd sp_ht];
        subplot('Position',pos)
        
        ax = gca; % get subplot axes
        % set position of sub-plot AXES in the figure (based on spacing determined above)
        ax.OuterPosition = [hspace*i+sp_wd*(i-1) 1-(vspace*j+sp_ht*j) sp_wd sp_ht];
        %set secondary font for elements of sub-plot
        set(gca,'fontsize',fs2,'fontname','Times');
    
        myq = [];
        hold on;
        myylim = 0;
        for ell = [3,4,5,7]%1:size(qtpts,2)
            p2=plot(ADALINEquantilecurve(:,1),ADALINEquantilecurve(:,ell),'-k');%,'Color',[.4 .4 .4]);
            p1=plot(RSPLINEquantilecurve(:,1),RSPLINEquantilecurve(:,ell),':k');
            myylimnew = max(max(RSPLINEquantilecurve(:,ell)),max(ADALINEquantilecurve(:,ell)));
            myylim = max([myylimnew myylim])/2;
            if ell == 1
                myq = [myq 10];
            end
            if ell == 2
                myq = [myq 15];
            end
            if ell == 3 
                myq = [myq 25];
            end
            if ell == 4
                myq = [myq 50];
            end
            if ell == 5
                myq = [myq 75];
            end
            if ell == 6
                myq = [myq 85];
            end
            if ell == 7
                myq = [myq 90];
            end
        end
        grid on
        %hold off;
        xlim([0 , max(ADALINEquantilecurve(:,1))])
        if j == 1 && i == 1
            ylim([0 , 30]);
            yticks([0  10  20  30 ])
        elseif j==1 && i == 2
            ylim([0 , 40]);
            yticks([0 10 20 30 40])
        elseif j==2 && i == 1
            ylim([0 , 40]);
            yticks([0 10 20 30 40])
        else
            ylim([0 , 50]);
            yticks([0 10 20 30 40 50])
        end   
        %ylim([0 , myylim]);
        xlabel('Total Oracle Calls $t$','Interpreter','latex','fontsize',fs);
        ylabel('Opt. Gap \%','Interpreter','latex','fontsize',fs);
        myqstr = '';
        for k = 1:length(myq)
            nextq = num2str(myq(k));
            if k<length(myq)
                myqstr = [myqstr nextq ', '];
            else
                myqstr = [myqstr nextq];
            end
        end
        titl_string = ['\fbox{\textbf{Ill-conditioned Quadratic, $d=$ ',dim,'}}'];%: ',myqstr,'}}'];
        titl_stringlgd = [myqstr ' %-iles'];
        %title(titl_string,'Interpreter','latex','fontsize',fs);
        if j==1 && i==1
            lgd=legend([p1 p2],{'R-SPLINE','ADALINE'},'Position',[0.336 0.8260 0.1549 0.0919]);
            title(lgd,titl_stringlgd)
            title(titl_string,'Interpreter','latex','fontsize',fs2);
        elseif j==1 && i==2
            %lgd=legend([p1 p2],{'R-SPLINE','ADALINE'},'Position',[0.8310 0.8413 0.1549 0.0919]) ;
            %title(lgd,titl_stringlgd)
            title(titl_string,'Interpreter','latex','fontsize',fs2);
        elseif j==2 && i==1
            %lgd=legend([p1 p2],{'R-SPLINE','ADALINE'},'Position',[0.3285 0.3498 0.1549 0.0919]); 
            %title(lgd,titl_string)
            title(titl_string,'Interpreter','latex','fontsize',fs2);
        elseif j==2 && i==2
            %lgd=legend([p1 p2],{'R-SPLINE','ADALINE'},'Position',[0.8310 0.3498 0.1549 0.0919]); 
            %title(lgd,titl_string)
            title(titl_string,'Interpreter','latex','fontsize',fs2);
        end
        %lgd.Position %to find out where the legend is; use with 'Location'
        %'northeast'
        
        % get rid of whitespace around plot
        % from https://www.mathworks.com/help/matlab/creating_plots/save-figure-with-minimal-white-space.html
        ax = gca;                       % get current sub-plot axes
        outerpos = ax.OuterPosition;    % get OuterPosition of sub-plot axes
        ti = ax.TightInset;             % get TightInset of sub-plot axes
        left = outerpos(1) + ti(1);     % make adjustments based on ticks, labels, etc...
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height]; % adjust the sub-plot size
   
    end
end

% the next steps set the output size of the figure
fig = gcf;                      % get the current figure
fig.PaperPositionMode = 'auto'; 
fig_pos = fig.PaperPosition;    % get the current PaperPosition
fig.PaperSize = [fig_pos(3) fig_pos(4)]; % Set PaperSize to width/height of PaperPosition

myfigname=strcat('Plot',myOracleName,'reps',strreps);
print(fig,myfigname,'-dpdf') % save to pdf file