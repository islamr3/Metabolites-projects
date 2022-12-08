%rabi=mdl.Coefficients.Estimate(2,1)
%https://www.mathworks.com/help/stats/linearmodel.plotresiduals.html

%Tukey’s outlier test
%https://www.mathworks.com/help/econ/archtest.html
%https://www.mathworks.com/help/stats/residuals.html
%%%%%%%mdl = fitlm(X,y)
%scatter(NormalizeAdnine,mdl1.Residuals.Raw)
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Male group analysis
data=readtable('MaleEGFR_Normal.csv');  %P-value Significant
%data=readtable('MaleEGFR_Micro.csv');
%data=readtable('MaleEGFR_Macro.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Male group analysis
%data=readtable('FemaleEGFR_Normal.csv');
%data=readtable('FemaleEGFR_Micro.csv');
%data=readtable('FemaleEGFR_Macro.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% All patient for Normal
%data=readtable('AllEGFR_Normal.csv');
%data=readtable('AllEGFR_Micro.csv');
%data=readtable('AllEGFR_Macro.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Tertile analysis%%%%%%%%%

%%%%%%%%%%%%%%%%%% Male Normal Gruup %%%%%%%%%%%%%%%%%%%
%%%Lower tertile
%data=readtable('LowTertile_Male_Normal.csv');

%%%%Uper tertile
%data=readtable('UPTertile_Male_Normal.csv'); %P-value is significant


%%%%%%%%%%%%%%%%%% Male Micro Gruup %%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%Uper tertile for all patients for Normal group
%data=readtable('UpTertile_Normal.csv');

%%%%%%%%%%%%Lower tertile for all patients for Normal group
%data=readtable('LowTertile_Normal.csv');


size(data);
j=1;
SlopData=[];
Y=[];
r11=1;
ra=[];
for rows=1:size(data,1)
    d=data(rows,1:13); %Extract each row with eGFR
    %Romiving data point from each row   
    K=1;
    t=1;
    for i=1:size(d,2)
        g=d.(i); %Extract each data point frm d
        if ~isnumeric(g)
            c=str2double(g); 
            if ~isnan(c)
                Y(t)=i;
                t=t+1;
            end
        else
            if ~isnan(g)
                Y(t)=i;
                t=t+1;
            else
            end
        end
        if (isnumeric(g)) %store each data point of g to p (removing NAN)
            p(j)=g;             
            j=j+1;
        else
           g1=str2double(g);
           if (isnan(g1))
               NN=K+1;
           else
               p(j)=g1;
               j=j+1;
           end
        end
    end
    p(find(isnan(p)))=[]; %Removing Nan from P
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %T=[1:size(p,2)]; %Crreate visit time %T is defined based on the visit (need correct)
     T=Y;
    %Apply model to each patient 
    if size(p,2)>=3
        [mdl] = fitlm(p,T);
        %[mdl] = fitlm(p,T,'RobustOpts','on');
        %mdl
        slop=mdl.Coefficients.Estimate(2,1);
        %Store slop to SlopData
        SlopData(r11,:)=slop; 
        r11=r11+1;
    else
         %Track the index with less than 3
        ra=[ra,rows];
    end
    %egFPt1=data(i,:)';
    p=[];
    Y=[];
    j=1;
end
%Extract adnine from data and preprocessing

adnine=data.(size(data,2));
adnine(ra)=[];
[idx]=find(isnan(adnine));
SlopData(idx)=[];
adnine(idx)=[];
%adnine(find(isnan(adnine)))=[]; %Removing Nan from P
%estimate estimate log adnine
NormalizeAdnine=log(adnine);
histogram(NormalizeAdnine)
%Normalize both adnine and slopData
%NormalizeAdnine=NormalizeAdnine/max(abs(NormalizeAdnine));
NormalizeAdnine = zscore(NormalizeAdnine)';

% % TF = isoutlier(NormalizeAdnine);
% % [id]=find(TF==1);
% % NormalizeAdnine(id)=[];

%NormalizeAdnine=NormalizeAdnine;
%SlopData=abs(SlopData/max(abs(SlopData)));
SlopData = zscore(SlopData)';
dataN=[NormalizeAdnine;SlopData];
[rankings,feaWeights] = SRCFS(dataN);
id=rankings(1:20);
SlopData(id)=[];
NormalizeAdnine(id)=[];
%SlopData=SlopData;
histogram(NormalizeAdnine);hold on;h=histogram(SlopData,'FaceColor', 'k');hold off;
legend('Adenine','eGFR slope')
xlabel('Patients Distribution', 'FontSize', 16);
ylabel('Frequency', 'FontSize', 16);
set(gca,'FontSize',12)


%NormalizeAdnine=Adnine/max(abs(Adnine));
%plot(SlopData);hold on;plot(NormalizeAdnine,'r');hold off;
%Distrivution estimation
% DisData=[NormalizeAdnine,SlopData];
% figure(1);
% for i=1:2
%     distributionfitRabi1(DisData(:,i));hold on;
% end
% hold off;
%estimate correlation between SlopData and NormalizeAdnine
%r = corr(SlopData,NormalizeAdnine,'rows','complete')
%r = corr(SlopData,NormalizeAdnine)
%r2=r*r % R-squared for correlation 

%x=[NormalizeAdnine,SlopData];
%[h,p,ci,stats] = ttest(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Create model with SlopData and Adnine
[mdl1] = fitlm(NormalizeAdnine,SlopData);
%[mdl1] = fitlm(NormalizeAdnine,SlopData,'RobustOpts','on')
mdl1
figure(2);
h=plot(mdl1,'Color', 'k', 'LineWidth', .5)
%figure(2);plot(mdl2)

% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
% The confidence bounds have 2 handles but only one of 
% the handles contains the legend string.  The first
% line below finds that object and then searches for 
% other objects in the plot that have the same linestyle
% and color. 
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);


dataHandle.Color = 'k'; 
fitHandle.Color = [0 0 0]; %orange 
set(cbHandles, 'Color', 'b', 'LineWidth', 2)


title('Data and model','FontSize', 16);
xlabel('log Adenine', 'FontSize', 16);
ylabel('eGFR Slope', 'FontSize', 16);
set(gca,'FontSize',12)
%text(x, y, 'Hey, look at this', 'FontSize', 24);






