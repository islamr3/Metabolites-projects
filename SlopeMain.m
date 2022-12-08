%rabi=mdl.Coefficients.Estimate(2,1)

load('egfrData.csv');
%Crreate visit time
T=[1:9]';
%remove adenine from egfrData
Adnine=egfrData(:,10);
%delete adnine colunm from egfrData
egfrData(:,10)=[];
SlopData=[];
for i=1:size(egfrData)
    %Extract each patient egrfData
    egFPt1=egfrData(i,:)';    %
    [mdl] = fitlm(egFPt1,T);
    %Access slope from the model
    slop=mdl.Coefficients.Estimate(2,1);
    %Store slope in SlopData
    SlopData(i,:)=slop;
end

%estimate correlation between Adenine SlopData
NormalizeAdnine=log(Adnine);
%NormalizeAdnine=Adnine/max(abs(Adnine));
plot(SlopData);hold on;plot(NormalizeAdnine,'r');hold off;

R = corrcoef(SlopData,NormalizeAdnine)
x=[NormalizeAdnine,SlopData];
h = ttest(x)

