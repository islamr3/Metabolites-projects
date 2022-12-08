

%Percentage=Number of case at each year/Total Population
%For Normal Albumuniria
load ("NA.mat")
B=A(:,2:4)
%No. of patients with 40% eGFR decline
figure(1)
subplot(2,1,1)
bar(B)
xlabel('year')
ylabel('40% eGFR decline')
legend('T1','T2','T3')
title("Only Positive events")



Cum = cumsum(B);
%Cumulative
subplot(2,1,2)
bar(Cum)
xlabel('year')
ylabel('Cum with 40% eGFR decline')
legend('T1','T2','T3')
title("Only Positive events with cumulative ways")


%Load data
clear all;
load("NAData.mat")

B=A(1:7,2:4)
%No. of patients with 40% eGFR decline
figure(2)
subplot(2,1,1)
bar(B)
xlabel('year')
ylabel('40% eGFR decline')
legend('T1','T2','T3')
title("No of positve patients divided by total number of patients at each tertile")

Cum = cumsum(B);
%Cumulative
subplot(2,1,2)
bar(Cum)
xlabel('year')
ylabel('Cum patients with 40% eGFR decline (in %)')
legend('T1','T2','T3')


