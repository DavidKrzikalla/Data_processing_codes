%% This script serves for evaluation of 3-point bending experiment according to ASTM D790
% Author: David Krzikalla, david.krzikalla@vsb.cz

% The script outline:
% 0) Import the data from Excel sheets
%	- import experimental data, specimens dimensions and supports span values
% 1) Process the data
%	- The data are manipulated to begin from zero force and then shifted accordingly (since at the begining the force is zero until the contact with specimen is established but the deflection is measured)
% 2) Get flexural stress, flexural strength, flexural strain according to standatd 12.2, 12.4, 12.8
% 3) Get flexural modulus according to standatd 12.9.1 using 'm' as slope of force-deflection linear part taken from 25 % up to 75 % maximal force region
% 4) Plot/write the results

%% Initiation
clc
clear all
close all

%% 0) Import the data
addpath('C:\Users\David Krzikalla\OneDrive - VSB-TUO\Testovani_materialu_Markforged');

% Import force/deflection data from experiments for five specimens 
experimental=xlsread('E_017_5mm_min.xlsx','Test Curve Data','A3:J5764');

% Import dimensions of specimens tested. Width, thickness
dimensions=xlsread('Ukoly a testovani.xlsx','Ohyb','D40:E44');

% Import supports span 
span=xlsread('Ukoly a testovani.xlsx','Ohyb','H4');


% Import data for comparison plot> experimental/numerical
addpath('C:\Users\David Krzikalla\OneDrive - VSB-TUO\Testovani_materialu_Markforged\FEA\Flexural_test\shell_solid1');
addpath('C:\Users\David Krzikalla\OneDrive - VSB-TUO\Testovani_materialu_Markforged\FEA\Flexural_test\shell2');
addpath('C:\Users\David Krzikalla\OneDrive - VSB-TUO\Testovani_materialu_Markforged\FEA\Flexural_test\solid2');
addpath('C:\Users\David Krzikalla\OneDrive - VSB-TUO\Testovani_materialu_Markforged\FEA\Flexural_test\solid3');

shell1=xlsread('shell_solid1_results.xlsx','List1','B2:C21'); % read shell 1 results, defl/force
solid1=xlsread('shell_solid1_results.xlsx','List1','F2:G21'); % read solid 1 results
shell2=xlsread('shell2_results.xlsx','List1','N2:O21'); % read shell2 results
solid2=xlsread('solid2_results.xlsx','List1','B2:C21'); % read solid 2 results
solid3=xlsread('solid3_results.xlsx','List1','N2:O21'); % read solid 3 results
%% 1) Process the data
% Locate the beginning of the loading (locate the zero force possitions within first 100 rows)

X=zeros(100,1);
R=zeros(size(experimental,2)/2,1);
for i=2:2:size(experimental,2)
    for j=1:100
        X(j,1)=experimental(j,i);
    end
    [R(i/2,1),~]= find(~X,1,'last'); % find last zero force, point when the contact is established 
end

% Create a new experimental results matrix with ommited zero force values
% and shifted back to zero
experimental_new=zeros(size(experimental,1),size(experimental,2));
for i=1:size(R,1)
    for j=1:size(experimental,1)-(R(i,1)-1) 
    experimental_new(j,i*2-1)=(experimental(R(i,1)+(j-1),i*2-1))-(experimental(R(i,1),i*2-1)); %defl
    experimental_new(j,i*2)=experimental(R(i,1)+(j-1),i*2); %force
    if experimental_new(j,i*2-1)>25
        break
    end
    end
end

%% 2) Get flexural stress, flexural strength, flexural strain
SigmaF_ULT=zeros(size(R,1),1); % initiate a vector for flexural strength, 
for i=1:size(R,1)
    SigmaF_ULT(i,1)=(3*max(experimental_new(:,i*2))*span)/(2*dimensions(i,1)*dimensions(i,2)^2);
end

SigF_EpsF=zeros(size(experimental_new,1),2*size(R,1)); % initiate a matrix for flexural stress and strain
for i=1:size(R,1)
    for j=1:size(experimental_new,1)
        SigF_EpsF(j,i*2-1)=(3*experimental_new(j,i*2)*span)/(2*dimensions(i,1)*dimensions(i,2)^2); % flexural stress
        SigF_EpsF(j,i*2)=(6*dimensions(i,2)*experimental_new(j,i*2-1))/span^2; % flexural strain
    end
end

%% 3) Get flexural modulus
P=zeros(size(R,1),3);
L=zeros(size(R,1),2);
m=zeros(size(R,1),1);
EF=zeros(size(R,1),1);
TOL_L1=0.1; % tolerance for finding the index in data matrix of the  60 % Pmax value
TOL_L2=0.2; % tolerance for finding the index in data matrix of the  10 % Pmax value
for i=1:size(R,1)
    P(i,1)=max(experimental_new(:,i*2)); % get max force
    P(i,2)=P(i,1)*0.6; % get 60 % of max force
    P(i,3)=P(i,1)*0.1; % get 10 % of max force
    
    [L(i,1),~]=find(abs(experimental_new(:,i*2)-P(i,2))<TOL_L1,1); % find index of the 60 % Pmax
    [L(i,2),~]=find(abs(experimental_new(:,i*2)-P(i,3))<TOL_L2,1); % find index of the 10 % Pmax
    % If there is an issue 'Assignment has more non-singleton rhs dimensions
    % than non-singleton subscripts' then just increase the tolerance
    % TOL_L1 and/or TOL_L2 by 0.1 until working
    
    m(i,1)=(experimental_new(L(i,1),i*2)-experimental_new(L(i,2),i*2))/(experimental_new(L(i,1),i*2-1)-experimental_new(L(i,2),i*2-1)); % get slope of the selected linear region of force-deflection curve

    EF(i,1)=(span^3*m(i,1))/(4*dimensions(i,1)*dimensions(i,2)^3); % get flexural modulus
end
Q=find(experimental_new(:,2)==P(1,1),1); % find index of the Pmax

T=zeros(1,size(R,1));
for i=2:2:size(experimental,2)
    [T(1,i/2),~]= find(experimental_new(2:size(experimental_new,1),i),1,'last'); % find first zero force at the end of the matrix (additional zeros to be exclueded from plotting) 
end

%% 4) Plot/write the results
SigmaF_ULT % write flexural strength for each specimen from 1
EF % write flexural modulus for each specimen from 1

figure % plot flexural stress vs strain for all specimens
plot(SigF_EpsF(1:T(1,1),2),SigF_EpsF(1:T(1,1),1)); 
hold on
for i=2:size(R,1)
	plot(SigF_EpsF(1:T(1,i),i*2),SigF_EpsF(1:T(1,i),i*2-1)); 
end	
xlabel('Flexural Strain [-]','FontSize',12)
ylabel('Flexural Stress [MPa]','FontSize',12)
legend({'Specimen 1','Specimen 2','Specimen 3','Specimen 4','Specimen 5'},'Location','southeast')
grid on
print('SS','-dpng');
print('SS','-dsvg');

figure % plot force vs mid-span deflection for all specimens
plot(experimental_new(1:T(1,1),1),experimental_new(1:T(1,1),2)); 
hold on
for i=2:size(R,1)
    plot(experimental_new(1:T(1,i),i*2-1),experimental_new(1:T(1,i),i*2)); 
end	
xlabel('Mid-Span Deflection [mm]')
ylabel('Force [N]')
legend({'Specimen 1','Specimen 2','Specimen 3','Specimen 4','Specimen 5'},'Location','southeast')
grid on
print('FD','-dpng');
print('FD','-dsvg');

figure % plot force vs mid-span deflection for experimental sp. 1 and numerical
plot(experimental_new(1:Q,1),experimental_new(1:Q,2),'-k'); % plot exprimental curve of the first specimen up to Pmax
hold on
plot(shell1(:,1),shell1(:,2),'--^');
plot(solid1(:,1),solid1(:,2),'--o');
plot(shell2(:,1),shell2(:,2),'--*');
plot(solid2(:,1),solid2(:,2),'--s');
plot(solid3(:,1),solid3(:,2),'--h');
xlabel('Mid-Span Deflection [mm]')
ylabel('Force [N]')
legend({'Experiment','Shell 1','Solid 1','Shell 2','Solid 2','Solid 3'},'Location','southeast')
grid on
print('exp_num','-dpng');
print('exp_num','-dsvg');