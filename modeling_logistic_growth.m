%model cell proliferation can explain notochord changes 
%until increase in cell size in vertebra is important later

%load manualcell.mat on desktop, data from 12/21/23, using stardist in
%qupath to count cells
%column + data:
% 57 is number of Nootochord cell in ivd level, 
% 58 is number of notochord cells in ivd level which is positive for edu
% 61 is number of notochord cells in vb level
%62 is number of notochord cells in ivd level which is positive for edu
%11 is stage
%t is thickness of cell
%load manualcell
th = 15;
manualcell.total_Noto_cell = manualcell.Noto_IVD_Total + manualcell.Noto_VB_total;


e14_stage_logic = (manualcell.Stage == 14);
embryo_logic = manualcell.Stage == 12;
e12_noto_total = mean(manualcell.total_Noto_cell(manualcell.Stage == 12 & manualcell.Noto_IVD_Total<100),'omitnan')

e14_noto_total = mean(manualcell.total_Noto_cell(e14_stage_logic & manualcell.Noto_IVD_Total>10),'omitnan')


e17_noto_total = mean(manualcell.total_Noto_cell(manualcell.Stage == 17),'omitnan')
e14_model = e14_noto_total + 3*(mean(manualcell.Noto_IVD_Edu(e14_stage_logic),'omitnan')+mean(manualcell.Noto_VB_Edu(e14_stage_logic),'omitnan'))
e12_model = e12_noto_total + 5*(mean(manualcell.Noto_IVD_Edu(manualcell.Stage == 12 & manualcell.Noto_IVD_Total<100),'omitnan')+mean(manualcell.Noto_VB_Edu(manualcell.Stage == 12 & manualcell.Noto_IVD_Total<100),'omitnan'))

%% estimate number of cells in 3d volume for e17 total, e14 total and e12 total,
% So the plan (see dfLab notebook nov2023.docx) is multiply each timepoint
% region (vert vs ivd) cell count by a geometric factor to expand a 2d
% slice to the 3d volume based on measured dimensions (radius, height)
% assume each 2d slice has thickness th
% for each cell count, total and edu: 
% For E14 ivd, assume sphere, multiply by 4r/3t
e14_ivd_total = manualcell.Noto_IVD_Total(e14_stage_logic).*((manualcell.Notochord_IVD_diameter(e14_stage_logic)/2)*4/(3*th));
e14_ivd_edu = manualcell.Noto_IVD_Edu(e14_stage_logic).*((manualcell.Notochord_IVD_diameter(e14_stage_logic)/2)*4/(3*th));

% For e14 vb, assume cylinder, multiply by r*pi/2t
e14_vb_total = manualcell.Noto_VB_total(e14_stage_logic).*((manualcell.Notochord_InferiorVB_diameter(e14_stage_logic)/2)*pi()/(2*th));
e14_vb_edu = manualcell.Noto_VB_Edu(e14_stage_logic).*((manualcell.Notochord_InferiorVB_diameter(e14_stage_logic)/2)*pi()/(2*th));

% For e12 ivd, assume cylinder, multiply by r*pi/2t
e12_ivd_total = manualcell.Noto_IVD_Total(manualcell.Stage == 12).*((manualcell.Notochord_IVD_diameter(manualcell.Stage == 12)/2)*pi()/(2*th));
e12_ivd_edu = manualcell.Noto_IVD_Edu(manualcell.Stage == 12).*((manualcell.Notochord_IVD_diameter(manualcell.Stage == 12)/2)*pi()/(2*th));
% For e12 vb, assume cylinder, multiple by r*pi/2t
e12_vb_total = manualcell.Noto_VB_total(manualcell.Stage == 12).*((manualcell.Notochord_InferiorVB_diameter(manualcell.Stage == 12)/2)*pi()/(2*th));
e12_vb_edu = manualcell.Noto_VB_Edu(manualcell.Stage == 12).*((manualcell.Notochord_InferiorVB_diameter(manualcell.Stage == 12)/2)*pi()/(2*th));
%for e17, assume sphere
e17_noto_total_volume = manualcell.total_Noto_cell(manualcell.Stage == 17).*((manualcell.Notochord_IVD_diameter(manualcell.Stage == 17)/2)*4/(3*th));


%% if ignoring 3d estimation
e14_ivd_total = manualcell.Noto_IVD_Total(e14_stage_logic);
e14_ivd_edu = manualcell.Noto_IVD_Edu(e14_stage_logic);

% For e14 vb
e14_vb_total = manualcell.Noto_VB_total(e14_stage_logic);
e14_vb_edu = manualcell.Noto_VB_Edu(e14_stage_logic);

% For e12 ivd
e12_ivd_total = manualcell.Noto_IVD_Total(manualcell.Stage == 12);
e12_ivd_edu = manualcell.Noto_IVD_Edu(manualcell.Stage == 12);
% For e12 vb
e12_vb_total = manualcell.Noto_VB_total(manualcell.Stage == 12);
e12_vb_edu = manualcell.Noto_VB_Edu(manualcell.Stage == 12);
%for e17
e17_noto_total_volume = manualcell.total_Noto_cell(manualcell.Stage == 17);


%-----

%estimates assuming linear growth (adding the number of prolfierating cells
%every day of growth
e14_vol_model = e14_ivd_total + 3*e14_ivd_edu + e14_vb_total +3*e14_vb_edu;
e12_vol_model = e12_ivd_total + 5*e12_ivd_edu + e12_vb_total +5*e12_vb_edu;


%divide number of cells by density at each timepoint to estimate notochord
%area. 

%% run the logistic model for every E12 and E14 data point
%for E14 estimate, run fro 3 days, for the E12 estimate, run for 5 days
%Define the parameters for the logistic model:
%define growth rate as proliferation rate
%need e12 growth rate and the number of noto_vb and noto_ivd cells
%next, use volume estimates from 2D, increase r linearly. 

K = 3*median(e17_noto_total_volume,'omitnan'); % carrying capacity
for ii = 1:length(e12_vb_edu)
    %r = 0.1; % growth rate
    r =(e12_ivd_edu(ii) + e12_vb_edu(ii))/(e12_ivd_total(ii)+e12_vb_total(ii));
    b(ii) = r;%Define the initial condition and time steps:
    y0 = [e12_ivd_total(ii)+e12_vb_total(ii)]; % initial population
    tspan = [0 5]; % time span
    dt = 0.1; % time step
    
    %Use ode45 to solve the ODE using Euler%s method:
    [t, y2] = ode45(@(t, y2) logistic_growth_notochord(t, y2, r, K), tspan, y0);
    y2(t==5);
    e17_e12model(ii) = y2(t==5);
%     %once I find the growth, replace NaN's with 0 to add ivd and vb estimates
%     %together
%     figure
%     plot(t, y2);
%     xlabel('Time');
%     ylabel('Population');
%     title('Logistic Cell Population Growth from E12');
    
end
mean(b,'omitnan') %0.1201
mean(e17_e12model,'omitnan')
clear r

K = 3*median(e17_noto_total_volume,'omitnan'); % carrying capacity
for ii = 1:length(e14_vb_edu)
    %r = 0.7; % growth rate
    r = (e14_ivd_edu(ii) + e14_vb_edu(ii))/(e14_ivd_total(ii)+e14_vb_total(ii));
    f(ii) = r;%Define the initial condition and time steps:
    y0 = [e14_ivd_total(ii)+e14_vb_total(ii)]; % initial population
    tspan = [0 3]; % time span
    dt = 0.1; % time step
    
    %Use ode45 to solve the ODE using Euler%s method:
    [t, y1] = ode45(@(t, y1) logistic_growth_notochord(t, y1, r, K), tspan, y0);
    y1(t == 3);
    e17_e14model(ii) = y1(t == 3);
%     %once I find the growth, replace NaN's with 0 to add ivd and vb estimates
%     %together
%     figure
%     plot(t, y1);
%     xlabel('Time');
%     ylabel('Population');
%     title('Logistic Cell Population Growth from E14');
    
end
mean(f,'omitnan') %0.3143
mean(e17_e14model,'omitnan')

%plot the logistic cell count estimates from e14 and e12 medians with the 3d e17 volume estimation
% figure
% y = [e17_e12model e17_e14model e17_noto_total_volume'];
% x = [1*ones(size(e17_e12model)) 2*ones(size(e17_e14model)) 3*ones(size(e17_noto_total_volume))'];
% s = swarmchart(x,y,5,'filled');
% hold on
% s.SizeData = 20;
% colormap(jet)
% s.XJitterWidth = 0.5;
% title('Logistic Growth Estimates and actual 3D estimate of Notochord Cell Number')
% ylabel('Cells')

%% calculate area given cell count estimates and measured 

%plot the logistic cell count estimates from e14 and e12 medians with the 3d e17 volume estimation
figure
manualcell.density_Noto(e14_stage_logic);
y = [e17_e12model./manualcell.density_Noto(manualcell.Stage == 12)' e17_e14model./manualcell.density_Noto(e14_stage_logic)' e17_noto_total_volume'./manualcell.density_Noto(manualcell.Stage == 17)'];

med_e12_model = median((e17_e12model'./manualcell.density_Noto(manualcell.Stage == 12)),'omitnan')
med_e14_model = median(e17_e14model./manualcell.density_Noto(e14_stage_logic)','omitnan')
med_e17_actual = median(e17_noto_total_volume'./manualcell.density_Noto(manualcell.Stage== 17)','omitnan')
std_e12_model = std((e17_e12model'./manualcell.density_Noto(manualcell.Stage == 12)),'omitnan')
std_e14_model = std(e17_e14model./manualcell.density_Noto(e14_stage_logic)','omitnan')
std_e17_actual = std(e17_noto_total_volume'./manualcell.density_Noto(manualcell.Stage== 17)','omitnan')
x = [1*ones(size(e17_e12model)) 2*ones(size(e17_e14model)) 3*ones(size(e17_noto_total_volume))'];
s = swarmchart(x,y,5,'filled');
hold on
s.SizeData = 20;
colormap(jet)
s.XJitterWidth = 0.5;
title('Logistic Growth Estimates and Actual Notochord Area')
ylabel('Area (um^2)')
%ss = gca;
xticks([1 2 3])
xticklabels({'Estimate from E12 Proliferation' , 'Estimate from E14 Proliferation' , 'E17 Actual Area'})
%ss.XTickLabel =  { 'Estimate from E12 Proliferation'  'E14 Proliferation Estimate'  'E17 Actual Area' }


d = errorbar([med_e12_model med_e14_model med_e17_actual],[std_e12_model std_e14_model std_e17_actual],'x')
d.LineWidth = 2;
d.Color = [0 0 0];

%stop and make plot area smaller be grabbing outside of box
print('-dtiff',['notochord_area_model_20240328','.tif']);
m00010('notochord_area_model_20240328')




%% plotting





figure
bar([e12_noto_total e14_noto_total e17_noto_total])
hold on
title('2d estimation of Notochord Cell number/level')
ylabel('number of cells')
d = gca;
d.FontSize = 15;
d.XTickLabel = {'e12','e14', 'e17'};

figure
daynames = ["E12.5" "E14.5" "E17.5"];
c = (manualcell.Level_number);
x = categorical(manualcell.Stage,[12,14,17],daynames);
y = manualcell.total_Noto_cell;
%c = (manualcell.Level_number);
s = swarmchart(x,y,5,c,'filled');
hold on
s.SizeData = 20;
colormap(jet)
s.XJitterWidth = 0.5;
title('2D Estimation of Notochord Cell Number')
ylabel('Cells')

figure
y = [[e12_ivd_total+e12_vb_total]' [e14_ivd_total+e14_vb_total]' e17_noto_total_volume'];
x = [ones(size(e12_vol_model))' 2*ones(size(e14_vol_model))' 3*ones(size(e17_noto_total_volume))'];
s = swarmchart(x,y,5,'filled');
hold on
s.SizeData = 20;
colormap(jet)
s.XJitterWidth = 0.5;
title('3D estimate of Notochord Cell Number')
ylabel('Cells')

figure
y = [e12_vol_model' e14_vol_model' e17_noto_total_volume'];
x = [ones(size(e12_vol_model))' 2*ones(size(e14_vol_model))' 3*ones(size(e17_noto_total_volume))'];
s = swarmchart(x,y,5,'filled');
hold on
s.SizeData = 20;
colormap(jet)
s.XJitterWidth = 0.5;
title('Linear Growth at e17 with 3D Estimation')
ylabel('Cells')






figure
d = bar([median(e12_vol_model,'omitnan') median(e14_vol_model,'omitnan') median(e17_noto_total_volume,'omitnan')])
hold on
title('3d, linear growth')
ylabel('number of cells')

d = gca;
d.FontSize = 15;
d.XTickLabel = {'e12 linear growth','e14 linear growth', 'actual e17'};



%% run logistic model on median
%Define the parameters for the logistic model:

%define growth rate as proliferation rate
r = 0.7; % growth rate
K = 3*median(e17_noto_total_volume,'omitnan'); % carrying capacity

%Define the initial condition and time steps:
y0 = [(median(e14_ivd_total,'omitnan')+median(e14_vb_total,'omitnan'))]; % initial population
tspan = [0 5]; % time span
dt = 0.1; % time step

%Use ode45 to solve the ODE using Euler%s method:
[t, y] = ode45(@(t, y) logistic_growth_notochord(t, y, r, K), tspan, y0);
%once I find the growth, replace NaN's with 0 to add ivd and vb estimates
%together
e17_e14model = y(t>3 & t<3.1)
figure
plot(t, y);
xlabel('Time');
ylabel('Population');
title('Logistic Cell Population Growth from E14, r = 0.7');

%Define the parameters for the logistic model:
%define growth rate as proliferation rate
r = 0.1; % growth rate
K = 1.5*median(e17_noto_total_volume,'omitnan'); % carrying capacity

%Define the initial condition and time steps:
y0 = [(median(e12_ivd_total,'omitnan')+median(e12_vb_total,'omitnan'))]; % initial population
tspan = [0 5]; % time span
dt = 0.1; % time step

%Use ode45 to solve the ODE using Euler%s method:
[t, y2] = ode45(@(t, y2) logistic_growth_notochord(t, y2, r, K), tspan, y0);
e17_e12model = y2(t==5)
%once I find the growth, replace NaN's with 0 to add ivd and vb estimates
%together
figure
plot(t, y2);

xlabel('Time');
ylabel('Population');
title('Logistic Cell Population Growth from E12, r = 0.1');

%plot the logistic estimates from e14 and e12 medians with the 3d e17 volume estimation
figure
y = [e17_e12model e17_e14model e17_noto_total_volume'];
x = [1 2 3*ones(size(e17_noto_total_volume))'];
s = swarmchart(x,y,5,'filled');
hold on
s.SizeData = 20;
colormap(jet)
s.XJitterWidth = 0.5;
title('Logistic Growth Estimates and actual estimate of Notochord Cell Number')
ylabel('Cells')
%% saving the data
%collect ydata from model projections in notochord_area_model_20240821
%which has the same data as notochord_area_model_20240328, just bigger
%font, and saved as .emf
%tt is product of gco when the data are highlihgted in the .fig
%i.e. and object which has the data in tt.YData
%the E12.5 is the first 44 points, the E14.5 is the second 44 points, and
%E17.5 is the last 44 points
E12_model_fromplot = tt.YData(1:44);
E14_model_fromplot = tt.YData(45:88);
E17_model_fromplot = tt.YData(89:132);
med_diffE12e17 = median([(E17_model_fromplot-E12_model_fromplot)./E12_model_fromplot],'omitnan')
stdev_diffE12e17 = std([(E17_model_fromplot-E12_model_fromplot)./E12_model_fromplot],'omitnan')

med_diffE14e17 = median([(E17_model_fromplot-E14_model_fromplot)./E14_model_fromplot],'omitnan')
stdev_diffE14e17 = std([(E17_model_fromplot-E14_model_fromplot)./E14_model_fromplot],'omitnan')
