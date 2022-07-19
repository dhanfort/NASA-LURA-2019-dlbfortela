% SOM for WWTP datasets derived from the SOM_DEMO2 of the SOM Toolbox.

% SOM_DEMO2 Basic usage of the SOM Toolbox.
% Contributed to SOM Toolbox 2.0, February 11th, 2000 by Juha Vesanto
% http://www.cis.hut.fi/projects/somtoolbox/

% Version 1.0beta juuso 071197 
% Version 2.0beta juuso 070200

clf reset;
figure(gcf)
echo on


clc

%    STEP 1: CONSTRUCT DATA
%    ======================

%    The SOM Toolbox has a special struct, called data struct, which
%    is used to group information regarding the data set in one
%    place.

%    Here, a data struct is created using function SOM_DATA_STRUCT.
%    First argument is the data matrix itself, then is the name 
%    given to the data set, and the names of the components
%    (variables) in the data matrix.


pause % Strike any key to read and construct data from directory file...


clc

% Read the data file from directory
sDdata = som_read_data('SOM_base1_data1_var4_oxygen_TEXT_tabDelimited.txt');


% Convert data into data structure workable in SOM Toolbox
%sDiris = som_data_struct(sDdata);
%sDdata = sDdata;

%%     Here are the histograms and scatter plots of the four variables.

% k=1;
% for i=1:10, 
%   for j=1:10, 
%     if i==j, 
%       subplot(10,10,k); 
%       hist(sDiris.data(:,i)); title(sDiris.comp_names{i})
%     elseif i<j, 
%       subplot(10,10,k); 
%       plot(sDiris.data(:,i),sDiris.data(:,j),'k.')
%       xlabel(sDiris.comp_names{i})
%       ylabel(sDiris.comp_names{j})
%     end
%     k=k+1;
%   end
% end


%     Actually, most SOM Toolbox functions
%     can also handle plain data matrices, but then one is without the
%     convenience offered by component names, labels and
%     denormalization operations.


pause % Strike any key to normalize the data...



%% SOM calcs

clc

%    STEP 2: DATA NORMALIZATION
%    ==========================

%    Since SOM algorithm is based on Euclidian distances, the scale of
%    the variables is very important in determining what the map will
%    be like. If the range of values of some variable is much bigger
%    than of the other variables, that variable will probably dominate
%    the map organization completely. 

%    For this reason, the components of the data set are usually
%    normalized, for example so that each component has unit
%    variance. This can be done with function SOM_NORMALIZE:


sDdata = som_normalize(sDdata,'var');

%    The function has also other normalization methods.

%    However, interpreting the values may be harder when they have
%    been normalized. Therefore, the normalization operations can be
%    reversed with function SOM_DENORMALIZE:


x = sDdata.data(1,:);

orig_x = som_denormalize(x,sDdata);


pause % Strike any key to to train the map...



%%

clc

%    STEP 3: MAP TRAINING
%    ====================

%    The function SOM_MAKE is used to train the SOM. By default, it
%    first determines the map size, then initializes the map using
%    linear initialization, and finally uses batch algorithm to train
%    the map.  Function SOM_DEMO1 has a more detailed description of
%    the training process.

%S=5srt(N)=346
%sMap = som_make(sDiris, 'msize', [21 10]); %qe=__, te=__
%sMapc = som_make(sDdata, 'msize', [20 10], 'shape', 'cyl'); %cylindrical topol
sMapr = som_make(sDdata, 'msize', [20 10], 'shape', 'sheet'); %cylindrical topol

%sMap = som_make(c, 'msize', [12 20], 'shape', 'toroid'); %toroid topol
%som_show_gui(sMap);



pause % Strike any key to continue...


%    The IRIS data set also has labels associated with the data
%    samples. Actually, the data set consists of 50 samples of three
%    species of Iris-flowers (a total of 150 samples) such that the
%    measurements are width and height of sepal and petal leaves. The
%    label associated with each sample is the species information:
%    'Setosa', 'Versicolor' or 'Virginica'.

%    Now, the map can be labelled with these labels. The best
%    matching unit of each sample is found from the map, and the
%    species label is given to the map unit. Function SOM_AUTOLABEL 
%    can be used to do this: 

%sMap = som_autolabel(sMap,sDiris,'vote');

pause % Strike any key to visualize the map...



%%

clc
%clf
%    STEP 4: VISUALIZING THE SELF-ORGANIZING MAP: SOM_SHOW
%    =====================================================

%    The basic visualization of the SOM is done with function SOM_SHOW.

%colormap(flipud(summer))
%colormap(summer)
%colormap(flipud(hot)) % we have been using this
%colormap(1-gray)
%som_show(sMap,'umat',[1 3 5 6 7 8 9 10 11 12 13],'comp', [1 3 5 6 7 8 9 10 11 12 13],'norm','d')
%som_show(sMapr,'umat','all','comp', [1 3 5 6 7 8 9 10 11 12 13],'norm','d')

som_show(sMapr,'umat','all')


%sMapc_d.codebook = som_denormalize(sMapc.codebook, sMapc);
sMapr_d.codebook = som_denormalize(sMapr.codebook, sMapr);
%Coc = som_unit_coords(sMapc.topol.msize,'hexa','sheet');
Cor = som_unit_coords(sMapr.topol.msize,'hexa','sheet');


%% Rectangulart SOM

wd = 300; %pixels
ht = 200; %pixels

% fig1 = figure('Name','Sheet SOM', 'Color','white');
% set(fig1, 'Position',  [50, 50, wd, ht])
% %subplot(1,5,1);
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,21),'Markersize',1,'Line','none');
% colormap(jet); colorbar;
% %cb=colorbar('eastoutside'); %set(cb,'position',[.15 .1 .1 .3]);
% %set(ax1,'ColorScale','log');
% %cb.Ruler.Scale = 'log';
% title('S_{NH4+}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% % print(fig_Snh4,'-dtiff','-r100','S_nh4_v3.tif')
% % saveas(fig_Snh4,'S_nh4_v2.tif')
% % export_fig('S_nh4_v4.tif', fig_Snh4)
% 
% fig2 = figure('Name','Sheet SOM', 'Color','white');
% set(fig2, 'Position',  [50, 50, wd, ht])
% %subplot(1,5,2)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,22),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('S_{NH3}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% 
% fig3 = figure('Name','Sheet SOM', 'Color','white');
% set(fig3, 'Position',  [50, 50, wd, ht])
% %subplot(1,5,3)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,23),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('S_{NO3}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig4 = figure('Name','Sheet SOM', 'Color','white');
% set(fig4, 'Position',  [50, 50, wd, ht])
% %subplot(1,5,4)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,24),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('S_{O2}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig5 = figure('Name','Sheet SOM', 'Color','white');
% set(fig5, 'Position', [50, 50, wd, ht])
% %subplot(1,5,5)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,25),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('S_{CO2}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% %
% %figure('Name','Sheet SOM', 'Color','white');
% 
% fig6 = figure('Name','Sheet SOM', 'Color','white');
% set(fig6, 'Position',  [50, 50, wd, ht])
% %subplot(1,5,1)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,26),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('S_{HCO3-}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig7 = figure('Name','Sheet SOM', 'Color','white');
% set(fig7, 'Position',  [50, 50, wd, ht])
% %subplot(1,5,2)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,27),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('S_{CO3^2-}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% 
% fig8 = figure('Name','Sheet SOM', 'Color','white');
% set(fig8, 'Position',  [50, 50, wd, ht])
% %subplot(1,5,3)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,28),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('S_{H+}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig9 = figure('Name','Sheet SOM', 'Color','white');
% set(fig9, 'Position',  [50, 50, wd, ht])
% %subplot(1,5,4)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,29),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('S_{OH-}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig10 = figure('Name','Sheet SOM', 'Color','white');
% set(fig10, 'Position',  [50, 50, wd, ht])
% %subplot(1,5,5)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,30),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('S_{X algae}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off


% Plot the parameter ranks
%
%figure('Name','Sheet SOM', 'Color','white');


% fig11 = figure('Name','Sheet SOM', 'Color','white');
% set(fig11, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,1)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,1),'Markersize',0.5,'Line','none');
% %colormap(jet(2)), colorbar
% %  cmap = colormap(winter(20)) ; %Create Colormap
% %  cbh = colorbar ; %Create Colorbar
% %  cbh.Ticks = linspace(1, 1, 20) ; %Create 8 ticks from zero to 1
% %  cbh.TickLabels = num2cell(1:20) ;    %Replace the labels of these 8 ticks with the numbers 1 to 8
% 
%  cmapdef = colormap(jet) ; %Define Colormap
%  cmap = cmapdef(1:20:end, :) ; %Find Values of colors corresponding to each point plotted
%  colorbar('YTickLabel', num2cell(1:20)) ;
% 
% title('u_{alg}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off

fig12 = figure('Name','Sheet SOM', 'Color','white');
set(fig12, 'Position',  [50, 50, wd, ht])
%subplot(4,5,2)
som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,2),'Markersize',0.5,'Line','none');
colormap(jet), colorbar
title('k_{resp}')%, view(-100,0)%, 
axis tight, axis equal
axis off
% 
% fig13 = figure('Name','Sheet SOM', 'Color','white');
% set(fig13, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,3)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,3),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('k_{death}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig14 = figure('Name','Sheet SOM', 'Color','white');
% set(fig14, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,4)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,4),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('K_{C}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig15 = figure('Name','Sheet SOM', 'Color','white');
% set(fig15, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,5)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,5),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('I_{CO2}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig16 = figure('Name','Sheet SOM', 'Color','white');
% set(fig16, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,6)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,6),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('K_{N}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig17 = figure('Name','Sheet SOM', 'Color','white');
% set(fig17, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,7)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,7),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('K_{O2}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig18 = figure('Name','Sheet SOM', 'Color','white');
% set(fig18, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,8)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,8),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('K_{pr}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig19 = figure('Name','Sheet SOM', 'Color','white');
% set(fig19, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,9)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,9),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('\tau')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig20 = figure('Name','Sheet SOM', 'Color','white');
% set(fig20, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,10)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,10),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('T_{opt}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig21 = figure('Name','Sheet SOM', 'Color','white');
% set(fig21, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,11)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,11),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('s')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig22 = figure('Name','Sheet SOM', 'Color','white');
% set(fig22, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,12)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,12),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('\alpha')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig23 = figure('Name','Sheet SOM', 'Color','white');
% set(fig23, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,13)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,13),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('\beta')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig24 = figure('Name','Sheet SOM', 'Color','white');
% set(fig24, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,14)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,14),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('\gamma')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig25 = figure('Name','Sheet SOM', 'Color','white');
% set(fig25, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,15)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,15),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('\delta')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig26 = figure('Name','Sheet SOM', 'Color','white');
% set(fig26, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,16)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,16),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('Ka_{O2}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig27 = figure('Name','Sheet SOM', 'Color','white');
% set(fig27, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,17)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,17),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('Ka_{CO2}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig28 = figure('Name','Sheet SOM', 'Color','white');
% set(fig28, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,18)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,18),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('Ka_{NH3}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig29 = figure('Name','Sheet SOM', 'Color','white');
% set(fig29, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,19)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,19),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('T_{act}')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% fig30 = figure('Name','Sheet SOM', 'Color','white');
% set(fig30, 'Position',  [50, 50, wd, ht])
% %subplot(4,5,20)
% som_grid(sMapr,'Coord',Cor,'Surf',sMapr_d.codebook(:,20),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('Irrad.')%, view(-100,0)%, 
% axis tight, axis equal
% axis off

%% Cylindrical SOM
% figure('Name','Cylindrical SOM', 'Color','white');
% 
% subplot(8,2,1);
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,9),'Markersize',0.5,'Line','none');
% colormap(jet);
% cb=colorbar('eastoutside'); %set(cb,'position',[.15 .1 .1 .3]);
% %set(ax1,'ColorScale','log');
% %cb.Ruler.Scale = 'log';
% title('Heat Duty (cal/sec)')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% 
% 
% subplot(8,2,2)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,1),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('IN Syngas-H_2 (kmol/h)')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,3)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,2),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('IN Syngas-CO (kmol/h)')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,4)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,3),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('IN Syngas-CH_4 (kmol/h)')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,5)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,4),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('IN Syngas-CO_2 (kmol/h)')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,6)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,5),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('IN Syngas-N_2 (kmol/h)')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,7)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,7),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('IN Syngas-H_2O (kmol/h)')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,8)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,6),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('IN O_2 (kmol/h)')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,9)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,8),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('T reactor (K)')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,10)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,10),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('OUT MoleFlow (kmol/h)')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,11)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,11),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('OUT CO-molFrac')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,12)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,12),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('OUT CH_4-molFrac')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,13)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,13),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('OUT CO_2-molFrac')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,14)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,14),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('OUT N_2-molFrac')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,15)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,15),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('OUT O_2-molFrac')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% subplot(8,2,16)
% som_grid(sMapc,'Coord',Coc,'Surf',sMapc_d.codebook(:,16),'Markersize',0.5,'Line','none');
% colormap(jet), colorbar
% title('OUT H_2O-molFrac')%, view(-100,0)%, 
% axis tight, axis equal
% axis off
% 
% 



%%
%
%som_show(sMap,'umat', 'all', 'empty', 'norm','d');
%som_show(sMap,'umat', 'all', 'comp','all','norm','d');
%
%    Notice that the names of the components are included as the
%    titles of the subplots. Notice also that the variable values
%    have been denormalized to the original range and scale.
%
%    The component planes ('PetalL', 'PetalW', 'SepalL' and 'SepalW')
%    show what kind of values the prototype vectors of the map units
%    have. The value is indicated with color, and the colorbar on the
%    right shows what the colors mean.
%
%    The 'U-matrix' shows distances between neighboring units and thus
%    visualizes the cluster structure of the map. Note that the
%    U-matrix visualization has much more hexagons that the
%    component planes. This is because distances *between* map units
%    are shown, and not only the distance values *at* the map units. 
%
%    High values on the U-matrix mean large distance between
%    neighboring map units, and thus indicate cluster
%    borders. Clusters are typically uniform areas of low
%    values. Refer to colorbar to see which colors mean high
%    values. In the IRIS map, there appear to be two clusters.
%
%pause % Strike any key to continue...
%
%    The subplots are linked together through similar position. In
%    each axis, a particular map unit is always in the same place. For
%    example:
%
%h=zeros(sMap.topol.msize); h(1,2) = 1;
%som_show_add('hit',h(:),'markercolor','r','markersize',0.5,'subplot','all')
%
%    the red marker is on top of the same unit on each axis. 

pause % Strike any key to continue...



%%

%clf

clc

%    STEP 4: VISUALIZING THE SELF-ORGANIZING MAP: SOM_SHOW_ADD
%    =========================================================

%    The SOM_SHOW_ADD function can be used to add markers, labels and
%    trajectories on top of SOM_SHOW created figures. The function
%    SOM_SHOW_CLEAR can be used to clear them away.

%    Here, the U-matrix is shown on the left, and an empty grid
%    named 'Labels' is shown on the right.

%som_show(sMap,'umat','all','empty','Labels')
%som_show(sMap,'umat','all')

%
%pause % Strike any key to add labels...
%
%    Here, the labels added to the map with SOM_AUTOLABEL function
%    are shown on the empty grid.
%
%som_show_add('label',sMap,'Textsize',8,'TextColor','r','Subplot',2)
%
%pause % Strike any key to add hits...
%
%    An important tool in data analysis using SOM are so called hit
%    histograms. They are formed by taking a data set, finding the BMU
%    of each data sample from the map, and increasing a counter in a
%    map unit each time it is the BMU. The hit histogram shows the
%    distribution of the data set on the map.
%
%    Here, the hit histogram for the whole data set is calculated
%    and visualized on the U-matrix.
%
%h = som_hits(sMap,sDiris);
%som_show_add('hit',h,'MarkerColor','w','Subplot',1)
%
%pause % Strike any key to continue...
%
%    Multiple hit histograms can be shown simultaniously. Here, three
%    hit histograms corresponding to the three species of Iris
%    flowers is calculated and shown. 
%
%    First, the old hit histogram is removed.
%
%som_show_clear('hit',1)
%
%    Then, the histograms are calculated. The first 50 samples in
%    the data set are of the 'Setosa' species, the next 50 samples
%    of the 'Versicolor' species and the last 50 samples of the
%    'Virginica' species. 
%
%h1 = som_hits(sMap,sDiris.data(1:50,:));
%h2 = som_hits(sMap,sDiris.data(51:100,:));
%h3 = som_hits(sMap,sDiris.data(101:150,:));
%
%som_show_add('hit',[h1, h2, h3],'MarkerColor',[1 0 0; 0 1 0; 0 0 1],'Subplot',1)
%
%    Red color is for 'Setosa', green for 'Versicolor' and blue for
%    'Virginica'. One can see that the three species are pretty well
%    separated, although 'Versicolor' and 'Virginica' are slightly
%    mixed up.
%
pause % Strike any key to continue...



%%

%clf

clc

%    STEP 4: VISUALIZING THE SELF-ORGANIZING MAP: SOM_GRID
%    =====================================================

%    There's also another visualization function: SOM_GRID.  This
%    allows visualization of the SOM in freely specified coordinates,
%    for example the input space (of course, only upto 3D space). This
%    function has quite a lot of options, and is pretty flexible.

%    Basically, the SOM_GRID visualizes the SOM network: each unit is
%    shown with a marker and connected to its neighbors with lines.
%    The user has control over: 
%     - the coordinate of each unit (2D or 3D)
%     - the marker type, color and size of each unit
%     - the linetype, color and width of the connecting lines
%    There are also some other options.

%pause % Strike any key to see some visualizations...

%    Here are four visualizations made with SOM_GRID: 
%     - The map grid in the output space.

%subplot(2,2,1)
%som_grid(sMap,'Linecolor','k')
%view(0,-90), title('Map grid')

%     - A surface plot of distance matrix: both color and 
%       z-coordinate indicate average distance to neighboring 
%       map units. This is closely related to the U-matrix.

%subplot(2,2,2)
%Co=som_unit_coords(sMap); U=som_umat(sMap); U=U(1:2:size(U,1),1:2:size(U,2));
%som_grid(sMap,'Coord',[Co, U(:)],'Surf',U(:),'Marker','none');
%view(-80,45), axis tight, title('Distance matrix')

%     - The map grid in the output space. Three first components
%       determine the 3D-coordinates of the map unit, and the size
%       of the marker is determined by the fourth component.
%       Note that the values have been denormalized.

%subplot(2,2,3)
%M = som_denormalize(sMap.codebook,sMap);
%som_grid(sMap,'Coord',M(:,1:3),'MarkerSize',M(:,4)*2)
%view(-80,45), axis tight, title('Prototypes')

%     - Map grid as above, but the original data has been plotted
%       also: coordinates show the values of three first components
%       and color indicates the species of each sample.  Fourth
%       component is not shown.

%subplot(2,2,4)
%som_grid(sMap,'Coord',M(:,1:3),'MarkerSize',M(:,4)*2)
%hold on
%D = som_denormalize(sDiris.data,sDiris); 
%plot3(D(1:50,1),D(1:50,2),D(1:50,3),'r.',...
%      D(51:100,1),D(51:100,2),D(51:100,3),'g.',...
%      D(101:150,1),D(101:150,2),D(101:150,3),'b.')
%view(-72,64), axis tight, title('Prototypes and data')

pause % Strike any key to continue...



%%

%    STEP 5: ANALYSIS OF RESULTS
%    ===========================

%    The purpose of this step highly depends on the purpose of the
%    whole data analysis: is it segmentation, modeling, novelty
%    detection, classification, or something else? For this reason, 
%    there is not a single general-purpose analysis function, but 
%    a number of individual functions which may, or may not, prove 
%    useful in any specific case.

%    Visualization is of course part of the analysis of
%    results. Examination of labels and hit histograms is another
%    part. Yet another is validation of the quality of the SOM (see
%    the use of SOM_QUALITY in SOM_DEMO1).

[qe_r,te_r] = som_quality(sMapr,sDdata)
%[qe_c,te_c] = som_quality(sMapc,sDiris)

%    People have contributed a number of functions to the Toolbox
%    which can be used for the analysis. These include functions for 
%    vector projection, clustering, pdf-estimation, modeling,
%    classification, etc. However, ultimately the use of these
%    tools is up to you.

%    More about visualization is presented in SOM_DEMO3.
%    More about data analysis is presented in SOM_DEMO4.

pause % Strike any key to continue...


%    CLUSTERING OF THE MAP
%    =====================

%f2=figure;

clc

%    Visual inspection already hinted that there are at least two
%    clusters in the data, and that the properties of the clusters are
%    different from each other (esp. relation of SepalL and
%    SepalW). For further investigation, the map needs to be
%    partitioned.

%    Here, the KMEANS_CLUSTERS function is used to find an initial
%    partitioning. The plot shows the Davies-Boulding clustering
%    index, which is minimized with best clustering.


%subplot(1,3,1)
%[c,p,err,ind] = kmeans_clusters(sMapr, 5); % find at most 7 clusters
%[c,p,err,ind] = kmeans_clusters(sMapc, 50);
%plot(1:length(ind),ind,'x-');
%[dummy,i] = min(ind);
%cl = p{i};

% figure('Name','SOM Clustering - Index', 'Color','white');
%     [c, p, err, ind] = kmeans_clusters(sMapr,10); % find clusterings
%     
%     plot(1:length(ind),ind,'x-');
%     [dummy,i] = min(ind); % select the one with smallest index
%     cl = p{i};
% figure('Name','SOM Clustering - Mapped', 'Color','white');
%     
%     som_show(sMapr,'color',{p{i},sprintf('%d clusters',i)}); % visualize
%     colormap(jet(i)), som_recolorbar % change colormap
  

%    The Davies-Boulding index seems to indicate that there are
%    two clusters on the map. Here is the clustering info
%    calculated previously and the partitioning result: 

%subplot(1,3,2)
%som_cplane(sM,Code,Dm)
%subplot(1,3,3)
%som_cplane(sM,cl)

%    You could use also function SOM_SELECT to manually make or modify
%    the partitioning.

%    After this, the analysis would proceed with summarization of the
%    results, and analysis of each cluster one at a time.
%    Unfortunately, you have to do that yourself. The SOM Toolbox does
%    not, yet, have functions for those purposes.

pause % Strike any key to continue...

echo off
warning on

