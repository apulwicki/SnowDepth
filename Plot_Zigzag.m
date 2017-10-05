%% Plotting zigzag data - all zigs for each G

load('ZigZagSWE.mat')

FS = 16; %fontsize
MinMaxC = [0 220]; ABC = ['a','b','c','d'];

for i = 1:3
   labels = categories(SWE(i).ZZ);
       pointsize = 30; clear s

   if length(labels) == 3
      for j = 1:length(labels)
          thatZZ = SWE(i).ZZ == labels(j);
          x = SWE(i).utm(thatZZ,1) - min(SWE(i).utm(thatZZ,1));     
          y = SWE(i).utm(thatZZ,2) - min(SWE(i).utm(thatZZ,2));   
          z = SWE(i).depth(thatZZ,1);
          s(j) = subplot(1,3,j);
            scatter(x, y, pointsize, z, 'filled');
                axis([0 45 0 45]); axis square; box on
                xlabel('Distance (m)'); ylabel('Distance (m)'); title(char(labels(j)));
                if j == length(labels); 
                caxis(MinMaxC);     c = colorbar;   c.Label.String = 'Snow depth (cm)';  %colormap(flipud(cool))
                end
            text(41, 42, ABC(j),'FontSize',FS,'Interpreter','Latex')
            text(5, 3, ['$$\bar{z}=$$',num2str(round(mean(z),0)),' cm'],'FontSize',FS,'Interpreter','Latex')
      end
      s1Pos = get(s(1),'position');   s3Pos = get(s(end),'position');   
      s3Pos(3:4) = s1Pos(3:4);        set(s(end),'position',s3Pos);
      fig=gcf;  fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
   
   elseif length(labels) == 4
      for j = 1:length(labels)
          thatZZ = SWE(i).ZZ == labels(j);
          x = SWE(i).utm(thatZZ,1) - min(SWE(i).utm(thatZZ,1));     
          y = SWE(i).utm(thatZZ,2) - min(SWE(i).utm(thatZZ,2));   
          z = SWE(i).depth(thatZZ,1);
          s(j) = subplot(2,2,j);
            scatter(x, y, pointsize, z, 'filled');
                axis([0 45 0 45]); axis square; box on
                xlabel('Distance (m)'); ylabel('Distance (m)'); title(char(labels(j)));
                if j == length(labels); 
                caxis(MinMaxC);     c = colorbar;   c.Label.String = 'Snow depth (cm)';  %colormap(flipud(cool))
                end
             text(41, 42, ABC(j),'FontSize',FS,'Interpreter','Latex')   
             text(5, 3, ['$$\bar{z}=$$',num2str(round(mean(z),0)),' cm'],'FontSize',FS,'Interpreter','Latex')

      end
      s1Pos = get(s(1),'position');   s3Pos = get(s(end),'position');   
      s3Pos(3:4) = s1Pos(3:4);        set(s(end),'position',s3Pos);
      fig=gcf;  fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 10 8];
   end 

   set(findall(fig,'-property','FontSize'),'FontSize',FS)
saveFIG(['ZigzagDepth_', char(SWE(i).glacier(1,1))])         
end

%dlmwrite('/home/glaciology1/Documents/QGIS/Data/TestPoints.csv',closest, 'delimiter', ',', 'precision', 9); %Write matrix with new waypoints to csv file for QGIS


%% Plotting zigzag data - compare location calc options
clf
    run OPTIONS.m
    options.ZZ = 3; %only zigzags
    options.ZigzagLocation  = 1; %location based on vertex
    run MAIN
    i = 3;
thatZZ = SWE(i).ZZ == SWE(i).ZZ(end,1);
minX = min(SWE(i).utm(thatZZ,1)); minY = min(SWE(i).utm(thatZZ,2));

x = SWE(i).utm(thatZZ,1) - minX;     
y = SWE(i).utm(thatZZ,2) - minY;   
plot(x, y,'.k','MarkerSize',10); hold on

    run OPTIONS.m
    options.ZZ = 3; %only zigzags
    options.ZigzagLocation  = 2; %location based on last point
    run MAIN
    i = 3;
x = SWE(i).utm(thatZZ,1) - minX;     
y = SWE(i).utm(thatZZ,2) - minY;   
plot(x, y,'ob','MarkerSize',6); hold on

	axis([0 45 0 45]); axis square; box on; ax = gca;
    ax.YTick = 0:10:40;
    xlabel('Distance (m)'); ylabel('Distance (m)'); 

fig=gcf;  fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 6 6];

      set(findall(fig,'-property','FontSize'),'FontSize',16)
filename = 'Zigzag_calOptions';
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')         
 

%% Boxplot to compare zigzags
% run OPTIONS.m
% options.ZZ = 3; %only zigzags
% run MAIN

clf
yaxis = [min([SWE(1).depth; SWE(2).depth; SWE(3).depth]) max([SWE(1).depth; SWE(2).depth; SWE(3).depth])];
text_dim = [.157 .61 .3 .3;   .446 .61 .3 .3;     .726 .61 .3 .3];

for i = 1:3
   s(i) = subplot(1,3,i);
   boxplot(SWE(i).depth,SWE(i).ZZ)
      ylabel('Depth (cm)');     ylim(yaxis)
      set(gca,'FontSize',10,'XTickLabelRotation',90)
      annotation('textbox',text_dim(i,:),'String',char(options.glacier(i)),'FitBoxToText','on','EdgeColor',[1 1 1]);
end
    fig=gcf;  fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 15 5];
      set(findall(fig,'-property','FontSize'),'FontSize',16)
saveFIG('Zigzag_Boxplot')
%%
% Glacier maps overlay
subplot(1,3,1)
    mapG4 = imread('/home/glaciology1/Documents/MastersDocuments/Methods/map_zigzaglocation_G4.jpeg');
     Nshow = imshow(mapG4);
subplot(1,3,2)
    mapG2 = imread('/home/glaciology1/Documents/MastersDocuments/Methods/map_zigzaglocation_G2.jpeg');
    Nshow = imshow(mapG2);
subplot(1,3,3)   
    mapG13 = imread('/home/glaciology1/Documents/MastersDocuments/Methods/map_zigzaglocation_G13.jpeg');
    Nshow = imshow(mapG13); 

      fig=gcf;  fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 15 5];
      set(findall(fig,'-property','FontSize'),'FontSize',16)
saveFIG('Zigzag_Location_Map')         

        
%% Length of measured vs vertex spacing

section_sum = [1];
for i = 1:size(ZZ.depth,1)-1
    if ~cellfun(@isequal,ZZ.depth(i,1),ZZ.depth(i+1,1))
        section_sum = [section_sum ; i+1];
    end
end
for i = 1:size(section_sum,1)-1
    section_sum(i,2) = cell2mat(ZZ.depth(section_sum(i+1,1)-1,2));
end
section_sum(end,:) = [];
display(strcat('Section sum = ', num2str(mean(section_sum(:,2)))))

vertex_spacing = [];
for i = 2:8
    vertex_spacing(i,1) = EuclideanDistance(cell2mat(ZZ.vertexcoord(i,1)),cell2mat(ZZ.vertexcoord(i,2)),...
        cell2mat(ZZ.vertexcoord(i+1,1)),cell2mat(ZZ.vertexcoord(i+1,2)));
end
vertex_spacing(1,:) = [];
vertex_spacing(4,:) = [];
display(strcat('Vertex spacing = ',num2str(mean(vertex_spacing(:,1)))))

plot(section_sum(:,2),'o')
hold on
plot([0 size(section_sum,1)],[mean(vertex_spacing(:,1)) mean(vertex_spacing(:,1))])
    ylabel('Distance (m)')
    legend('Measured spacing','Vertex spacing', 'Location','best')

%% Probably density function

load TopoSWE.mat

%Histogram with individal zigzags divided into glaciers
figure(3)
clf; clc
for i = 1:3
        name    = char(options.glacier(i)); 
        zz      = categories(SWEzz(i).ZZ);
        
        subplot(3,1,i)
        for j = 1:length(zz)
            ZZdata = SWEzz(i).swe(SWEzz(i).ZZ==char(zz(j)));
            ZZdata = (ZZdata-mean(ZZdata))/mean(ZZdata)*100;%/std(ZZdata);
            
            display([zz(j), num2str(chi2gof(ZZdata)), num2str(std(SWEzz(i).swe(SWEzz(i).ZZ==char(zz(j)))))])
                if j ==1;   bins    = round(sqrt(length(SWEzz(i).swe)));
                            edges   = linspace(-50,50,bins); end
                            N = histcounts(ZZdata,edges);   
            plot((edges(:,1:end-1)+edges(:,2:end))/2,N/sum(N),'LineWidth',2); hold on 
            xlabel('SWE variability (%)');     ylabel('Probability')
            title(name)
            ax = gca; ax.XTick = [-40:20:40];
        end
            legend(zz)
 %            xlim([-6 6]);         ylim([0 60]);
end
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13) 
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 8 12];
saveFIG('ZigzagPDF');

%% Probability Density Function -> Fitted
figure
%G4
data = SWEzz(1).swe;
chi2gof(data(SWEzz(1).ZZ=='G04\_Z2A'))
chi2gof(data(SWEzz(1).ZZ=='G04\_Z3A'))
chi2gof(data(SWEzz(1).ZZ=='G04\_Z5B'))
% 
% [d pd] = allfitdist(data(SWEzz(1).ZZ=='G04\_Z2A')); d(1)
% [d pd] = allfitdist(data(SWEzz(1).ZZ=='G04\_Z3A')); d(1)
% [d pd] = allfitdist(data(SWEzz(1).ZZ=='G04\_Z5B')); d(1)

[pdca,gn,~] = fitdist(data,'Normal','By',SWEzz(1).ZZ);
Z2A = pdca{1};      Z3A = pdca{2};      Z5B = pdca{3};
x_values = 0:0.01:1;
Z2Apdf = pdf(Z2A,x_values);     Z3Apdf = pdf(Z3A,x_values);     Z5Bpdf = pdf(Z5B,x_values);
    subplot(3,1,1)
    plot(x_values,Z2Apdf,'LineWidth',2); hold on
    plot(x_values,Z3Apdf,'LineWidth',2); hold on
    plot(x_values,Z5Bpdf,'LineWidth',2);
    legend(gn,'Location','NorthEast')
    title('G4'); xlabel('SWE (m)'); ylabel('Probability density (m^{-1})');

%G2
data = SWEzz(2).swe;
chi2gof(data(SWEzz(2).ZZ=='G02\_Z3B'))
chi2gof(data(SWEzz(2).ZZ=='G02\_Z5C'))
chi2gof(data(SWEzz(2).ZZ=='G02\_Z7A'))

% [d, pd] = allfitdist(data(SWEzz(2).ZZ=='G02\_Z3B')); d(1)


[pdca,gn,~] = fitdist(data,'Normal','By',SWEzz(2).ZZ);
x_values = 0:0.01:1;
Z3B = pdca{1};                  Z5C = pdca{2};                  Z7A = pdca{3};
Z3Bpdf = pdf(Z3B,x_values);     Z5Cpdf = pdf(Z5C,x_values);     Z7Apdf = pdf(Z7A,x_values);
    subplot(3,1,2)
    plot(x_values,Z3Bpdf,'LineWidth',2); hold on
    plot(x_values,Z5Cpdf,'LineWidth',2); hold on
    plot(x_values,Z7Apdf,'LineWidth',2);
    legend(gn,'Location','NorthEast')
    title('G2')
    
%G13
data = SWEzz(3).swe;
    chi2gof(data(SWEzz(3).ZZ=='G13\_Z3B'))
    chi2gof(data(SWEzz(3).ZZ=='G13\_Z4C'))
    chi2gof(data(SWEzz(3).ZZ=='G13\_Z5A'))
    chi2gof(data(SWEzz(3).ZZ=='G13\_Z7C'))

data = SWEzz(3).swe;
[pdca,gn,~] = fitdist(data,'Normal','By',SWEzz(3).ZZ);
x_values = 0:0.01:1;
Z3B = pdca{1};              Z4C = pdca{2};              Z5A = pdca{3};              Z7C = pdca{4};
Z3Bpdf = pdf(Z3B,x_values); Z4Cpdf = pdf(Z4C,x_values);	Z5Apdf = pdf(Z5A,x_values); Z7Cpdf = pdf(Z7C,x_values); 
    subplot(3,1,3)
    plot(x_values,Z3Bpdf,'LineWidth',2); hold on
    plot(x_values,Z4Cpdf,'LineWidth',2); hold on
    plot(x_values,Z5Apdf,'LineWidth',2); hold on
    plot(x_values,Z7Cpdf,'LineWidth',2);
    legend(gn,'Location','NorthEast')
    title('G13')

        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13) 
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 10 10];

%% Chi squared normal distribution test
figure(1); clf
figure(2); clf
 %Prob that normal and ZZ dist are different (>0.05 means ZZ data is
 %normal)
clc
for g = 1:3
    allZZ = categories(SWEzz(g).ZZ);
    for z = 1:length(allZZ)
        tempZZ = allZZ(z);
        data = SWEzz(g).swe(SWEzz(g).ZZ == tempZZ);
             cv = std(data)/mean(data);
        data = (data-mean(data));%/std(data);
        [~,chiP] = chi2gof(data);
           
        display([tempZZ, num2str(chiP), num2str(cv)])
    end
end

 % Glacier wide ZZ data compiled
figure(1); clf
fullZZ = []; 
for g = 1:3
    allZZ = categories(SWEzz(g).ZZ);
    gZZ = [];
    for z = 1:length(allZZ)
        tempZZ = allZZ(z);
        if ~strcmp(tempZZ,'G02\_Z7A')
        data = SWEzz(g).swe(SWEzz(g).ZZ == tempZZ);
        data = (data-mean(data));%/mean(data)*100;%/std(data);
             
        gZZ = [gZZ; data];
        fullZZ = [fullZZ; data];
        end
    end
    
    glacier = char(options.glacier(g));
%     [~,chiP] = chi2gof(gZZ);
%     display([glacier, ' ', num2str(chiP)])
        display([glacier,' 2 sigma ', num2str(2*nanstd(gZZ))]);
if g ==1;     bins    = round(sqrt(length(gZZ))); end   
   edges   = linspace(-0.2,0.2,bins);
   N       =  histcounts(data,edges);
figure(1);   
 plot([min(edges) (edges(:,1:end-1)+edges(:,2:end))/2 max(edges)],[0 N/sum(N) 0],...
       'Color',options.RGB(g,:),'LineWidth',3); hold on
%      figure(2);
%     normALL = normpdf(gZZ,mean(gZZ),std(gZZ));
%     [gZZ, I] = sort(gZZ);
%     plot(gZZ,normALL(I),'LineWidth',2,'Color',options.RGB(g,:)); hold on
end
    [~,chiP] = chi2gof(fullZZ);
    %display(['All zigzags ', num2str(chiP)])
     figure(1)
%      N       =  histcounts(fullZZ,edges);
% figure(1);   plot((edges(:,1:end-1)+edges(:,2:end))/2,N/sum(N),'k--','LineWidth',2); hold on 
          xlabel('Winter balance (m w.e.)');     ylabel('Probability')
%         xlim([-6 6]);         
            
        legend('Glacier 4','Glacier 2','Glacier 13')
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18) 
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 10 6];
saveFIG('ZigzagPDF_G');

 
 %Normal Dist
  figure(2);
normALL = normpdf(fullZZ,mean(fullZZ),std(fullZZ));
[fullZZ, I] = sort(fullZZ);
plot(fullZZ,normALL(I),'LineWidth',2,'Color','k','LineStyle','--')
        legend([options.glacier,{'All zigzags'}]);     
        xlabel('SWE about local mean (m w.e.)');     ylabel('Freqency')
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13) 
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 11 8];

 
 %Confidence Interval
[mu,s,muci,sci] = normfit(fullZZ);
muci

%% Kruskal Wallis Test - between zigzags

clc; close all

for g = 1:3
    allZZ = categories(SWEzz(g).ZZ);
    for z = 1:length(allZZ)
        tempZZ = allZZ(z); I = SWEzz(g).ZZ == tempZZ;
        D = (SWEzz(g).swe(I)-mean(SWEzz(g).swe(I)))/mean(SWEzz(g).swe(I));
        data = [data, [D; nan(S-length(D),1)]];
    end
end
data = data(:,2:end);
%data(:,6) = [];

[kw, ~, stats] = kruskalwallis(data)
multcompare(stats)

%% How many zigzag points do you need? -> % error from mean
clf
S = 200;
data = nan(S,1);
for g = 1:3
clear P Esub 
    allZZ = categories(SWEzz(g).ZZ);
    for z = 1:length(allZZ)
        tempZZ = allZZ(z); I = SWEzz(g).ZZ == tempZZ;
        Dfull = SWEzz(g).swe(I);
        Dmean = mean(Dfull);
    
     %pull subset
    N = 1:2:40;
    for i = 1:length(N)
        for rep = 1:1000;
        r = randi([i length(Dfull)],1, N(i));
        Esub(i,rep) = abs(mean(Dfull(r))-Dmean)/Dmean*100;
        end
    end
subplot(1,3,g)
    P(z) = plot(N,mean(Esub,2),'LineWidth', 3); hold on
    plot(N,mean(Esub,2)+2*std(Esub,[],2), 'Color', P(z).Color); hold on
    end
    plot([0 max(N)],[5 5],'k--')
    legend(P, allZZ)
    ylabel('Error (%)');    xlabel('Sample size');
end


