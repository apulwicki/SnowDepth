function [  ] = saveFIG_IGS( filename, Ncol, height)

global options

 fig=gcf;    

 if     Ncol == 1
     width = 8.6;
 elseif Ncol == 2
     width = 17.8;
 end

fig.PaperUnits = 'centimeters';
fig.PaperSize  = [width height];
    
set(findall(fig,'-property','FontSize'),'FontSize',9)
set(findall(fig,'-property','FontName'),'FontName','Arial')

%print(['/home/glaciology1/Documents/MastersDocuments/Paper I/', filename],'-dtiff'); 
print(['/home/glaciology1/Documents/MastersDocuments/Paper I/', filename],'-dpdf','-fillpage'); 
%print(['/Users/Alexandra/Documents/SFU/MastersDocuments/Paper I/', filename],'-dpdf','-fillpage'); 
% print([options.path1, filename],'-dpdf','-fillpage'); 
% print([options.path2, filename],'-dpdf','-fillpage');  
% print([options.path3, filename],'-dpdf','-fillpage');  

end