function [  ] = saveFIG( filename, fontN, size )

global options

 fig=gcf;    

if nargin == 1
    fontN = 18;
elseif nargin == 3  
    if strcmp(size,'3G')
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 6];
    end
end    
    
set(findall(fig,'-property','FontSize'),'FontSize',fontN)

print([options.path1, filename],'-dpng','-r0'); 
print([options.path2, filename],'-dpng','-r0');  
print([options.path3, filename],'-dpng','-r0');  

end

