function [  ] = saveFIG( filename, size )

global options


if nargin == 1
    size = [];
elseif nargin == 2  
    if strcmp(size,'3G')
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 6];
    end
end    
    

print([options.path1, filename],'-dpng','-r0'); 
print([options.path2, filename],'-dpng','-r0');  
print([options.path3, filename],'-dpng','-r0');  

end

