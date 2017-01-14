function [  ] = saveFIG( filename )

global options

print([options.path1, filename],'-dpng','-r0'); 
print([options.path2, filename],'-dpng','-r0');  
print([options.path3, filename],'-dpng','-r0');  

end

