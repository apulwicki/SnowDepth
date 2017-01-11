
%2D derivative
h = 40; % step size
f = topo_full_ns.G2.elevation;

order = 1;
Yn = diff(f,order,1)/h;
Ye = diff(f,order,2)/h;

figure(1)
    subplot(1,3,1)
    imagesc(f); colorbar
    
    subplot(1,3,2)
    imagesc(Yn); colorbar
    
    subplot(1,3,3)
    imagesc(Ye); colorbar
    
    
%% Central Difference

h = 40; % step size
f = topo_full_ns.G2.elevation;

%Smoothing grid
sizeG = 11; C = (sizeG+1)/2;
grid = ones(C);
for i = 1:C;
    for j = 1:C
        grid(i,j) = sqrt(((C-i)*h)^2+((C-j)*h)^2);
    end
end
grid = [grid, fliplr(grid(:,1:C-1))];     grid = [grid; flipud(grid(1:C-1,:))];
grid = 1./grid;     grid(C,C) = 1;
grid = grid/(sum(grid(:)));

% f_filt = filter2(grid,f);
 f_filt = imgaussfilt(f,0.75);

% f_filt = nan(size(f));
% for i = C:size(f,1)-C
%    for j =  C:size(f,2)-C
%        T    = grid.*f(i-(C-1):i+(C-1),j-(C-1):j+(C-1));
%        f_filt(i,j) = nansum(T(:));
%    end
% end
% f_filt(f_filt<min(f(:))) = NaN;

%Stencil N = 3
    Fmn = f_filt(1:end-2,:);  Fme = f_filt(:,1:end-2);
    F0n = f_filt(2:end-1,:);  F0e = f_filt(:,2:end-1);  
    Fpn = f_filt(3:end,:);    Fpe = f_filt(:,3:end);
Y3n_1 = (Fpn-Fmn)/(2*h);            Y3e_1 = (Fpe-Fme)/(2*h);
Y3n_2 = (Fmn-2*F0n+Fpn)/h^2;        Y3e_2 = (Fme-2*F0e+Fpe)/h^2;

figure(2); clf
    subplot(3,2,1)
    imagesc(f); colorbar
        title('Original')
    
    subplot(3,2,2)
        imagesc(f_filt); colorbar; caxis([1900 3100])
        title(['Smoothed, window size = ',num2str(sizeG)]);
    
    subplot(3,2,3)
        imagesc(Y3n_1); colorbar; caxis([-1 1.5])
        title('1st Derivative - N')
    
    subplot(3,2,4)
        imagesc(Y3e_1); colorbar; caxis([-1 1.5])
        title('1st Derivative - E')

    subplot(3,2,5)
        imagesc(Y3n_2); colorbar; caxis([-0.01 0.01])
        title('2nd Derivative - N')
    
    subplot(3,2,6)
        imagesc(Y3e_2); colorbar; caxis([-0.01 0.01])
        title('2nd Derivative - E')
 
figure(1)
    imagesc((Y3n_2(:,3:end)+Y3e_2(3:end,:))/2); colorbar; caxis([-0.01 0.01])

        
% %Stencil N = 5
%     Fm1n = f(2:end-3,:);  Fm1e = f(:,2:end-3);
%     Fp1n = f(4:end-1,:);  Fp1e = f(:,4:end-1);
%     Fm2n = f(1:end-4,:);  Fm2e = f(:,1:end-4);
%     Fp2n = f(5:end,:);    Fp2e = f(:,5:end);
% Y5n = (Fm2n-8*Fm1n+8*Fp1n-Fp2n)/(12*h);    Y5e = (Fm2e-8*Fm1e+8*Fp1e-Fp2e)/(12*h);
% 
% 
% figure(3); clf
%     subplot(1,3,1)
%     imagesc(f); colorbar
%     
%     subplot(1,3,2)
%     imagesc(Y5n); colorbar
%     
%     subplot(1,3,3)
%     imagesc(Y5e); colorbar    
%     
    
    
    
    
    