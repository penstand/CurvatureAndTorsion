%% Loading data
pixelSize = 11/64; %micron
frameRate = 1000;

cellPath = './Human/20170104-Normal-510nm-1KHz-64x/Exp10Cell1/';
load(fullfile(cellPath, 'xArc.txt'))
load(fullfile(cellPath, 'yArc.txt'))
load(fullfile(cellPath, 'zArc.txt'))

range = 1:600;% frame range to be analyzed
%% convert the unit from pixel to micron

x = xArc(range,:)*pixelSize;
y = yArc(range,:)*pixelSize;
z = zArc(range,:)*pixelSize;

[m,n] = size(x);

kappa = nan(m,n);
tau = nan(m,n);
%% here is for smoothing
fX = nan(m,n);
fY = nan(m,n);
fZ = nan(m,n);
for i=1:m
    fX(i,:) = powersmooth2(x(i,:),2,1e4);% you can also try other smoothing methods such as 'rlowess'
    fY(i,:) = powersmooth2(y(i,:),2,1e4);
    fZ(i,:) = powersmooth2(z(i,:),2,2e4);
end

%% curvature and torsion from 'frenet_robust'method
for i=1:m
    eInd=n;
    for j=n:-1:1
        if ~isnan(fX(i,j))
            eInd=j;
            break;
        end
    end
    tX=fX(i,1:eInd);
    tY=fY(i,1:eInd);
    tZ=fZ(i,1:eInd);
    
    winLen = 30; 
    weight = 0.1;
    [ka,ta]=frenet_robust([tX;tY;tZ],winLen,weight);

    kappa(i,1:eInd)=ka;
    tau(i,1:eInd)=ta;
end

%% fix reverse sign of curvature
skip_end=50;
nframe = m;
for iframe=2:nframe
    ind=isfinite(kappa(iframe,:)) & isfinite(kappa(iframe-1,:));
    ind=ind & ( (1:n)>skip_end ) & ( (1:n)<n-skip_end+1 );
    if norm(kappa(iframe,ind)-kappa(iframe-1,ind)) > norm(kappa(iframe,ind)+kappa(iframe-1,ind))
    %if kappa(iframe,ind).*kappa(iframe-1,ind)<0
        kappa(iframe,:)= -kappa(iframe,:);
        if norm(kappa(iframe,ind)+kappa(iframe-1,ind)) < 0.5 
          warning(sprintf('problem at frame %d\n',iframe));
        end
    end
end

%% plot curvature and torsion kymogrpah
figure
imagesc(kappa')
%colormap(othercolor('RdYlBu8')) %one need the "othercolor" package to use
%...the color LUT 'RdYlBu8'
len10um = 10/pixelSize;
xlim([0, 600])
yticks([1:5]*len10um)
yticklabels([10 20 30 40 50])
colorbar
caxis([-0.4,0.4])

figure
imagesc(tau')
%colormap(othercolor('RdYlBu8'))
yticks([1:5]*len10um)
yticklabels([10 20 30 40 50])
colorbar
caxis([-0.5,0.5])
xlim([0, 600])


