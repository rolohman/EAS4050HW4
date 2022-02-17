
%==================================================

clear all

load arrows.scec3;
load coast;
load faults;

lat1 = arrows(:,1);
lon1 = arrows(:,2); %input 0->360 longitude


for i = 1:length(lon1)
  [E(i),N(i)]=utm2ll3(lon1(i),lat1(i),11,1); %zone 11, wgs84
end

for i = 1:length(coast_lon)
  [coastE(i),coastN(i)]=utm2ll3(coast_lon(i),coast_lat(i),11,1); %zone 11, wgs84
end
for i = 1:length(faults_lon)
  [faultsE(i),faultsN(i)]=utm2ll3(faults_lon(i),faults_lat(i),11,1); %zone 11, wgs84
end

UE = arrows(:,3)/1e3; %from mm to m
UN = arrows(:,4)/1e3; %from mm to m


dx=100e3;
gridy = 3.384e+06:dx:4.100e+06;
gridx = 1.989e+04:dx:8.840e+05;

nx=size(gridx,2);
ny=size(gridy,2);

[x,y] = meshgrid(gridx,gridy);

alpha = 50e3;

m1 = [];

for i=1:nx*ny
    W   = exp(-((E'-x(i)).^2 +(N'-y(i)).^2)/(2*(alpha^2)));
    ind = find(sqrt(((E'-x(i)).^2 +(N'-y(i)).^2)) < 2*alpha);
    if (size(ind,1) < 4)
        m1 = [m1; NaN(1,8)];
    else
        sig =  diag([W(ind); W(ind)]);
        d =sig*[UE(ind); UN(ind)];
        G = sig*[ones(size(UE(ind))) zeros(size(UE(ind))) E(ind)'-x(i) N(ind)'-y(i) zeros(size(UE(ind)))  zeros(size(UE(ind))); zeros(size(UE(ind))) ones(size(UE(ind))) zeros(size(UE(ind))) zeros(size(UE(ind))) E(ind)'-x(i) N(ind)'-y(i)];
        m2 = inv(G'*sig*G)*G'*sig*d;
        vgradi  = [m2(3) m2(4); m2(5) m2(6)];
        straini = 0.5*(vgradi + vgradi');
        rotati  = 0.5*(vgradi - vgradi');
        dA = trace(straini)/2;
        dS = sqrt(straini(1,2).^2 + (straini(1,1)-straini(2,2)).^2/4);
        [v d] = eig(straini);
        phi = atan2(v(2),v(1))*180/pi;
        m1 = [m1; dA dS rotati(3) phi straini(1) straini(3) straini(2) straini(4)];
    end
end

ext1 = reshape(m1(:,1),size(gridy,2),size(gridx,2));
shear = reshape(m1(:,2),size(gridy,2),size(gridx,2));
rot = reshape(m1(:,3),size(gridy,2),size(gridx,2));
ang = reshape(m1(:,4),size(gridy,2),size(gridx,2));
exx = reshape(m1(:,5),size(gridy,2),size(gridx,2));
eyy = reshape(m1(:,8),size(gridy,2),size(gridx,2));
exy = reshape(m1(:,6),size(gridy,2),size(gridx,2));


figure
plot(E/1e3,N/1e3,'k.')
hold on
plot(coastE/1e3,coastN/1e3,'k')
plot(faultsE/1e3,faultsN/1e3,'r')
set(gca,'xtick',[],'ytick',[],'position',[0 0 1 1])
axis([100 750 3500 4050])
print('-dpng','basic.png')



figure
pf=pcolor(gridx/1e3,gridy/1e3,shear);
pf.CData=NaN(size(x));
pf.EdgeColor='b';
hold on
shading faceted
plot(E/1e3,N/1e3,'k.')
plot(coastE/1e3,coastN/1e3,'k')
plot(faultsE/1e3,faultsN/1e3,'r')
set(gca,'xtick',[],'ytick',[],'position',[0 0 1 1])
axis([100 750 3500 4050])
print('-dpng','grid.png')


figure
pcolor(gridx/1e3,gridy/1e3,shear)
hold on
shading flat
plot(E/1e3,N/1e3,'k.')
plot(coastE/1e3,coastN/1e3,'k')
plot(faultsE/1e3,faultsN/1e3,'r')
caxis([0 1e-7])
set(gca,'xtick',[],'ytick',[],'position',[0 0 1 1])
axis([100 750 3500 4050])
print('-dpng','shear.png')



dx=50e3;
gridy = 3.384e+06:dx:4.100e+06;
gridx = 1.989e+04:dx:8.840e+05;

nx=size(gridx,2);
ny=size(gridy,2);

[x,y] = meshgrid(gridx,gridy);

alpha = 25e3;

m1 = [];

for i=1:nx*ny
    W   = exp(-((E'-x(i)).^2 +(N'-y(i)).^2)/(2*(alpha^2)));
    ind = find(sqrt(((E'-x(i)).^2 +(N'-y(i)).^2)) < 2*alpha);
    if (size(ind,1) < 4)
        m1 = [m1; NaN(1,8)];
    else
        sig =  diag([W(ind); W(ind)]);
        d =sig*[UE(ind); UN(ind)];
        G = sig*[ones(size(UE(ind))) zeros(size(UE(ind))) E(ind)'-x(i) N(ind)'-y(i) zeros(size(UE(ind)))  zeros(size(UE(ind))); zeros(size(UE(ind))) ones(size(UE(ind))) zeros(size(UE(ind))) zeros(size(UE(ind))) E(ind)'-x(i) N(ind)'-y(i)];
        m2 = inv(G'*sig*G)*G'*sig*d;
        vgradi  = [m2(3) m2(4); m2(5) m2(6)];
        straini = 0.5*(vgradi + vgradi');
        rotati  = 0.5*(vgradi - vgradi');
        dA = trace(straini)/2;
        dS = sqrt(straini(1,2).^2 + (straini(1,1)-straini(2,2)).^2/4);
        [v d] = eig(straini);
        phi = atan2(v(2),v(1))*180/pi;
        m1 = [m1; dA dS rotati(3) phi straini(1) straini(3) straini(2) straini(4)];
    end
end


shear = reshape(m1(:,2),size(gridy,2),size(gridx,2));

figure
pcolor(gridx/1e3,gridy/1e3,shear)
hold on
shading flat
plot(E/1e3,N/1e3,'k.')
plot(coastE/1e3,coastN/1e3,'k')
plot(faultsE/1e3,faultsN/1e3,'r')
caxis([0 1e-7])
set(gca,'xtick',[],'ytick',[],'position',[0 0 1 1])
axis([100 750 3500 4050])
print('-dpng','shear_highres.png')



dx=10e3;
gridy = 3.384e+06:dx:4.100e+06;
gridx = 1.989e+04:dx:8.840e+05;

nx=size(gridx,2);
ny=size(gridy,2);

[x,y] = meshgrid(gridx,gridy);

alpha = 5e3;

m1 = [];

for i=1:nx*ny
    W   = exp(-((E'-x(i)).^2 +(N'-y(i)).^2)/(2*(alpha^2)));
    ind = find(sqrt(((E'-x(i)).^2 +(N'-y(i)).^2)) < 2*alpha);
    if (size(ind,1) < 4)
        m1 = [m1; NaN(1,8)];
    else
        sig =  diag([W(ind); W(ind)]);
        d =sig*[UE(ind); UN(ind)];
        G = sig*[ones(size(UE(ind))) zeros(size(UE(ind))) E(ind)'-x(i) N(ind)'-y(i) zeros(size(UE(ind)))  zeros(size(UE(ind))); zeros(size(UE(ind))) ones(size(UE(ind))) zeros(size(UE(ind))) zeros(size(UE(ind))) E(ind)'-x(i) N(ind)'-y(i)];
        m2 = inv(G'*sig*G)*G'*sig*d;
        vgradi  = [m2(3) m2(4); m2(5) m2(6)];
        straini = 0.5*(vgradi + vgradi');
        rotati  = 0.5*(vgradi - vgradi');
        dA = trace(straini)/2;
        dS = sqrt(straini(1,2).^2 + (straini(1,1)-straini(2,2)).^2/4);
        [v d] = eig(straini);
        phi = atan2(v(2),v(1))*180/pi;
        m1 = [m1; dA dS rotati(3) phi straini(1) straini(3) straini(2) straini(4)];
    end
end

shear = reshape(m1(:,2),size(gridy,2),size(gridx,2));


figure
pcolor(gridx/1e3,gridy/1e3,shear)
hold on
shading flat
plot(E/1e3,N/1e3,'k.')
plot(coastE/1e3,coastN/1e3,'k')
plot(faultsE/1e3,faultsN/1e3,'r')
caxis([0 1e-7])
set(gca,'xtick',[],'ytick',[],'position',[0 0 1 1])
axis([100 750 3500 4050])
print('-dpng','shear_higherres.png')




dx=10e3;
gridy = 3.384e+06:dx:4.100e+06;
gridx = 1.989e+04:dx:8.840e+05;

nx=size(gridx,2);
ny=size(gridy,2);

[x,y] = meshgrid(gridx,gridy);

alpha = 50e3;

m1 = [];

for i=1:nx*ny
    W   = exp(-((E'-x(i)).^2 +(N'-y(i)).^2)/(2*(alpha^2)));
    ind = find(sqrt(((E'-x(i)).^2 +(N'-y(i)).^2)) < 2*alpha);
    if (size(ind,1) < 4)
        m1 = [m1; NaN(1,8)];
    else
        sig =  diag([W(ind); W(ind)]);
        d =sig*[UE(ind); UN(ind)];
        G = sig*[ones(size(UE(ind))) zeros(size(UE(ind))) E(ind)'-x(i) N(ind)'-y(i) zeros(size(UE(ind)))  zeros(size(UE(ind))); zeros(size(UE(ind))) ones(size(UE(ind))) zeros(size(UE(ind))) zeros(size(UE(ind))) E(ind)'-x(i) N(ind)'-y(i)];
        m2 = inv(G'*sig*G)*G'*sig*d;
        vgradi  = [m2(3) m2(4); m2(5) m2(6)];
        straini = 0.5*(vgradi + vgradi');
        rotati  = 0.5*(vgradi - vgradi');
        dA = trace(straini)/2;
        dS = sqrt(straini(1,2).^2 + (straini(1,1)-straini(2,2)).^2/4);
        [v d] = eig(straini);
        phi = atan2(v(2),v(1))*180/pi;
        m1 = [m1; dA dS rotati(3) phi straini(1) straini(3) straini(2) straini(4)];
    end
end

shear = reshape(m1(:,2),size(gridy,2),size(gridx,2));


figure
pcolor(gridx/1e3,gridy/1e3,shear)
hold on
shading flat
plot(E/1e3,N/1e3,'k.')
plot(coastE/1e3,coastN/1e3,'k')
plot(faultsE/1e3,faultsN/1e3,'r')
caxis([0 1e-7])
set(gca,'xtick',[],'ytick',[],'position',[0 0 1 1])
axis([100 750 3500 4050])
print('-dpng','shear_higherres_smooth.png')






dx=10e3;
gridy = 3.384e+06:dx:4.100e+06;
gridx = 1.989e+04:dx:8.840e+05;

nx=size(gridx,2);
ny=size(gridy,2);

[x,y] = meshgrid(gridx,gridy);

alpha = 100e3;

m1 = [];

for i=1:nx*ny
    W   = exp(-((E'-x(i)).^2 +(N'-y(i)).^2)/(2*(alpha^2)));
    ind = find(sqrt(((E'-x(i)).^2 +(N'-y(i)).^2)) < 2*alpha);
    if (size(ind,1) < 4)
        m1 = [m1; NaN(1,8)];
    else
        sig =  diag([W(ind); W(ind)]);
        d =sig*[UE(ind); UN(ind)];
        G = sig*[ones(size(UE(ind))) zeros(size(UE(ind))) E(ind)'-x(i) N(ind)'-y(i) zeros(size(UE(ind)))  zeros(size(UE(ind))); zeros(size(UE(ind))) ones(size(UE(ind))) zeros(size(UE(ind))) zeros(size(UE(ind))) E(ind)'-x(i) N(ind)'-y(i)];
        m2 = inv(G'*sig*G)*G'*sig*d;
        vgradi  = [m2(3) m2(4); m2(5) m2(6)];
        straini = 0.5*(vgradi + vgradi');
        rotati  = 0.5*(vgradi - vgradi');
        dA = trace(straini)/2;
        dS = sqrt(straini(1,2).^2 + (straini(1,1)-straini(2,2)).^2/4);
        [v d] = eig(straini);
        phi = atan2(v(2),v(1))*180/pi;
        m1 = [m1; dA dS rotati(3) phi straini(1) straini(3) straini(2) straini(4)];
    end
end

shear = reshape(m1(:,2),size(gridy,2),size(gridx,2));


figure
pcolor(gridx/1e3,gridy/1e3,shear)
hold on
shading flat
plot(E/1e3,N/1e3,'k.')
plot(coastE/1e3,coastN/1e3,'k')
plot(faultsE/1e3,faultsN/1e3,'r')
caxis([0 1e-7])
set(gca,'xtick',[],'ytick',[],'position',[0 0 1 1])
axis([100 750 3500 4050])
print('-dpng','shear_higherres_toosmooth.png')








