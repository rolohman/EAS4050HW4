%Note: I like to line up the = signs in matlab, it just makes things look
%prettier.


tmp        = csvread('rates.csv',1); %ignore first row, since it has text column labels
GPS_lon    = tmp(:,1); %Longitude
GPS_lat    = tmp(:,2); %Latitude
GPS_x      = tmp(:,3); %Easting (in meters, UTM zone 11)
GPS_y      = tmp(:,4); %Northing (in meters, UTM zone 11)
GPS_UE     = tmp(:,5); %Displacement rate in the east direction, mm/yr 
GPS_UN     = tmp(:,6); %Displacement rate in the north direction, mm/yr
nd         = length(tmp); %how many points do we have?

%faults and coastlines for plotting.
tmp        = csvread('faults.csv',1);
faults_lon = tmp(:,1);
faults_lat = tmp(:,2);

tmp        = csvread('coasts.csv',1);
coast_lon  = tmp(:,1);
coast_lat  = tmp(:,2);


%the original "map" locations have values like 4e6 meters - this just makes
%the plots a little cleaner.  Doesn't affect any of the math at all.
GPS_x=GPS_x-mean(GPS_x);
GPS_y=GPS_y-mean(GPS_y);

GPS_UE     = GPS_UE/1000;  %convert to meters/yr, since our units of distance are in meters
GPS_UN     = GPS_UN/1000;

%It can be helpful to add a "scale" arrow of known size - so here I
%"augment" the lon, lat and displacement rate vectors by to add in a new
%value sitting somewhere offshore where I will also place a text label.

aug_lon = [GPS_lon;-120]; 
aug_lat = [GPS_lat;32.5];
aug_UE  = [GPS_UE;0.05];
aug_UN  = [GPS_UN;0];

%Below - note that the "quiver" command has a step where it automatically
%scales the arrows so that they "fit" on the screen in a reasonable way.
%If you make the "augmented" value above really big (like, 0.2) you will
%see that all the other ones will appear smaller.

figure('name','THIS IS FIGURE 1')
plot(coast_lon,coast_lat,'k','linewidth',[2])
hold on
plot(faults_lon,faults_lat,'r')
quiver(aug_lon,aug_lat,aug_UE,aug_UN,'color',[0.1 0.5 0.6])
text(-120.1,32.6,'50mm/yr')
axis('image')
axis([-121 -115 32 36])
title('Data')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
%%%Now we invert for the big matrix L.  The positions we use here have to
%%%be in meters (GPS_x and GPS_y), not in latitude/longitude.  But it
%%%doesn't matter that we removed that average value.

%remember this discussion from class?
%G=[ T1 T2 U11 U12 U21 U22]


%We are first doing this for all the data at once. (next week look at spatial variation)
%First, predefine vectors of ones and  of zeros, since we use them a bunch.
z    = zeros(nd,1); 
o    = ones(nd,1);

%The first "nd" rows are the related to the eastward components of
%displacement, the next nd rows are related to the northward component of
%displacement.  You can see that in our data, too.  Note which "chunks" of
%G are zeroed out.

G    = [o z GPS_x GPS_y z z;z o z z GPS_x GPS_y];
data = [GPS_UE;GPS_UN];
Gg   = inv(G'*G)*G';  %This is that master relationship that holds for when you have more data than unknowns.

model   = Gg*data;          %inversion -> model
%now lets divide up the individual parts of the model in symmetric and
%antisymmetric and plot them separately

T = model(1:2); %translation
L = [model(3:4)';model(5:6)'];
W = 1/2*(L-L');
E = 1/2*(L+L');


%below  I'm plotting in the projected coordinates, GPS_x and GPS_y, instead
%of lat/lon.  Could do either way. A benefit of units in meters is that
%they are the same size in both directions, whereas a unit of latitude is
%not usually the same distance as a unit of longitude.
figure     
subplot(2,2,1)
quiver(GPS_x,GPS_y,synthL(1,:)',synthL(2,:)') %I have to take the transpose here, since GPS_x is np x 1 and the synthetics are 1 x np
axis image %this requires that both axes have the same scale, i.e., the same length on the axis is the same "real" distance.
xlabel('East (m)')
ylabel('North (m)')
title('Full matrix L')

subplot(2,2,2)
quiver(GPS_x,GPS_y,synthE(1,:)',synthE(2,:)')
axis image
xlabel('East (m)')
ylabel('North (m)')
title('Just symmetric component, E')

subplot(2,2,3)
quiver(GPS_x,GPS_y,synthW(1,:)',synthW(2,:)')
axis image
xlabel('East (m)')
ylabel('North (m)')
title('Just anti-symmetric component, W')
%Now all together, with better x, y plot scale and a legend
subplot(2,2,4)
%The "1.5" here fixes the scaling of the arrows to a set number - 
%this let's me combine quiver commands ofor different matrices and know that they'll be scaled to the same degree.
quiver((GPS_x-min(GPS_x))/1e3,(GPS_y-min(GPS_y))/1e3,synthL(1,:)',synthL(2,:)',1.5,'k') 
hold on
quiver((GPS_x-min(GPS_x))/1e3,(GPS_y-min(GPS_y))/1e3,synthE(1,:)',synthE(2,:)',1.5,'r')
quiver((GPS_x-min(GPS_x))/1e3,(GPS_y-min(GPS_y))/1e3,synthW(1,:)',synthW(2,:)',1.5,'b')
axis image
xlabel('East (km)')
ylabel('North (km)')
legend('L','E','W')
title('All components')


%now the eigenvector stuff
[v,e]  = eig(E); %matrix of eigenvectors (u) and eigenvalues (diag e)

%plot principle components with arrows

figure
quiver(GPS_x/1e3,GPS_y/1e3,GPS_UE,GPS_UN,1.5,'color',[0.8 0.8 0.8])
hold on, axis image
xlabel('East (km)')
ylabel('North (km)')
%plot arrows pointing in or out, depending on the sign of that eigenvalue
for i=1:2
    if(e(i,i)>0)
        quiver(0,0,v(1,i)*e(i,i),v(2,i)*e(i,i),1e9,'r','linewidth',5,'maxheadsize',20)
        quiver(0,0,-v(1,i)*e(i,i),-v(2,i)*e(i,i),1e9,'r','linewidth',5,'maxheadsize',20)
    else
        quiver(v(1,i)*e(i,i)*1e9,v(2,i)*e(i,i)*1e9,-v(1,i)*e(i,i),-v(2,i)*e(i,i),1e9,'b','linewidth',5,'maxheadsize',20)
        quiver(-v(1,i)*e(i,i)*1e9,-v(2,i)*e(i,i)*1e9,v(1,i)*e(i,i),v(2,i)*e(i,i),1e9,'b','linewidth',5,'maxheadsize',20)
    end
end

