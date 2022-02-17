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
faults_x =   tmp(:,1); %
faults_y =   tmp(:,2);
tmp        = csvread('coasts.csv',1);
coast_lon  = tmp(:,1);
coast_lat  = tmp(:,2);


%the original "map" locations have values like 4e6 meters - this just makes
%the plots a little cleaner.  Doesn't affect any of the math at all.
GPS_x=GPS_x-mean(GPS_x);
GPS_y=GPS_y-mean(GPS_y);

GPS_UE     = GPS_UE/1000;  %convert to meters/yr, since our units of distance are in meters
GPS_UN     = GPS_UN/1000;

%%%Now we are going to make a grid of points, going from a range of x and y
%%%values that contains most of the region.  You can change these numbers
%%%if you'd like!
dx=10e3; %grid spacing in meters (could have a dy and dx but we are making them the same)
gridy = min(GPS_y):dx:max(GPS_y); %this makes a vector of y values, spaced by dx
gridx = min(GPS_x):dx:max(GPS_x);

%how many x and y points do we have?  We use that later.
nx=size(gridx,2); 
ny=size(gridy,2);

%This makes a grid of x and y values from our two vectors
[x,y] = meshgrid(gridx,gridy);

%smoothing distance - this allows us to use points that are not just within
%our single grid box, while still making the ones in the box more
%"important".  Change this up!
alpha = 50e3; %in meters

m = zeros(6,nx*ny); %initialize an empty model vector 
count = zeros(size(x)); %keep a vector of how many points we end up using for each box

for i=1:nx*ny %loop over all the nx*np boxes
    
    %Calculate the distance from ALL GPS points to your grid box center.
    dists   = sqrt((GPS_x-x(i)).^2 +(GPS_y-y(i)).^2);  % the .^2 symbol makes matlab square each element of the vector separately
    weights = exp(-(dists.^2)/(2*alpha^2));            % this is the equation for a Gaussian curve, with width alpha
    ind     = find(or(dists<dx,dists < 2*alpha));      % only use points that are within 2*alpha from the center of our box, or the box size itself
    nd      = size(ind,1);                             % number of "good" points
    
    if (nd < 6)   %too few nearby points left!
        m(:,i)   = [NaN(6,1)]; %just set this to not-a-number
    else
        %make a diagonal matrix of weights for the points we are using,
        %repeated since we have E and N components
        sig      =  diag([weights(ind); weights(ind)]);
        
        %data vector of just the points we are using
        d        = [GPS_UE(ind); GPS_UN(ind)];
        z        = zeros(nd,1); 
        o        = ones(nd,1);
        %build G matrix
        G        = [o z GPS_x(ind) GPS_y(ind) z z; z o z z GPS_x(ind) GPS_y(ind)];
        Gg       = inv(G'*sig*G)*G'*sig;
        m(:,i)   = Gg*d;
        count(i) = nd;
    end
end

%now we want to pull out the extension and shear components
volumetric = nan(size(x));
shear      = nan(size(x));
angle      = nan(size(x));

for i=1:nx*ny
    L  = [m(3,i) m(4,i); m(5,i) m(6,i)];
    if(~isnan(L))
        symm = 0.5*(L + L');
        asym  = 0.5*(L - L');
        
        volumetric(i) = trace(symm);
        shear(i)      = sqrt(symm(1,2)^2+(symm(1,1)-symm(2,2))^2/4);
        [v,d]         = eig(symm);
        angle(i)      = atan2(v(2),v(1))*180/pi;
        roti(i)       = asym(2);
        
    end
end
    %         dA = trace(straini)/2;
%         dS = sqrt(straini(1,2).^2 + (straini(1,1)-straini(2,2)).^2/4);
%         [v d] = eig(straini);
%         phi = atan2(v(2),v(1))*180/pi;
%         m = [m; dA dS rotati(3) phi straini(1) straini(3) straini(2) straini(4)];
%     end
% end
% 
% ext1 = reshape(m(:,1),size(gridy,2),size(gridx,2));
% shear = reshape(m(:,2),size(gridy,2),size(gridx,2));
% rot = reshape(m(:,3),size(gridy,2),size(gridx,2));
% ang = reshape(m(:,4),size(gridy,2),size(gridx,2));
% exx = reshape(m(:,5),size(gridy,2),size(gridx,2));
% eyy = reshape(m(:,8),size(gridy,2),size(gridx,2));
% exy = reshape(m(:,6),size(gridy,2),size(gridx,2));
% 
% 
% %We are first doing this for all the data at once. (next week look at spatial variation)
% %First, predefine vectors of ones and  of zeros, since we use them a bunch.
% z    = zeros(nd,1); 
% o    = ones(nd,1);
% 
% %The first "nd" rows are the related to the eastward components of
% %displacement, the next nd rows are related to the northward component of
% %displacement.  You can see that in our data, too.  Note which "chunks" of
% %G are zeroed out.
% 
% G    = [o z GPS_x GPS_y z z;z o z z GPS_x GPS_y];
% data = [GPS_UE;GPS_UN];
% Gg   = inv(G'*G)*G';  %This is that master relationship that holds for when you have more data than unknowns.
% 
% model   = Gg*data;          %inversion -> model
% %now lets divide up the individual parts of the model in symmetric and
% %antisymmetric and save them separately
% 
% T = model(1:2); %translation
% L = [model(3:4)';model(5:6)'];
% W = 1/2*(L-L');
% E = 1/2*(L+L');
% 
% 
% %below  I'm plotting in the projected coordinates, GPS_x and GPS_y, instead
% %of lat/lon.  Could do either way. A benefit of units in meters is that
% %they are the same size in both directions, whereas a unit of latitude is
% %not usually the same distance as a unit of longitude.
% figure     
% subplot(2,2,1)
% quiver(GPS_x,GPS_y,synthL(1,:)',synthL(2,:)') %I have to take the transpose here, since GPS_x is np x 1 and the synthetics are 1 x np
% axis image %this requires that both axes have the same scale, i.e., the same length on the axis is the same "real" distance.
% xlabel('East (m)')
% ylabel('North (m)')
% title('Full matrix L')
% 
% subplot(2,2,2)
% quiver(GPS_x,GPS_y,synthE(1,:)',synthE(2,:)')
% axis image
% xlabel('East (m)')
% ylabel('North (m)')
% title('Just symmetric component, E')
% 
% subplot(2,2,3)
% quiver(GPS_x,GPS_y,synthW(1,:)',synthW(2,:)')
% axis image
% xlabel('East (m)')
% ylabel('North (m)')
% title('Just anti-symmetric component, W')
% %Now all together, with better x, y plot scale and a legend
% subplot(2,2,4)
% %The "1.5" here fixes the scaling of the arrows to a set number - 
% %this let's me combine quiver commands ofor different matrices and know that they'll be scaled to the same degree.
% quiver((GPS_x-min(GPS_x))/1e3,(GPS_y-min(GPS_y))/1e3,synthL(1,:)',synthL(2,:)',1.5,'k') 
% hold on
% quiver((GPS_x-min(GPS_x))/1e3,(GPS_y-min(GPS_y))/1e3,synthE(1,:)',synthE(2,:)',1.5,'r')
% quiver((GPS_x-min(GPS_x))/1e3,(GPS_y-min(GPS_y))/1e3,synthW(1,:)',synthW(2,:)',1.5,'b')
% axis image
% xlabel('East (km)')
% ylabel('North (km)')
% legend('L','E','W')
% title('All components')
% 
% 
% %now the eigenvector stuff
% [v,e]  = eig(E); %matrix of eigenvectors (u) and eigenvalues (diag e)
% 
% %plot principle components with arrows
% 
% figure
% quiver(GPS_x/1e3,GPS_y/1e3,GPS_UE,GPS_UN,1.5,'color',[0.8 0.8 0.8])
% hold on, axis image
% xlabel('East (km)')
% ylabel('North (km)')
% %plot arrows pointing in or out, depending on the sign of that eigenvalue
% for i=1:2
%     if(e(i,i)>0)
%         quiver(0,0,v(1,i)*e(i,i),v(2,i)*e(i,i),1e9,'r','linewidth',5,'maxheadsize',20)
%         quiver(0,0,-v(1,i)*e(i,i),-v(2,i)*e(i,i),1e9,'r','linewidth',5,'maxheadsize',20)
%     else
%         quiver(v(1,i)*e(i,i)*1e9,v(2,i)*e(i,i)*1e9,-v(1,i)*e(i,i),-v(2,i)*e(i,i),1e9,'b','linewidth',5,'maxheadsize',20)
%         quiver(-v(1,i)*e(i,i)*1e9,-v(2,i)*e(i,i)*1e9,v(1,i)*e(i,i),v(2,i)*e(i,i),1e9,'b','linewidth',5,'maxheadsize',20)
%     end
% end
% 
