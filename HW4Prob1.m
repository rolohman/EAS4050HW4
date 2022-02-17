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
faults_x =   tmp(:,3); 
faults_y =   tmp(:,4);

tmp        = csvread('coasts.csv',1);
coast_lon  = tmp(:,1);
coast_lat  = tmp(:,2);

%the original "map" locations have values like 4e6 meters - this just makes
%the plots a little cleaner.  Doesn't affect any of the math at all.
faults_x  = faults_x-mean(GPS_x);
faults_y  = faults_y-mean(GPS_y);
GPS_x     = GPS_x-mean(GPS_x);
GPS_y     = GPS_y-mean(GPS_y);

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
    
    %make a diagonal matrix of weights for the points we are using,
    %repeated since we have E and N components
    sig      =  diag([weights(ind); weights(ind)]);
    
    %data vector of just the points we are using
    d        = [GPS_UE(ind); GPS_UN(ind)];
    z        = zeros(nd,1);
    o        = ones(nd,1);
    %build G matrix
    G        = [o z GPS_x(ind) GPS_y(ind) z z; z o z z GPS_x(ind) GPS_y(ind)];
   
    if(rcond(G'*sig*G)<eps) %this means the matrix is not invertible, probably because you have < 6 points or points that are very close together
        m(:,i)=NaN;
    else
        Gg       = inv(G'*sig*G)*G'*sig;
        m(:,i)   = Gg*d;
        count(i) = nd;
    end
end

%now we want to pull out the extension and shear components
volumetric = nan(size(x));
shear      = nan(size(x));

for i=1:nx*ny
    L  = [m(3,i) m(4,i); m(5,i) m(6,i)];
    if(~isnan(L))
        symm = 0.5*(L + L');
        asym  = 0.5*(L - L');
        
        volumetric(i) = trace(symm);
        shear(i)      = sqrt(symm(1,2)^2+(symm(1,1)-symm(2,2))^2/4);    
    end
end
  
figure
pcolor(x,y,volumetric)
hold on
shading flat
plot(faults_x,faults_y,'r')
plot(GPS_x,GPS_y,'k.')
caxis([-1e-7 1e-7])