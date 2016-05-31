clear all;

scrsz = get(0,'ScreenSize');

scrsz(3) = 1920;

[filename, pathname] = uigetfile( '*.bak','Pick a file');

savefig = 0;

read_bak(pathname,filename)

filename = strcat(pathname,filename);

filename = strrep(filename,'.bak','');

run(strcat(filename,'.m'));

if ( strcmp(DISTRIBUTION,'lognormal'))
    
    [ ~,mom ] = lognormal_moments( mu,sigma,gas_volume_fraction,...
        solid_mass_fractions,d1,d2,rho1,rho2);
    
    % return
    
    
    if savefig
        
        export_fig(strcat(filename,'_sizes.pdf'));
        
    end
    
end

importfile(strcat(filename,'.col'))
prova = data;
%  pause
% 
% prova = importdata(strcat(filename,'.col'),' ',1);

moments = importdata(strcat(filename,'.mom'));

n_part = moments(1);
n_mom = moments(2);

moments =  (reshape(moments(3:end)',n_part*(n_mom+1)+1,[]))';

z_levels = size(moments,1);

prova = reshape(prova,z_levels,[]);

if ( n_mom > 1 ) 
    
    moments = reshape(moments(:,2:end),[],n_part,n_mom+1);
    
    set_mom = moments(:,:,end);
    moments = moments(:,:,1:end-1);
    
else
        
    moments = moments(:,2:n_mom);
    
end

z = prova(:,1)/1000;
r_mt = prova(:,2);
r = prova(:,2)/1000;
x = prova(:,3)/1000;
y = prova(:,4)/1000;
rho_mix = prova(:,5);
temp = prova(:,6);
w = prova(:,7);
mag_u = prova(:,8);
atm_mass_fraction = prova(:,9);
wvapour_mass_fraction = prova(:,10);

gas_mass_fraction = atm_mass_fraction + wvapour_mass_fraction;
for i=1:n_part,

    solid_partial_mass_fraction(:,i) = prova(:,10+i);

end


rho_atm = prova(:,10+n_part+1);
mfr = prova(:,10+n_part+2);
% ta =  prova(:,9+n_part+3);

n_z = length(z);

   
rho_rel = rho_mix - rho_atm;

%%

%figure('Position',[1 scrsz(4)/2 2*scrsz(3)/3 2*scrsz(4)/3])
figure
subplot(1,3,1)

plot(atm_mass_fraction,z, wvapour_mass_fraction,z,gas_mass_fraction,z);

xlabel('Gas mass fraction');
ylabel('Height (km)');
box on;

%%


subplot(1,3,2)

hold all;



for i=1:n_part,

    solid_mass_fraction(:,i) = solid_partial_mass_fraction(:,i) ...
        .* ( 1.D0 - gas_mass_fraction(:) );

%    plot(solid_partial_mass_fraction(:,i),z,'--');
    plot(solid_mass_fraction(:,i),z,'-');
    
           
end

xlabel('Particles mass fractions');
ylabel('Height (km)');
box on;

%%


subplot(1,3,3)

hold all;

solid_mass_loss(1,1:n_part) = 0.D0
for i=1:n_part,

    solid_mass_flux(:,i) = solid_mass_fraction(:,i) .* rho_mix(:) * pi ...
        .* r_mt(:).^2 .* mag_u(:);

       
    solid_mass_loss(2:z_levels,i) = -diff(solid_mass_flux(:,i)./ ...
        solid_mass_flux(1,i) );

%    plot(solid_mass_loss(:,i),z,'-');

    solid_mass_loss_cum(:,i) = ( 1.D0 - solid_mass_flux(:,i) ./ ...
        solid_mass_flux(1,i) );

    plot(solid_mass_loss_cum(:,i),z,'-');

end


xlabel('Particles mass lost fractions');
ylabel('Height (km)');
box on;

if savefig

    export_fig(strcat(filename,'_mom.pdf'));

end

%%

figure

if ( n_part == 1 )
    
    subplot(1,4,1)
    
%    plot(-log2(1000*moments(:,4)./moments(:,3)),z);
    plot(moments(:,2)./moments(:,1),z);
    
    xlabel({'Volume avg. mean';'diameter (\phi)'});
    ylabel('z (m)');
    ylabel('Height (km)');
    box on;
    
    
else
    
    for i=1:n_part,
        
        subplot(2,n_part,i)
        
        % plot(-log2(1000*moments(:,i,4)./moments(:,i,3)),z);
        plot(moments(:,i,2)./moments(:,i,1),z);
        
        xlabel({'\mu (\phi)'});
        box on;
               
    end
    
    ylabel('Height (km)');

end



%%


if ( n_part == 1 )
    
    subplot(1,4,2)
    
    k1 = log( 1.D3 * moments(:,2) ./ moments(:,1) );
    k2 = log( 1.D3 * moments(:,4) ./ moments(:,3) );
    
    mu_test = - 0.25D0 * ( 5*k2(:) - k1(:) ) / log(2);
    sigma_test = sqrt( 0.5d0 * ( k2(:)-k1(:) ) ) / log(2);
    
    sigma_phi = sqrt( -(moments(:,2)./moments(:,1)).^2 + (moments(:,3)./moments(:,1)) );
    plot(sigma_phi,z);
    
%    plot(mu_test,z,'-',mu_test+sigma_test,z,'--',mu_test-sigma_test,z,'--')
    
    xlabel('\sigma (\phi)');
    ylabel('Height (km)');
    box on;

    subplot(1,4,3)

    skew_phi = 2.D0 * (moments(:,2)./moments(:,1)).^3 ...
            - 3.D0 * (moments(:,2)./moments(:,1)) .* ...
            (moments(:,3) ./ moments(:,1)) + (moments(:,4)./moments(:,1));
 
    plot(skew_phi,z);

    xlabel('Skew (\phi)');
    ylabel('Height (km)');
    box on;


    subplot(1,4,4)
    plot(100*solid_mass_loss_cum(:,i),z,'-');
    xlabel('Solid Mass Flux lost (%)');
    ylabel('Height (km)');
    box on;
   
else
    
    
    for i=1:n_part,
        
        subplot(2,n_part,n_part+i)
        
        k1(:,i) = log( 1.D3 * moments(:,i,2) ./ moments(:,i,1) );
        k2(:,i) = log( 1.D3 * moments(:,i,4) ./ moments(:,i,3) );
        
        mu_test(:,i) = - 0.25D0 * ( 5*k2(:,i) - k1(:,i) ) / log(2);
        sigma_test(:,i) = sqrt( 0.5d0 * ( max(0,k2(:,i)-k1(:,i) )) ) / log(2);
        
        plot(mu_test(:,i),z,'-',mu_test(:,i)+sigma_test(:,i),z,'--',...
            mu_test(:,i)-sigma_test(:,i),z,'--')
        
        xlabel('phi');
        box on;
        
    end
    
    ylabel('Height (km)');

end




%%

%figure('Position',[1 scrsz(4)/2 2*scrsz(3)/3 2*scrsz(4)/3])
figure
subplot(2,2,1)

% plot(x,z,x+r,z,x-r,z);
plot(x+r,z);
xlabel('Radius (km)');
ylabel('Height (km)');

%xlim([ 0 r(ceil(0.999*size(r,1))) ]);
hold all;

subplot(2,2,2)

plot(w,z);
xlabel('Velocity (m/s)');
ylabel('Height (km)');
hold all;


subplot(2,2,3)
% plot(temp,z,ta,z);
% xlabel('Temperature (K)');
% ylabel('Height (km)');
% hold all;

plot(rho_mix,z);
xlabel('Mixture density (kg/m^3)');
ylabel('Height (km)');
hold all;


subplot(2,2,4)

plot(rho_rel,z);
xlabel('Relative density (kg/m^3)');
ylabel('Height (km)');
hold all;



if savefig

    export_fig(strcat(filename,'_profiles.pdf'));

end

%%

%figure('Position',[1 scrsz(4)/2 2*scrsz(3)/3 2*scrsz(4)/3])
figure
plot3(x,y,z);

hold all;

box on;

angle = linspace(0,2*pi,50);

x_plume = cos(angle);
y_plume = sin(angle);

z_max = max(z);
z_min = min(z);

n_sect = 50;

zeta_grid = linspace(z_min,z_max*0.99,n_sect);

for i=2:length(x),
    
   l_seg(i) = sqrt( (x(i)-x(i-1))^2 + (y(i)-y(i-1))^2 + (z(i)-z(i-1))^2 ); 
    
end

s_axis(1) = 0.0;
s_axis(2:length(x)) = cumsum( l_seg(2:length(x)) );

s_grid = linspace(0,max(s_axis)*0.99,n_sect);

for i=1:n_sect,
    
   ind = find(z-zeta_grid(i)>0, 1,'first');
   ind = find(s_axis-s_grid(i)>0, 1,'first');
  
   vect(1) = x(ind) - x(ind-1);
   vect(2) = y(ind) - y(ind-1);
   vect(3) = z(ind) - z(ind-1);
    
   vect = vect / norm(vect,2);

   vect0(1) = 0;
   vect0(2) = 0;
   vect0(3) = 1;
   
   v = cross(vect0,vect);
   
   s = norm(v,2);
   
   c = dot(vect0,vect);
   
   mat_v = zeros(3,3);
   mat_v(2,1) = v(3);
   mat_v(1,2) = -v(3);
   
   mat_v(3,1) = -v(2);
   mat_v(1,3) = v(2);
   
   mat_v(2,3) = -v(1);
   mat_v(3,2) = v(1);

   R = eye(3) + mat_v + mat_v^2 * ( 1-c ) / s^2;
   
   plume = [ r(ind)*x_plume ;r(ind)*y_plume; zeros(1,50) ];
   plume_rotated = R*plume;
   
   plot3(x(ind)+plume_rotated(1,:),y(ind)+plume_rotated(2,:),...
       z(ind)+plume_rotated(3,:));
   
end

axis equal;

xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');


if savefig

    export_fig(strcat(filename,'_plume.pdf'));

end

%%
grid_size = 40;

xgrid = linspace(min(xlim),max(xlim),grid_size);
ygrid = linspace(min(ylim),max(ylim),grid_size);
zgrid = linspace(min(zlim),max(zlim),grid_size);

[X,Y,Z] = meshgrid(xgrid,ygrid,zgrid);

C = zeros(grid_size,grid_size,grid_size);
I = zeros(grid_size,grid_size,grid_size);

for i=1:length(xgrid),
    
    for j=1:length(ygrid),

        for k=1:length(zgrid),

            [C(i,j,k),I(i,j,k)] = min( (x-xgrid(i)).^2 + (y-ygrid(j)).^2 + (z-zgrid(k)).^2 );

        end
        
    end
    
end
  
C = sqrt(C);

B = solid_mass_fraction(I) .*exp( -C.^2 ./ r(I) );

% vtkwrite(strcat(filename,'.vtk'),'structured_grid',X,Y,Z,'scalars','mass_fraction',B)

% figure('Position',[1 scrsz(4)/2 2*scrsz(3)/3 2*scrsz(4)/3])
% 
% hold all;
% 
% for i=1:n_part,
% 
%     plot(rho_mix.*r.^2*pi.*solid_partial_mass_fraction(:,i),z);
% 
% end
% 
% 
% 
% xlabel('Particles mass per unit length (m)');
% ylabel('Height (km)');
% box on;



%%

% mom0 = zeros(n_z,n_part,n_mom);
% 
% for i_mom=1:n_mom-1,
% 
%     mom0(1:n_z,1:n_part,i_mom) = moments(1:n_z,1:n_part,i_mom+1) ./ ...
%         moments(1:n_z,1:n_part,1);
% 
% end
% 
% mu(1:n_z,1:n_part,1) = 0;
% 
% mu(1:n_z,1:n_part,2) = - mom0(1:n_z,1:n_part,1).^2 + mom0(1:n_z,1:n_part,2);
% 
% mu(1:n_z,1:n_part,3) = 2.D0 * mom0(1:n_z,1:n_part,1).^3 ...
%     - 3.D0 * mom0(1:n_z,1:n_part,1) .* mom0(1:n_z,1:n_part,2) ...
%     + mom0(1:n_z,1:n_part,3);
% 
% mu(1:n_z,1:n_part,4) = - 3.D0 * mom0(1:n_z,1:n_part,1).^4 ...
%             + 6.D0 * mom0(1:n_z,1:n_part,1).^2 .* mom0(1:n_z,1:n_part,2) ...
%             - 4.D0 * mom0(1:n_z,1:n_part,1) .* mom0(1:n_z,1:n_part,3) ...
%             + mom0(1:n_z,1:n_part,4);
% 
% mean = mu(1:n_z,1:n_part,1);        
%         
% variance = mu(1:n_z,1:n_part,2);
% 
% skewness = mu(1:n_z,1:n_part,3) ./ mu(1:n_z,1:n_part,2).^(3.d0/2.d0);
% 
% kurtosis = mu(1:n_z,1:n_part,4) ./ ( mu(1:n_z,1:n_part,2).^2 );
% 
% i_z = 1;
% i_part = 1;
% 
% mom = {mean(i_z,i_part),variance(i_z,i_part),skewness(i_z,i_part),...
%     kurtosis(i_z,i_part)};
% 
% figure
% 
% [r,type] = pearsrnd(mom{:},10000,1);
% 
% [Fi,xi] = ecdf(r);
% hold all;
% stairs(xi,Fi);
% 
% i_z = ceil(n_z/2);
% 
% mom = {mean(i_z,i_part),variance(i_z,i_part),skewness(i_z,i_part),...
%     kurtosis(i_z,i_part)};
% 
% [r,type] = pearsrnd(mom{:},10000,1);
% 
% [Fi,xi] = ecdf(r);
% stairs(xi,Fi);

z_norm = (z-min(z))/(max(z)-min(z));
w_norm = w / max(w);

first_der = ( z_norm(2:end) - z_norm(1:end-1) ) ./ ( w_norm(2:end) - w_norm(1:end-1) );
sec_der = ( first_der(2:end) - first_der(1:end-1) ) ...
    ./ ( 0.5 * ( w_norm(3:end) - w_norm(1:end-2 ) ) );

first_der_central = ( z_norm(3:end) - z_norm(1:end-2) ) ./ ( w_norm(3:end) - w_norm(1:end-2) );

k = sec_der ./ ( ( 1.0 + first_der_central.^2 ).^(1.5d0) );

max_k = max(k)

Radius = 1 / max_k