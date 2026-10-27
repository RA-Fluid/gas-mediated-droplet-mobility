function Corrs = RA_Droplet_Correlations0(~, rho_c, mu_c, rho_d, mu_d, D_d)

%% main inputs

% gravitational acceleration
g = 9.81; % [m/s^2]

%% density difference & viscosity ratio

delta_rho = abs(rho_c - rho_d) ; % [kg/m3]
mu_s = mu_d / mu_c ;

%% constant of the equation

c = ( (4/3) * delta_rho * rho_c * g * (mu_c^-2) ) * (D_d.^3) ; % []
options = optimoptions('fsolve', 'Display', 'off');

%% Schiller & Naumann: rigid sphere only

Re_0 = D_d*100 ;
Re = zeros(size(c)) ;

for i = 1 : numel(c)
    fun = @(Re)RA_Schiller_Naumann(Re, c(i));
    Re(i) = fsolve(fun, Re_0(i), options) ;
end

c_d = (24./Re) .* (1+0.15.*(Re.^0.687)) ;
U_d = Re .* mu_c ./ (D_d .* rho_c) ; % [m/s]
U_d = U_d .* 1e3 ; % [mm/s]

% remove nonsense
ins = find(abs (U_d-250) < 250) ;

% Schiller-Numann is valid only for Re < 800
iv = find(Re < 800) ;

ic = intersect(ins, iv) ;

Schiller_Numann.D_d = D_d(ic)*1e3 ;
Schiller_Numann.Re = Re(ic) ;
Schiller_Numann.c_d = c_d(ic) ;
Schiller_Numann.U_d = U_d(ic) ;

%% Rybszynski: viscous sphere

Re_0 = 0.1 ;
Re = zeros(size(c)) ;

for i = 1 : numel(c)
    fun = @(Re)RA_Rybszynski(Re, c(i), mu_s);
    Re(i) = fsolve(fun, Re_0, options) ;
end

c_d = (24./Re) .* ( ( (2/3) + mu_s ) / (1 + mu_s ) ) ;
U_d = Re .* mu_c ./ (D_d .* rho_c) ; % [m/s]
U_d = U_d .* 1e3 ; % [mm/s]

% Rybszynski is valid only for Re < 1
iv = find(Re < 1) ;

Rybszynski.D_d = D_d(iv)*1e3 ;
Rybszynski.Re = Re(iv) ;
Rybszynski.c_d = c_d(iv) ;
Rybszynski.U_d = U_d(iv) ;

%% Feng and Michaelides: viscous sphere

Re_0 = 0.1 ;
Re = zeros(size(c)) ;

for i = 1 : numel(c)
    fun = @(Re)RA_Feng_Michaelides(Re, c(i), mu_s);
    Re(i) = fsolve(fun, Re_0, options) ;
end

c_d_0 = (48./Re) .* ( 1 + 2.21*(Re.^-0.5) - 2.14./Re ) ;
c_d_2 = 17 .* (Re.^(-2/3)) ;
c_d_inf = (24./Re) .* (1 + Re.^(2/3) / 6) ;

if mu_s <= 2

    c_d = 0.5*(2-mu_s).*c_d_0 + ( (4*mu_s)/(6+mu_s) ) .* c_d_2 ;

else
    c_d = (4/(mu_s+2)) .* c_d_2 + ((mu_s-2)/(mu_s+2)).*c_d_inf ;
end

U_d = Re .* mu_c ./ (D_d .* rho_c) ; % [m/s]
U_d = U_d .* 1e3 ; % [mm/s]

% Feng and Michaelides is valid only for 5 <= Re <= 1000
iv = find( abs( Re - 1005/2 ) <= 1000 - (1005/2) ) ;

Feng_Michaelides.D_d = D_d(iv)*1e3 ;
Feng_Michaelides.Re = Re(iv) ;
Feng_Michaelides.c_d = c_d(iv) ;
Feng_Michaelides.U_d = U_d(iv) ;

%% Feng and Michaelides low Re: viscous sphere

Re_0 = 0.1 ;
Re = zeros(size(c)) ;

for i = 1 : numel(c)
    fun = @(Re)RA_Feng_Michaelides_LowRe(Re, c(i), mu_s);
    Re(i) = fsolve(fun, Re_0, options) ;
end

a = (3*mu_s+2)/(1+mu_s) ;
c_d = (8./Re) * a .* (1 + 0.05 * a .* Re ) - 0.01*a.*Re.*log(Re) ;

U_d = Re .* mu_c ./ (D_d .* rho_c) ; % [m/s]
U_d = U_d .* 1e3 ; % [mm/s]

% is valid only for R <= 5
iv = find( Re <= 5 ) ;

Feng_Michaelides_lowRe.D_d = D_d(iv)*1e3 ;
Feng_Michaelides_lowRe.Re = Re(iv) ;
Feng_Michaelides_lowRe.c_d = c_d(iv) ;
Feng_Michaelides_lowRe.U_d = U_d(iv) ;

%% Feng and Michaelides The most general: viscous sphere + slip

beta = [0, 0.1, 1, 10, 100, 10000] ;

Re_0 = 0.01 ;
Re = zeros(length(c), length(beta)) ;

for i = 1 : numel(c)
    for j = 1 : numel(beta)

        fun = @(Re_s)RA_Feng_Michaelides_Full(Re_s, c(i), mu_s, beta(j));
        Re(i,j) = fsolve(fun, Re_0, options) ;

    end
end

a = (2 + 3*mu_s + 6*beta*mu_s)./(1+mu_s+3*beta*mu_s) ;

c_d = (8./Re) .* (a + Re.*(a.^2/16) + (a.^3/160).*( Re.^2 .* log(Re/2) ) ) ;

U_d = Re .* mu_c ./ (D_d .* rho_c) ; % [m/s]
U_d = U_d .* 1e3 ; % [mm/s]

%is valid only for Re <= 2
Feng_Michaelides_Full = struct([]) ;
for j = 1 : numel(beta)
    iv = find( Re(:,j) <= 2 ) ;

    Feng_Michaelides_Full(j).beta = beta(j) ;
    Feng_Michaelides_Full(j).D_d = D_d(iv)*1e3 ;

    Feng_Michaelides_Full(j).Re = Re(iv,j) ;
    Feng_Michaelides_Full(j).c_d = c_d(iv,j) ;
    Feng_Michaelides_Full(j).U_d = U_d(iv,j) ;
end

%% save outputs

Corrs.SolidSphere.Schiller_Numann = Schiller_Numann ;

Corrs.ViscousSphere.Rybszynski = Rybszynski ;
Corrs.ViscousSphere.Feng_Michaelides = Feng_Michaelides ;
Corrs.ViscousSphere.Feng_Michaelides_lowRe = Feng_Michaelides_lowRe ;
Corrs.ViscousSphere.Feng_Michaelides_Full = Feng_Michaelides_Full ;

end

function f = RA_Schiller_Naumann(Re, c)
% Schiller-Naumann correlation for RIGID SPHERES only
% vaid only for Re < 800

c_d = (24/Re) * (1+0.15*(Re^0.687)) ;
f = c_d - c / Re^2 ;

end

function f = RA_Rybszynski(Re, c, mu_s)
% Rybszynski correlation for VISCOUS SPHERES
% vaid only for Re << 1

c_d = (24/Re) * ( ( (2/3) + mu_s ) / (1+mu_s ) ) ;
f = Re^2 * c_d - c ;

end

function f = RA_Feng_Michaelides(Re, c, mu_s)
% Feng and Michaelides correlation for VISCOUS SPHERES
% vaid only for 5 <= Re <= 1000

c_d_0 = (48/Re) * ( 1 + 2.21*(Re^-0.5) - 2.14/Re ) ;
c_d_2 = 17 * (Re^(-2/3)) ;
c_d_inf = (24/Re) * (1 + Re^(2/3) / 6) ;

if mu_s <= 2

    c_d = 0.5*(2-mu_s)*c_d_0 + ( (4*mu_s)/(6+mu_s) ) * c_d_2 ;

else
    c_d = (4/(mu_s+2)) * c_d_2 + ((mu_s-2)/(mu_s+2))*c_d_inf ;
end

f = Re^2 * c_d - c ;

end

function f = RA_Feng_Michaelides_LowRe(Re, c, mu_s)
% Suggested correlation by Feng and Michaelides for VISCOUS SPHERES
% when R <= 5

a = (3*mu_s+2)/(1+mu_s) ;

c_d = (8/Re) * a * (1 + 0.05 * a * Re ) - 0.01*a*Re*log(Re) ;

f = Re^2 * c_d - c ;

end

function f = RA_Feng_Michaelides_Full(Re, c, mu_s, beta)
% The most general correlation by Feng and Michaelides for VISCOUS SPHERES
% with slip coefficient beta (beta = 0 no slip and b = inf perfect slip)
% for small but finite Reynolds numbers
% IN THE ORIGINAL EQUATION, RE IS DEFINED BASED ON RADIUS; SO I MODIFIED
% THE EQUATION FOR RE DEFINED BASED ON DIAMETER

a = (2 + 3*mu_s + 6*beta*mu_s)/(1+mu_s+3*beta*mu_s) ;

c_d = (8/Re) * (a + (a^2/16)*Re + (a^3/160)*( Re^2 * log(Re/2) ) ) ;

f = Re^2 * c_d - c ;

end