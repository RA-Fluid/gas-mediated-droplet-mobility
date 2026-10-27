function mu = RA_DynViscosity_Correction_0(T, name)

switch name
    case 'Air'
        mu = ( 1.729 + 0.0049 * T ) * 1e-2 ; % [mPa.s]
    case 'Water'
        % mu = XSteam('my_PT', 1.01325, T)*1e3 ;
        mu = 1.3*exp(-0.0395*T)+0.4777*exp(-0.0061*T) ; % [mPas] % from Engineering toolbox; fitted the line myself using cftool
    case 'Hexadecane'
        mu = 0.6311 * exp( 239.8917 / (T + 106.2064) ) ; % from viscometry measurements and fitting the Vogel–Fulcher–Tammann (VFT) model
    case 'SO-10cst'
        mu = 0.1057 * exp( 969 / (T + 189.143) ) ; % from viscometry measurements and fitting the Vogel–Fulcher–Tammann (VFT) model        
    case 'SO-20cst'
        mu = 0.1652 * exp( 1261 / (T + 229.1) ) ; % from viscometry measurements and fitting the Vogel–Fulcher–Tammann (VFT) model
end