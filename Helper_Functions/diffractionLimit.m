%% CHECKKKK

function  diffractionLimit(lambda)
    
    sp = systemParameters();
    
    % diffraction limit:
    d = lambda / (2*sp.NA);
    disp("d = %.2f", d);

    D = 2*sp.f*sp.NA;
end