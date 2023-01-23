function [interferogram] = simulateOCTinterferogram(wavelength, intensityWavlength, numberReflectors, depth)

    %% Linearize spectrum with respect to wavenumber
    wavenumber = 2*pi./wavelength;
    wavnumberLinear = linspace(max(wavenumber),min(wavenumber),numel(wavelength));
    deltaK = wavnumberLinear(1) - wavnumberLinear(2);

    % Resample the spectrum
    intensityWavenumber = interp1(wavenumber, intensityWavlength, wavnumberLinear);

    %% Simulate the interferogram 
    
    numLayer = size(numberReflectors);
    reflector_signal = 0;

    % Initialize the frequency response function
    H = 0;
    
    for layer = 1:numLayer(2)-1

        % Calculate the refelection coefficient 
        reflector_n = numberReflectors(layer+1)-(numberReflectors(layer))/...
                    (numberReflectors(layer+1)+numberReflectors(layer));


        reflector_signal = reflector_signal + numberReflectors(layer)...
                    * (depth(layer+1) - depth(layer));

        H = H + reflector_n*exp(1i*2.*wavnumberLinear*reflector_signal);

    end
    
    % Calculate the cross-crorelation term
    interferogram = zeros(1,numel(wavelength));
    
    
% min normalized interferogram 
realH = real(H) - min(real(H));
    for layer=1:length(wavnumberLinear)

        interferogram(layer) = intensityWavenumber(layer)...
                    *realH(layer);

    end

end

