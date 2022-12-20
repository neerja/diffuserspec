% createSourceSpectrum.m
% Neerja Aggarwal
% June 8th, 2022
% Purpose: import the spf2 spectrum. find peak (using either max function
% or findpeak function).  store the wavelength and intensity. 

cd transferMatrix_gtWavelengths

files = dir('*.spf2')

N = length(files)

wvcal = zeros(N,1)
wvint = zeros(N,1)

for k1 = 1:N
    filename = files(k1).name;
    % import
    spectrum = readSPF2(filename);
    % get max val
    [val,idx] = max(spectrum(:,2));
    wvcal(k1) = spectrum(idx,1);
    wvint(k1) = val;
end

wvint = wvint./max(wvint)

figure
plot(wvcal,wvint)
xlabel('Wavelength (nm)')
ylabel('Intensity')