function [ThermistorFiltered] = FilterThermistor(Thermistor,varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('samplerate', 500, @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
SampleRate = params.Results.samplerate;

% default settings for filtering
fband = [0.1 30];
Np    = 4; % filter order
[b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients

% delaing with nans
while any(isinf(Thermistor))
    x1 = find(isinf(Thermistor),1,'first') - 1;
    x2 = x1 + find(~isinf(Thermistor(x1+1:end)),1,'first');
    Thermistor(x1:x2) = linspace(Thermistor(x1),Thermistor(x2),(x2-x1+1));
end

% % delaing with nans
% if any(isinf(Thermistor))
%     x1 = find(isinf(Thermistor),1,'first') - 1;
%     x2 = find(isinf(Thermistor),1,'last')  + 1;
%     Thermistor(x1:x2) = linspace(Thermistor(x1),Thermistor(x2),(x2-x1+1));
% end

% filtering - thermistor
ThermistorFiltered = filtfilt(b,a,Thermistor);

end