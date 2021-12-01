function playtone(freq,amplitude,duration)
%Plays a simple tone.
%freq = frequency of tone (in Hz).
%sf = sampling frequency (in Hz).
%amplitude = sound amplitude (dimensionless).
%duration = sound duration (in seconds).

if nargin==0
    freq=1000;
    amplitude=0.08;
    duration=0.3;
end

sf=freq*10; %sampling frequency
t = 0:1/sf:duration;
sound_vector = amplitude*sin(2*pi*freq*t);
sound(sound_vector,sf)