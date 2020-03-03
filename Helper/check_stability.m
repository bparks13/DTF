function isStable=check_stability(filtering)
%% isStable=check_stability(filtering)
%
%  Automatically goes through the entire filtering struct and checks that each set of
%  numerators and denominators are stable
%
%   Inputs:
%    - filtering: Struct containing filtering parameters
%    -- hpf: If defined, should be a struct containing num and den for a high pass filter. 
%           Default is to use a 4th order butterworth filter, with a cutoff of 1 Hz
%    -- comb: If defined, should be a struct containing num and den for a comb filter. 
%           Default is to use a 60 Hz comb filter, with a qFactor of 35. If this field is
%           empty, the default comb is used
%    -- lpf: If defined, should be a struct containing num and den for a low pass filter. 
%           No default
%    -- notch: If defined, should be a struct containing num and den for a notch filter. 
%           No default. Can contain more than one set of num and den for multiple notches
%    -- ma: If defined, should be an integer defining the number of samples to run
%           through the moving average
%
%   Outputs:
%    - isStable: Boolean defining whether the entire struct is stable or not. If any
%       individual filter is unstable the error message will define which one failed
%

isStable=true;

if isstruct(filtering)
    if isfield(filtering,'hpf')
        if ~isstable(filtering.hpf.num,filtering.hpf.den)
            isStable=false;
            warning('High pass filter is unstable');
            return
        end
    end
    
    if isfield(filtering,'comb')
        if ~isstable(filtering.comb.num,filtering.comb.den)
            isStable=false;
            warning('Comb filter is unstable');
            return
        end
    end

    if isfield(filtering,'lpf')
        if ~isstable(filtering.lpf.num,filtering.lpf.den)
            isStable=false;
            warning('Low pass filter is unstable');
            return
        end
    end

    if isfield(filtering,'notch')
        for i=1:length(filtering.notch)
            if ~isstable(filtering.notch.num,filtering.notch.den)
                isStable=false;
                warning('Notch filter #%d is unstable',i);
                return
            end
        end
    end

    if isfield(filtering,'ma')
       if ~isstable(filtering.ma.num,filtering.ma.den)
            isStable=false;
            warning('Moving average filter is unstable');
            return
        end
    end
end

end