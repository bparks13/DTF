function c=load_cortical_contact_coordinates(subjID)
%% c=load_cortical_contact_coordinates(subjID)
%
%  Load the coordinates of the cortical contacts as they are found in freesurfer
%
%   Inputs:
%    - subjID: String containing the ID of the subject (either their absolute ID or the
%       relative ID, 'ET_CL_002' or 'S1', respectively)
%
%   Outputs:
%    - c: Matrix of coordinates corresponding to the XYZ coordinates of the cortical
%       contacts found on the CT scan. Size is [n x 3], where n is the number of contacts
%       and the three dimensions are in order, XYZ
%

if strcmp(subjID,'ET_CL_002') || strcmp(subjID,'S1')
    c = [-37.16, -35.27, 47.67;
        -37.33, -45.39, 45.47;
        -36.50, -55.23, 42.75;
        -35.16, -64.96, 39.48];
elseif strcmp(subjID,'ET_CL_004') || strcmp(subjID,'S2')
    c = [-42.94, -33.01, 56.25;
        -44.82, -43.21, 55.20;
        -46.38, -53.34, 52.94;
        -48.31, -62.87, 51.57];
elseif strcmp(subjID,'ET_OR_018') || strcmp(subjID,'S3')
    c = [-35.84, 3.21, 59.86;
        -34.76, -4.59, 66.71;
        -32.69, -12.86, 71.00;
        -29.04, -21.05, 74.39;
        -25.23, -30.02, 77.24;
        -21.06, -39.46, 78.29];
end

end