function sAlgName= algName(iAlg)
if iAlg == 0,
    % Matched filter
    sAlgName = 'matched';
elseif iAlg == 2,
    % Zero forcing
    sAlgName = 'zforcing';
else
     sAlgName = '';
end