function[modelTAC,K,RSS] = SIME_convolve(inFcn,wb,inFcn_t,tac,tac_t,model,Vnd,weights,Settings)

% Performs the convolution between a metabolite corrected
% arterial input function and a impuls response function (IRF)defined by the
% kinetic model. The convolution is performed numerically, using the matlab
% command "filter". The model TAC is returned as well as a vector K
% containing all rate constants included in the IRF.
% _________________________________________________________________________
%                                                   Martin Schain, CU, 2017

%% Do the convolution
options = optimset('display','off'); options.MaxFunEvals = 1000;
switch model
                
    case '2TCM'
        startParameters = [.01 .01 .01 .01]; %[K1 k2 k3 k4], unit = 1/min
        lowerBound      = [0 0 0 0];
        upperBound      = [inf inf inf inf];
        Vnd             = [];
        
    case '2TCMfixedVnd'
        if nargin < 6, error('Fixed value for Vnd is not provided'), end
        startParameters = [.1 .1 .1]; %[k2 k3 k4], unit = 1/min        
        lowerBound      = [0 0 0];
        upperBound      = [inf inf inf];        
end

% If vB should be fitted, add it to the start paramaters. 
if ~isscalar(Settings.vB)    
    startParameters(end+1) = 0.05;
    lowerBound(end+1) = 0;
    upperBound(end+1) = inf;
end

% Check if the input function has even time steps
if std(diff(inFcn_t)) > 2
    error('Not equally sampled inFcn')
end

[K,RSS]  = lsqnonlin(@convWithIRF_fit,startParameters,lowerBound,upperBound,options,inFcn,wb,inFcn_t,tac,tac_t,weights,model,Vnd,Settings);

%% Build the model TAC
irf  = SIME_getIRF(K,inFcn_t,model,Vnd);
stepSize = inFcn_t(2) - inFcn_t(1);
if isscalar(Settings.vB)
    vB = Settings.vB;
    K(end+1) = vB;
else
    vB = K(end);
end    
model_int = (1-vB)*stepSize*filter(inFcn,1,irf) + vB*wb;
modelTAC  = interp1(inFcn_t,model_int,tac_t,'pchip');

    
function err = convWithIRF_fit(K,inFcn,wb,inFcn_t,tac,tac_t,weights,model,Vnd,Settings)

    irf       = SIME_getIRF(K,inFcn_t,model,Vnd);
    stepSize  = inFcn_t(2)-inFcn_t(1);
    
    if isscalar(Settings.vB)
        vB = Settings.vB;
    else
        vB = K(end);
    end
    
    model_int = (1-vB)*stepSize*filter(inFcn,1,irf) + vB*wb;
    err       = ( interp1(inFcn_t,model_int,tac_t,'pchip') - tac ).*weights(:);