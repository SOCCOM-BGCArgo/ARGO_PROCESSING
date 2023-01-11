function [Estimates,Uncertainties]=ESPER_Mixed(varargin)
% Averages estimates from ESPER_NN and ESPER_LIR.  See either subfunction
% for comments or enter 'help ESPER_NN' into the command line for a guide
% on the use of this function.  The input arguments are the same for this
% function and for both subfunctions.
    if nargout==1
        [NNEsts]=ESPER_NN(varargin{:});
        [LIREsts]=ESPER_LIR(varargin{:});
        EstTypes=fieldnames(NNEsts);
    elseif nargout==2
        [NNEsts,NNUncEsts]=ESPER_NN(varargin{:});
        [LIREsts,LIRUncEsts]=ESPER_LIR(varargin{:});
        EstTypes=fieldnames(NNEsts);
    end
    for Property=1:size(EstTypes,1)
        Estimates.(EstTypes{Property,1})=(NNEsts.(EstTypes{Property,1})+LIREsts.(EstTypes{Property,1}))/2;
        if nargout==2;
            Uncertainties.(EstTypes{Property,1})=min(cat(3,NNUncEsts.(EstTypes{Property,1}),LIRUncEsts.(EstTypes{Property,1})),[],3);
        end
    end
end