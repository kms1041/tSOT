function tSOTsolution = tSOT_v1(model, data, productRxn)
% Transcriptomics-based strain optimization tool [Minsuk Kim et al, Biotech Bioeng, Under revision]
%
%INPUTS
% model             COBRA model structure
% data              Gene expression data structure
%       genes                   Cell containing gene names
%       transcriptomics         Gene expression levels
%       genes and transcriptomics should have same length
% productRxn        Exchange reaction for target product
%
%%OUTPUTS
% tSOTsolution      tSOT solution structure (overexpression targets and expected yield improvement)

    IMAT_LOWER_QUANTILE = 0.5;
    IMAT_UPPER_QUANTILE = 0.75;
    eps_param = 0.0001;
    
    lower_threshold = quantile(data.transcriptomics, IMAT_LOWER_QUANTILE);
    upper_threshold = quantile(data.transcriptomics, IMAT_UPPER_QUANTILE);
    scale_value = abs(model.lb(strcmp(model.rxns, 'EX_glc(e)')));
    eps_param = eps_param * scale_value;
    
    iMATresult = call_iMAT(model, data.genes, data.transcriptomics, lower_threshold, upper_threshold, eps_param);
    offStateRxns = setdiff(iMATresult.lowlyExpressedRxns, iMATresult.upregulatedRxns);
    
    referenceModel = buildReferenceModel(model, offStateRxns);
    referenceModel = changeObjective(referenceModel, productRxn, 1);

    [targetRxns, yieldImprovements] = findTargetRxns(referenceModel, model, offStateRxns);
    
    tSOTsolution.targetRxns = targetRxns;
    tSOTsolution.yieldImprovements = yieldImprovements;

end


function referenceModel = buildReferenceModel(model, offStateRxns)
    referenceModel = model;
    
    for i = 1:length(offStateRxns)
        referenceModel = changeRxnBounds(referenceModel, offStateRxns(i), 0, 'b');
    end
end

function [targetRxns, yieldImprovements] = findTargetRxns(referenceModel, model, offStateRxns)
    targetRxns = {};
    yieldImprovements = [];
    
    solution = optimizeCbModel(referenceModel);
    referenceYield = solution.f;
    
    for i = 1:length(offStateRxns)
        targetUpregulatedModel = referenceModel;
        rxnID = find(strcmp(model.rxns, offStateRxns(i)));
        targetUpregulatedModel.lb(rxnID) = model.lb(rxnID);
        targetUpregulatedModel.ub(rxnID) = model.ub(rxnID);
        solution = optimizeCbModel(targetUpregulatedModel);
        if solution.f/referenceYield > 1.001
            targetRxns = [targetRxns; offStateRxns(i)];
            yieldImprovements = [yieldImprovements; (solution.f/referenceYield - 1) * 100];
        end
    end
end


function iMATresult = call_iMAT(model, gene_names, gene_exp, lower_threshold, upper_threshold, eps_param)
% Implements iMAT as defined in [Shlomi et al, Nat Biotech, 2008].
% Adapted from the implementation provided in the cobra toolbox.
%
% INPUTS
%       model - cobra model
%       gene_names - genes ids
%       gene_exp - genes expression
%       lower_threshold - lower expression threshold
%       upper_threshold - upper expression threshold
%       eps_param - flux activation threshold
%
% OUTPUTS
%       model_exp - model integrated with transcriptomics data
%
% Original author: Daniel Machado, 2013, [Machado and Herrgard, PLOS Comput Biol, 2014]
%
% Modified for tSOT simulation : Minsuk Kim, 2014

    discrete_levels = zeros(size(gene_exp));
    discrete_levels(gene_exp > upper_threshold) = 1;
    discrete_levels(gene_exp < lower_threshold) = -1;
    reaction_levels = gene_to_reaction_levels(model, gene_names, discrete_levels, @min, @max);
    RHindex = find(reaction_levels > 0);
    RLindex = find(reaction_levels < 0);
    
    [~, iMATresult] = shlomi(model, RHindex, RLindex, eps_param);

end

function [v_sol, iMATresult] = shlomi(model, RHindex, RLindex, eps_param)
% Implementation from the cobra toolbox (createTissueSpecificModel.m)
% Modified for tSOT simulation : Minsuk Kim, 2014

    lowlyExpressedRxns = model.rxns(RLindex);
    upregulatedRxns = {};
    
    S = model.S;
    lb = model.lb;
    ub = model.ub;

    % Creating A matrix
    A = sparse(size(S,1)+2*length(RHindex)+2*length(RLindex),size(S,2)+2*length(RHindex)+length(RLindex));
    [m,n,s] = find(S);
    for i = 1:length(m)
        A(m(i),n(i)) = s(i); %#ok<SPRIX>
    end

    for i = 1:length(RHindex)
        A(i+size(S,1),RHindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1),i+size(S,2)) = lb(RHindex(i)) - eps_param; %#ok<SPRIX>
        A(i+size(S,1)+length(RHindex),RHindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1)+length(RHindex),i+size(S,2)+length(RHindex)+length(RLindex)) = ub(RHindex(i)) + eps_param; %#ok<SPRIX>
    end

    for i = 1:length(RLindex)
        A(i+size(S,1)+2*length(RHindex),RLindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1)+2*length(RHindex),i+size(S,2)+length(RHindex)) = lb(RLindex(i)); %#ok<SPRIX>
        A(i+size(S,1)+2*length(RHindex)+length(RLindex),RLindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1)+2*length(RHindex)+length(RLindex),i+size(S,2)+length(RHindex)) = ub(RLindex(i)); %#ok<SPRIX>
    end

    % Creating csense
    csense1(1:size(S,1)) = 'E';
    csense2(1:length(RHindex)) = 'G';
    csense3(1:length(RHindex)) = 'L';
    csense4(1:length(RLindex)) = 'G';
    csense5(1:length(RLindex)) = 'L';
    csense = [csense1 csense2 csense3 csense4 csense5];

    % Creating lb and ub
    lb_y = zeros(2*length(RHindex)+length(RLindex),1);
    ub_y = ones(2*length(RHindex)+length(RLindex),1);
    lb = [lb;lb_y];
    ub = [ub;ub_y];

    % Creating c
    c_v = zeros(size(S,2),1);
    c_y = ones(2*length(RHindex)+length(RLindex),1);
    c = [c_v;c_y];

    % Creating b
    b_s = zeros(size(S,1),1);
    lb_rh = lb(RHindex);
    ub_rh = ub(RHindex);
    lb_rl = lb(RLindex);
    ub_rl = ub(RLindex);
    b = [b_s;lb_rh;ub_rh;lb_rl;ub_rl];

    % Creating vartype
    vartype1(1:size(S,2),1) = 'C';
    vartype2(1:2*length(RHindex)+length(RLindex),1) = 'B';
    vartype = [vartype1;vartype2];

    MILPproblem.A = A;
    MILPproblem.b = b;
    MILPproblem.c = c;
    MILPproblem.lb = lb;
    MILPproblem.ub = ub;
    MILPproblem.csense = csense;
    MILPproblem.vartype = vartype;
    MILPproblem.osense = -1;
    MILPproblem.x0 = [];

    params.timeLimit = 100;
    solution = solveCobraMILP(MILPproblem, params);

    x = solution.cont;
    for i = 1:length(x)
        if abs(x(i)) < 1e-6
            x(i,1) = 0;
        end
    end

    v_sol = x;
    
    for i = 1:length(model.rxns)
        if length(find(RLindex == i))
            if ~(solution.int(find(RLindex == i) + length(RHindex)) == 1)
                upregulatedRxns = [upregulatedRxns; model.rxns(i)];
            end
        end
    end

    iMATresult.lowlyExpressedRxns = lowlyExpressedRxns;
    iMATresult.upregulatedRxns = upregulatedRxns;

end

function reaction_levels = gene_to_reaction_levels( model, genes, levels, f_and, f_or )
% Convert gene expression levels to reaction levels using GPR associations.
% Level is NaN if there is no GPR for the reaction or no measured genes.
%
% INPUTS
%       model - cobra model
%       genes - gene names
%       levels - gene expression levels
%       f_and - function to replace AND
%       f_or - function to replace OR
%
% OUTPUTS
%       reaction_levels - reaction expression levels
%
% Author: Daniel Machado, 2013

    reaction_levels = zeros(length(model.rxns), 1);

    for i = 1:length(model.rxns)
        level = eval_gpr(model.grRules{i}, genes, levels, f_and, f_or);
        reaction_levels(i) = level;
    end

end

function [result, status] = eval_gpr(rule, genes, levels, f_and, f_or)
% Evaluate the expression level for a single reaction using the GPRs.
% Note: Computes the expression level even if there are missing measured
% values for the given rule. This implementation is a modified version of
% an implementation provided in [Lee et al, BMC Sys Biol, 2012]

    EVAL_OK = 1;
    PARTIAL_MEASUREMENTS = 0;
    NO_GPR_ERROR = -1;
    NO_MEASUREMENTS = -2;
    MAX_EVALS_EXCEEDED = -3;

    MAX_EVALS = 1000;
    NONETYPE = 'NaN';

    NUMBER = '[0-9\.\-e]+';
    MAYBE_NUMBER = [NUMBER '|' NONETYPE];

    expression = rule;
    result = NaN;
    status = EVAL_OK;

    if isempty(expression)
        status = NO_GPR_ERROR;
    else
        rule_genes = setdiff(regexp(expression,'\<(\w|\-)+\>','match'), {'and', 'or'});
        
        total_measured = 0;
        
        for i = 1:length(rule_genes)
            j = find(strcmp(rule_genes{i}, genes));
            if isempty(j)
                level = NONETYPE;
            else
                level = num2str(levels(j));
                total_measured = total_measured + 1;
            end
            expression = regexprep(expression, ['\<', rule_genes{i}, '\>'], level );
        end
        
        
        if total_measured == 0
            status = NO_MEASUREMENTS;
        else
            if total_measured < length(rule_genes)
                status = PARTIAL_MEASUREMENTS;
            end
            
            maybe_and = @(a,b)maybe_functor(f_and, a, b);
            maybe_or = @(a,b)maybe_functor(f_or, a, b); 
            str_wrapper = @(f, a, b)num2str(f(str2double(a), str2double(b)));

            counter = 0;
            
            while isnan(result)

                counter = counter + 1;
                if counter > MAX_EVALS
                    status = MAX_EVALS_EXCEEDED;
                    break
                end

                try 
                    result = eval(expression);            
                catch e   
                    paren_expr = ['\(\s*(', MAYBE_NUMBER,')\s*\)'];
                    and_expr = ['(',MAYBE_NUMBER,')\s+and\s+(',MAYBE_NUMBER,')'];
                    or_expr = ['(',MAYBE_NUMBER,')\s+or\s+(',MAYBE_NUMBER,')'];

                    expression = regexprep(expression, paren_expr, '$1');
                    expression = regexprep(expression, and_expr, '${str_wrapper(maybe_and, $1, $2)}');
                    expression = regexprep(expression, or_expr, '${str_wrapper(maybe_or, $1, $2)}');
                end
            end
            
        end
    end

end

function c = maybe_functor(f, a, b)
    
    if isnan(a) && isnan(b)
        c = nan;
    elseif ~isnan(a) && isnan(b)
        c = a;
    elseif isnan(a) && ~isnan(b)
        c = b;
    else 
        c = f(a,b);
    end
end