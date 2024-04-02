function [solution] = calcATPAndNADHMin(starchmodel)
    starchpos = find(ismember(starchmodel.rxns,'TEH_Starch'));
    atpasepos = find(ismember(starchmodel.rxns,'ATPase'));
    dehogpos = find(ismember(starchmodel.rxns,'DEHOG'));
    res = optimizeCbModel(starchmodel,'min');  
    starchFixed = changeRxnBounds(starchmodel,'TEH_Starch',res.v(starchpos),'b');
    
    % Set objective to dehog
    starchFixed = changeObjective(starchFixed,'DEHOG',1);    
    
    % Check, if it works without free NADPH
    nadphCheck = changeRxnBounds(starchFixed,'DEHOG',0, 'l');    
    res = optimizeCbModel(nadphCheck,'max');
    if(res.stat == 1) 
        % This is possible, so there is surplus reductant
        %Set ATP objective
        atpcheck = changeObjective(nadphCheck,'ATPase',1);
        % test, what's the max ATP amount that can be produced
        res = optimizeCbModel(atpcheck,'max');
        if res.f > -1e-9
            % We have surplus (or at least no requirement for) ATP as well. 
            % now, don't allow free ATP to be produced
            dehogMinModel = changeRxnBounds(nadphCheck,'ATPase',0,'l');
            mindehogRes = optimizeCbModel(dehogMinModel,'max');            
            atpMinModel = changeRxnBounds(dehogMinModel,'DEHOG',mindehogRes.v(dehogpos),'b');
            atpMinModel = changeObjective(atpMinModel,'ATPase',1);
            solution = optimizeCbModel(atpMinModel,'max');
        else
            % we do need ATP, when no NADPH is available
            % so this is the requirement.
            solution = res;
        end
    else
        % We lack reductant, so optimize the reductant usage
        res = optimizeCbModel(starchFixed,'max');
        % now, set it, and optimize the ATP usage
        dehogFixed = changeRxnBounds(starchFixed,'DEHOG',res.v(dehogpos),'b');
        dehogFixed = changeObjective(dehogFixed,'ATPase',1);
        solution = optimizeCbModel(dehogFixed,'max');    
    end
end