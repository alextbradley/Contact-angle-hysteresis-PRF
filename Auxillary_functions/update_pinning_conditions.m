function is_pinned = update_pinning_conditions(is_pinned, sol)
%updates the pinning conditions based on the previous integration, whose
%solution structure is given by sol. This function is simply an exhaustive
%list.

if all(is_pinned == [0,0]) %both pinned
    if sol.ie == 1 %theta_l passes thru theta_a: set xl to advancing
        is_pinned(1) = 2;
    elseif sol.ie == 2  %theta_- passes thru theta_r: set xl to receding
        is_pinned(1) = 1;
    elseif sol.ie == 3 %theta_+ passes thru theta_a: set xl to advancing
        is_pinned(2) = 2;
    elseif sol.ie == 4 %theta_+ passes thru theta_r: set xl to receding
        is_pinned(1) = 1;
    end
    
elseif all(is_pinned == [1,0]) %xl receding, xu pinned
    if sol.ie == 1 %xl speed thru 0: set xl to pinned
        is_pinned(1) = 0;   
    elseif sol.ie == 3 %theta_+ passes thru theta_a: set xl to advancing
        is_pinned(2) = 2;
    elseif sol.ie == 4 %theta_+ passes thru theta_r: set xl to receding
        is_pinned(1) = 1;
    end
    
elseif all(is_pinned == [2,0]) %xl advancing, xu pinned
    if sol.ie == 1 %xl speed thru 0: set xl to pinned
        is_pinned(1) = 0;   
    elseif sol.ie == 3 %theta_+ passes thru theta_a: set xl to advancing
        is_pinned(2) = 2;
    elseif sol.ie == 4 %theta_+ passes thru theta_r: set xl to receding
        is_pinned(1) = 1;
    end
    
elseif all(is_pinned == [0,1]) %xl pinned, xu receding
    if sol.ie == 1 %theta_l passes thru theta_a: set xl to advancing
        is_pinned(1) = 2;
    elseif sol.ie == 2  %theta_- passes thru theta_r: set xl to receding
        is_pinned(1) = 1;
    elseif sol.ie == 3 %xu speed thru 0, so xu now pinned
        is_pinned(2) = 0;
    end   
    
elseif all(is_pinned == [0,2]) %xl pinned, xu advancing
    if sol.ie == 1 %theta_l passes thru theta_a: set xl to advancing
        is_pinned(1) = 2;
    elseif sol.ie == 2  %theta_l passes thru theta_r: set xl to receding
        is_pinned(1) = 1;
    elseif sol.ie == 3 %xu velocity thru 0: set xu to pinned
        is_pinned(2) = 0;
    end
    
elseif all(is_pinned == [1,1]) %xl receding, xu receding
    if sol.ie == 1 %xl speed thru 0: set xl to pinned
        is_pinned(1) = 0;
    elseif sol.ie == 3 %xu speed thru 0, so xu now pinned
        is_pinned(2) = 0;
    end   
        
        
elseif all(is_pinned == [1,2]) %xl receding, xu advancing
    if sol.ie == 1 %xl speed thru 0, so xl now pinned
        is_pinned(1) = 0;
    elseif sol.ie == 3 %xu speed thru 0, so xu now pinned
        is_pinned(2) = 0;
    end     
  
elseif all(is_pinned == [2,1]) %xl advancing, xu receding
    if sol.ie == 1 %xl speed thru 0: set xl to pinned
        is_pinned(1) = 0;
    elseif sol.ie == 3 %xu speed thru 0, so xu now pinned
        is_pinned(2) = 0;
    end
    
elseif all(is_pinned == [2,2]) %xl advancing, xu advancing
    if sol.ie == 1 % xl speed thru 0, so xl now pinned
        is_pinned(1) = 0;
    elseif sol.ie == 3 %xu speed thru 0, so xu now pinned
        is_pinned(2) = 0;
    end
end
