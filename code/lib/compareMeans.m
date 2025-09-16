function [nameAlgo,p_ouput] = compareMeans(feature1,feature2)
%Protocol

%First, we evaluated whether features values of 2 kind of data presented normal 
%distribution and variance using Kolmogorov–Smirnov and F?Snedecor tests, respectively. 
%In case that data presented normal distribution but not an equal variance,
%we employed a two?tailed Welch test to compare the means from both groups. 
%In case that data presented not normal distribution, 
%we employed Wilcoxon test to compare the medians from both groups. 
%We employed a two?tailed Student's t ?test to compare the means in the cases
%where both distribution and variance of the two sets of data were similar
    

    %kolmogorov smirnov test if n>50 -> 
    %h=0 -> normal distribution; h=1 -> reject the null hypothesis at p-value ->0.05*
    %else: shapiro-wilk test

    if length(feature1)<50
        [hKS1,pks1] = swtest(feature1);
    else
%             modFeat1 = (feature1-mean(feature1))/std(feature1);
%             [hKS1,pks1] = kstest(modFeat1);
        [hKS1,pks1] = lillietest(feature1);
    end
    if length(feature2)<50
        [hKS2,pks2] = swtest(feature2);
    else
%             modFeat2 = (feature2-mean(feature2))/std(feature2);
%             [hKS2,pks2] = kstest(modFeat2);
        [hKS2,pks2] = lillietest(feature2);
    end
   
    
    
    if hKS1 == 0 && hKS2 ==0
    
        %   H = VARTEST2(X,Y) performs an F test of the hypothesis that two
        %   independent samples, in the vectors X and Y, come from normal
        %   distributions with the same variance, against the alternative that
        %   they come from normal distributions with different variances.
        %   The result is H=0 if the null hypothesis ("variances are equal")
        %   cannot be rejected at the 5% significance level, or H=1 if the null
        %   hypothesis can be rejected at the 5% level.  X and Y can have
        %   different lengths.
        hVar = vartest2(feature1,feature2);
        
        if hVar==0
            %t-test. Normal distribution and same variances
            [h_ttest, p_ttest]=ttest2(feature1,feature2);
            nameAlgo = 't-test';
            p_ouput = p_ttest;
        else
            %Satterthwaite's approximation. This test is sometimes
            %called Welch’s t-test. Normal distribution and different
            %variances.
            [h_wtest,p_wtest] = ttest2(feature1, feature2, 'Vartype','unequal');
            nameAlgo = 'Welch';
            p_ouput = p_wtest;
        end
    
    else
        % two-sided Wilcoxon rank sum test. ranksum tests the null hypothesis 
        %that data in x and y are samples from continuous distributions 
        %with equal medians, against the alternative that they are not.
        
        %The Wilcoxon rank sum test is equivalent to the Mann-Whitney U test.
        [p_Wilcoxon,h_Wilcoxon] = ranksum(feature1,feature2);
        nameAlgo = 'Wilcoxon';
        p_ouput = p_Wilcoxon;
    end

end

