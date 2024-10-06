# Bayesian-Variable-Selection

# 1. Horseshoe prior 

The horse-shoe prior is a strong prior, used for variable selection, intending local-global shrinkage. Refer to this paper for more information 

https://link.springer.com/article/10.1007/s10260-023-00727-9
  


# 2. Bayesian Lasso

Lasso is indeed a popular technique for variable selection , but its functional form makes optimization of the likelihood difficult, even in the bayesian setting, hence a hiearchical model was employed to achieve the same likelihood (as for lasso), but with a greater flexibility of sampling from the  rather "well known" posterior distribution obtained as a result. This also makes it easy for obtaining the hyper parameter "lambda", which was previously obtained using cross validation. Refer to this paper for more information 

https://tandfonline.com/doi/abs/10.1198/016214508000000337

# 3. Kuo-Mallik
 
 Along with finding which variable is important, this method also lets you know what is the probability of inclusion of each variable , is in the model. Based on these probabilities, the researcher can now judge if the variable is at all important or not.
 
 This refers to the paper : Kuo, L. and Mallick, B., 1998. Variable selection for regression models. Sankhyā: The Indian Journal of Statistics, Series B, pp.65-81.
 
 Link to the paper : https://www.jstor.org/stable/25053023

# 4. Slab and Spike

Variable selection through local-global shrinkage , is what slab and spike prior is used for. Refer to this paper for more information 

https://academic.oup.com/biomet/article-abstract/97/2/465/219397?redirectedFrom=fulltext
