##' @title
##' PmatrixPosterior
##'
##' @description
##' Estimate a posterior distribution of P-matrices for a given population and linear model
##'
##' @param data matrix of individual data
##' @param model character string with model specifications, referring to columns in info
##' @param info fixed effects
##' @param modelpath path to stan file
##' @param full.output should output full fitted model object?
##'                    defaults to FALSE (outputs only posterior P)
##' @param ... possibly extra parameters for stan model fitting
##'
##' @return array with P-matrix posterior distribution
##'
##' @importFrom rstan stan extract
##' @importFrom evolqg RandomMatrix
##'
##' @author Guilherme Garcia
##'
PmatrixPosterior <- function(data, model, info,
                             modelpath = '../Stan/pmatrix_lm.stan',
                             full.output = FALSE, ...)
    {
        ## type III sum of squares
        options(contrasts = c('contr.sum', 'contr.poly'))
        
        initialConditions <-
            function(i, traits, effects)
                list('Omega_P' = chol(RandomMatrix(traits)),
                     'sigma_P' = rchisq(traits, 1),
                     'beta' = matrix(rnorm(traits * effects, 0, 1), nrow = traits))

        if(model == 'NONE')
            form <- as.formula('~ 1')
        else
            form <- as.formula(paste('~', model))

        modmat <- model.matrix(form, data = info)

        mod.input <- list('N' = nrow(data),
                          'J' = ncol(modmat),
                          'K' = ncol(data),
                          'Y' = data,
                          'X' = modmat)

        modelfit <-
            stan(file = modelpath,
                 data = mod.input,
                 pars = c('beta', 'P'), chains = 1,
                 init = list(initialConditions(1, 
                                               traits = mod.input $ K,
                                               effects = mod.input $ J)), ...)
        
        posterior <- rstan::extract(modelfit)

        if(full.output) return(posterior) else return(posterior $ P)
    }
