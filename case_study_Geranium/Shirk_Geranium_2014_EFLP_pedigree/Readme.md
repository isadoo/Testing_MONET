# Geranium case study

- We keep the two experiments separate becuase they are a different set of individuals, so we use a different set of M and Theta_P matrices. 

- We also have two models
    -- one that is "country-only"
    `Y ~ country + (1 | gr(pop, cov = two.Theta.P)) + (1 | gr(ind, cov = M))`
    This one is the analogue of their Q_CT - F_CT

    -- one that has interaction
    `Y ~ country + (1 | gr(pop, cov = two.Theta.P, by = country)) + (1 | gr(ind, cov = M))
    V_AB is estimated for each country, but there is a single V_AW for both.
    This one is the analogue to their Q_SC - F_SC

- We also have two scripts per model, one of them is the basic version, and the other has `_marginal`
    -- marginal: 
    We have this verison because it's faster, and it works for Gaussian traits, which is the case for many of the traits here, and when M is a block diagonal for population which is the case for this metapopulation.

    M is block-diagonal by population, so each population block is rotated to the eigenbasis of M, where observations are independent with variance sd_A^2 * lambda + sigma^2; w = U' 1 carries the intercept, the country effect and the population effect. This removes the individual effects from the sampler. Integrating them out is exact only when they and the residuals are both gaussian on the scale of the observations, so the Poisson (boltA-C) and Student-t (totC, totN) traits are fitted with the full model.


