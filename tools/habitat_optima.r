# Habitat information used as covariate by MONET w/ habitat (monet_whabitat) and by the Driftsel H-test.
# This "habitat" is the selection optimum of each population in the quantinemo simulations,
# in population order (pop 1 to 18), so it must be the same order as the populations in the prepared data.
#
# usage: source("tools/habitat_optima.r")
#        habitat_optima <- get_habitat_optima(Population_structure, correlation, wdiff)

#The optima depend on wdiff. Values copied from the quanti_stab_sel_optima of the quantinemo .ini files.
#(review 2026: the wdiff 10 values were added; before, the wdiff 4.6 values were used for everything)
get_habitat_optima <- function(Population_structure, correlation, wdiff) {

  if (wdiff == "10") {
    if (Population_structure == "hierarchical") {
      if (correlation == "grouped") {
        return(rep(c(-5, -3, -1, 1, 3, 5), each = 3))
      } else if (correlation == "swapped") {
        return(rep(c(-5, -2.5, 0, -2.5, 0, 2.5, 0, 2.5, 5), 2))
      }
    } else if (Population_structure == "SS") {
      if (correlation == "cline") {
        return(c(-5.0000000, -4.4117647, -3.8235294, -3.2352941,
                 -2.6470588, -2.0588235, -1.4705882, -0.8823529,
                 -0.2941176, 0.2941176, 0.8823529, 1.4705882, 2.0588235,
                 2.6470588, 3.2352941, 3.8235294, 4.4117647, 5.0000000))
      } else if (correlation == "sine") {
        return(c(-5.0000000, -4.1731656, -2.7778512, -0.9549150, 1.0355339,
                 2.9644661, 4.6939394, 5.8268344, 6.1731656, 5.8268344, 4.6939394,
                 2.9644661, 1.0355339, -0.9549150, -2.7778512, -4.1731656, -5.0000000, -4.8268344))
      }
    } else if (Population_structure == "IM") {
      return(c(-5.00, -4.41, -3.82, -3.24, -2.65, -2.06, -1.47, -0.88, -0.29,
               0.29,  0.88,  1.47,  2.06,  2.65, 3.24,  3.82,  4.41,  5.00))
    }
    stop(paste0("Invalid parameter combination: Population_structure=", Population_structure,
                ", correlation=", correlation, ", wdiff=", wdiff))
  }

  if (wdiff != "4p6") {
    stop(paste0("No habitat optima defined for wdiff=", wdiff))
  }

  if (Population_structure == "hierarchical") {
    if (correlation == "grouped") {
      return(rep(c(-2.3, -1.38, -0.46, 0.46, 1.38, 2.3), each = 3))
    } else if (correlation == "swapped") {
      return(rep(c(-2.3, -1.15, 0, -1.15, 0, 1.15, 0, 1.15, 2.3), 2))
    }
  } else if (Population_structure == "SS") {
    if (correlation == "cline") {
      return(c(-2.3000000, -2.0294118, -1.7588235, -1.4882353,
               -1.2176471, -0.9470588, -0.6764706, -0.4058824,
               -0.1352941, 0.1352941, 0.4058824, 0.6764706, 0.9470588,
               1.2176471, 1.4882353, 1.7588235, 2.0294118, 2.3000000))
    } else if (correlation == "sine") {
      return(c(-2.3000000, -1.9196562, -1.2778076, -0.4392609, 0.4763456, 
              1.3636544, 2.1592261, 2.6803438, 2.9196562, 2.6803438, 2.1592261, 
              1.3636544, 0.4763456, -0.4392609, -1.2778076, -1.9196562, -2.3000000, -2.2203438))
    }
  } else if (Population_structure == "IM") {
    # For IM, correlation is NA, so we don't check it
    return(c(-2.30, -2.03, -1.76, -1.49, -1.22, -0.95, -0.68, -0.41, -0.14,  
             0.14,  0.41,  0.68,  0.95,  1.22, 1.49,  1.76,  2.03,  2.30))
  }
  
  stop(paste0("Invalid parameter combination: Population_structure=", Population_structure,
              ", correlation=", correlation))
}

