###########################################################################
#' @title Calculate coancestries
#'
#' @description Given genetic data it provides coancestry for between and within populations
#' 
#' @usage calculate_coancestries(genetic_data_parents, genotyped_parent_populations, 
#'   genetic_data_F1 = NA, population_individual_id = NA, 
#'   column_individual = "id", column_population = "population_id", 
#'   pedigree = NA, all_parents_genotyped = FALSE)
#'
#' @param genetic_data_parents A data frame or file name containing genetic data for the parental generation. Must be dosage format.
#' 
#' @param genotyped_parent_populations A vector of population IDs for the parental generation.
#'
#' @param genetic_data_F1 A data frame or file name containing genetic data for the F1 generation. Must be dosage format.
#' 
#' @param population_individual_id (Optional) A data frame mapping individuals to their respective populations. First parent ids, then F1 ids.
#' If not provided, a balanced design is assumed.
#'
#' @param column_individual The name of the column containing individual IDs. Default is "id".
#'
#' @param column_population The name of the column containing population IDs. Default is "population_id".
#'
#' @param pedigree (Optional) A data frame containing pedigree information with at least three columns: 
#' individual ID, dam ID, and sire ID.
#'
#' @param all_parents_genotyped Logical; if TRUE, indicates that the genotyped individuals in "genotyped_parent_populations" are the parents of F1 individuals
#' 
#' @return A list containing:
#' \item{The.M}{A kinship-based relatedness matrix adjusted for F1 individuals.}
#' \item{Theta.P}{An adjusted population coancestry matrix.}
#'
#' @details This function calculates kinship and coancestry measures from genetic data 
#' or pedigree-based information, adjusting for population structure and F1 generation characteristics.
#'
#' @author Isabela do O \email{isabela.doo@@unil.ch}
#' @references
#' - Goudet & Weir (2023)
#'
#' @export

calculate_coancestries <- function(genetic_data_parents,
                                   genotyped_parent_populations, 
                                   genetic_data_F1 = NA, 
                                   population_individual_id = NA, 
                                   column_individual = "id", 
                                   column_population = "population_id", 
                                   pedigree = NA,
                                   all_parents_genotyped = FALSE) {

    # Fix for R CMD check global variable notes
    dam <- dam_pop <- sire <- sire_pop <- id <- population_id <- n <- sire_n <- dam_n <- NULL



    parent_dosage <- genetic_data_parents
    

    # Load F1 data and get kinship ------------------------------------------------------ 
    if (!identical(pedigree, NA)) {
        # Use kinship2::kinship function
        # The kinship2::kinship function returns kinship coefficients (values 0-0.5)
        # We need to convert to a relatedness matrix (values 0-1) by multiplying by 2
        kinship_F1_raw <- kinship2::kinship(id = pedigree$id, 
                                             dadid = pedigree$sire, 
                                             momid = pedigree$dam)
        # Convert kinship coefficients to relatedness matrix and ensure proper format
        kinship_F1 <- as.matrix(kinship_F1_raw * 2)
        rownames(kinship_F1) <- colnames(kinship_F1) <- pedigree$id

        pop_ids_F1 <- pedigree$dam_pop

        F1_id <- pedigree$id
    } else if (!identical(genetic_data_F1, NA)) {
        F1_dosage <- genetic_data_F1
       
        matching_matrix_F1 <- hierfstat::matching(F1_dosage) 
        kinship_F1 <- hierfstat::beta.dosage(matching_matrix_F1, MATCHING = TRUE)
        F1_id <- population_individual_id[,column_individual]
    } else {
        stop("Missing F1 data, either include pedigree or genetic data")
    }
   
    # ----------------------------------------------------------------------------

    # Determine parental population IDs and F1 population IDs (which includes P) --------------------
    if (!identical(pedigree, NA)) {
        population_names <- unique(pedigree$dam_pop)
        parent_pop_id <- genotyped_parent_populations #may differ from pedigree
        #Calculate population sizes for The.M
        dam_df <- pedigree %>%
            select(id = dam, population_id = dam_pop) %>%
            distinct()
        sire_df <- pedigree %>%
            select(id = sire, population_id = sire_pop) %>%
            distinct()
        id_df <- pedigree %>%
            select(id, population_id = dam_pop) %>%
            distinct()

        parent_pop_id_pedigree <- c(sire_df$population_id, dam_df$population_id)
        F1_pop_id <- id_df$population_id

        dam_per_pop <- dam_df %>%
                count(population_id)
        sire_per_pop <- sire_df %>% 
                count(population_id)
        F1_per_pop <- id_df %>%
                count(population_id)
        sire_per_pop <- sire_per_pop %>%
                rename(sire_n = n)
        dam_per_pop <- dam_per_pop %>%
            rename(dam_n = n)
        parent_per_pop <- full_join(sire_per_pop, dam_per_pop, by = "population_id") %>%
                mutate(
                    sire_n = replace_na(sire_n, 0),
                    dam_n = replace_na(dam_n, 0),
                    parent_n = sire_n + dam_n
                )
        
        population_sizes_F1 <- F1_per_pop$n
        population_sizes_P <- parent_per_pop$parent_n #On pedigree (may differ from genotyped)
        population_sizes <- population_sizes_P + population_sizes_F1 #On pedigree (may differ from genotyped)
    } else if (!identical(population_individual_id, NA)) {
        
            population_names <- unique(population_individual_id[, column_population])
            parent_pop_id <- genotyped_parent_populations
            F1_pop_id <- population_individual_id[, column_population] 
            #Calculate population sizes for The.M
            population_sizes_P <- as.data.frame(table(parent_pop_id))$Freq
            population_sizes_F1 <- as.data.frame(table(F1_pop_id))$Freq
            population_sizes <- population_sizes_P + population_sizes_F1
         
    } else {
            stop("Population data missing. \n
            Provide either population_individual_id \n")
    }
    

    cat("There are ", length(unique(parent_pop_id)), " populations. \n")

    if (exists("pop_ids_F1")) {
        if (!identical(pop_ids_F1,F1_pop_id)) {
            warning("Mismatch between ids in F1 genetic data or pedigree and population identification.\n")
        }
    } else {
        pop_ids_F1 <- F1_pop_id
    }



    #Theta_P calculation -----------------------------------------------------
    #Matching matrix and kinship for parents
    matching_matrix_parents <- hierfstat::matching(parent_dosage)
    kinship_parents <- hierfstat::beta.dosage(matching_matrix_parents, MATCHING = TRUE)
    fst_founders <- hierfstat::fs.dosage(matching_matrix_parents, pop = parent_pop_id, matching = TRUE)

    cat("Calculating Theta.P \n")
    min_Fst <- min(hierfstat::mat2vec(fst_founders$FsM))
    Theta_P <- (fst_founders$FsM - min_Fst) / (1 - min_Fst)
    cat("Theta.P calculated with dimensions", dim(Theta_P), "\n")
    #Theta.P must be positive definite. Corrected here, before The.M is built from it
    #(as in the MONET package). nearPD without corr = TRUE: Theta.P is a coancestry
    #matrix, not a correlation matrix, so its diagonal must not be forced to 1.
    eigenvalues_ThetaP <- eigen(Theta_P)$values
    if (any(eigenvalues_ThetaP < 0)) {
        cat("Theta.P matrix not positive definite. \n")
        cat("Minimum eigenvalue is ", min(eigenvalues_ThetaP), "\n")
        Theta_P <- as.matrix(Matrix::nearPD(Theta_P)$mat)
        eigenvalues_ThetaP_new <- eigen(Theta_P)$values
        cat("Theta.P matrix corrected. New min eigenvalue: ", min(eigenvalues_ThetaP_new), "\n")
    }
    # --------------------------------------------------------------------------


    # The.M calculation -------------------------------------------------------
    cat("Calculating The.M \n")
    
    # Check if population sizes are correct ----------------------
    total_individuals <- sum(population_sizes)
    total_individuals_F1 <- sum(population_sizes_F1)
    #cat("population sizes including F1 and Founders are ", population_sizes, "\n")
    #cat("total individuals are ", total_individuals, "\n")
    cat("population sizes of F1", population_sizes_F1, "\n")
    cat("total individuals are ", total_individuals_F1, "\n")

    

    #Calculate mean kinship per population of parental (founder) population
    unique_pops <- unique(F1_pop_id) #same number of populations in F1 and in P
    #Should be the same as population_names

    

    
    #Explaining following If statement:
    #If you don't have genetic info from the parents, 
    #then you cannot infer their relationship, 
    #so you have to assume it is zero 
    #when you use the pedigree to estimate relatedness in F1
    if(all_parents_genotyped == FALSE && !identical(pedigree, NA)){
        adjusted_kinship_F1 <- kinship_F1
    } else {
        #Adjust kinship for F1 individuals - we use the mean kinship of the founders to standardize M
        adjusted_kinship_F1 <- matrix(0,nrow=nrow(kinship_F1),ncol=nrow(kinship_F1))
        last_individual_counted <- 0
        mean_kinship_per_population <- sapply(unique_pops, function(pop) {
        pop_indices_P <- which(parent_pop_id == pop)
        mean(hierfstat::mat2vec(kinship_parents[pop_indices_P, pop_indices_P]))
        })
        
        for (pop in 0:(length(unique_pops) - 1)) {  #Version before May 14th 2025 allowed for dosage of phenotyped individuals to include founders. 
        # Find indices for this population's F1
            cat("pop: ", pop, "\n")

            #Adjusting F1 individuals' kinship matrix
            start_idx <- last_individual_counted + 1
            end_idx <- last_individual_counted + population_sizes_F1[pop+1]

                
            adjusted_kinship_F1[start_idx:end_idx, start_idx:end_idx] <- 
                (kinship_F1[start_idx:end_idx, start_idx:end_idx] - mean_kinship_per_population[pop+1]) / 
                (1 - mean_kinship_per_population[pop+1])

            #Update the last individual count tracker
            last_individual_counted <- end_idx
            }

    }

    
    The.M <- matrix(0, nrow = sum(population_sizes_F1), ncol = sum(population_sizes_F1))
    row.names(The.M) <- colnames(The.M) <- F1_id


    for (pop in unique_pops) {
    #Grab the indices for individuals from current population
        pop_indices <- which(pop_ids_F1 == pop)
        dim_current_population <- length(pop_indices)
        cat("There are", dim_current_population, "individuals in population", pop, "for calculating the M matrix\n")
        kin_block <- hierfstat::kinship2grm(adjusted_kinship_F1)[pop_indices, pop_indices]
        block_adjusted <- kin_block * (1 - Theta_P[pop, pop])
        The.M[pop_indices, pop_indices] <- block_adjusted
    }

    #The M must be positive definite (nearPD without corr = TRUE, as in the MONET package)
    eigenvalues <- eigen(The.M)$values
    if (any(eigenvalues < 0)) {
        cat("M matrix not positive definite. \n")
        cat("Minimum eigenvalue is ", min(eigenvalues), "\n")
        The.M <- as.matrix(Matrix::nearPD(The.M)$mat)
        eigenvalues_new <- eigen(The.M)$values
        cat("M matrix corrected. New min eigenvalue: ", min(eigenvalues_new), "\n")
    }

    cat("The.M calculated with dimensions ", dim(The.M), "\n")

    return(list(The.M = The.M, Theta.P = Theta_P))
}
