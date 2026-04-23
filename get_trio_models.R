# ===============================================
# Step 1: Prepare trio data
# ===============================================
# Your trio indices
trio_indices <- cor_M1.2$pos_trios_row # replace your trio indices
# trio_indices <- rownames(trios)

# ==============================================
# Step 2: Load the data from get_trio_data.R
# ==============================================

# The same loading from get_trio_data.R

# ==============================================
# Step 3: Write a function: get_trio_models()
# ==============================================

# Loop over all trios and collect results
## write a function:
get_trio_models <- function(trio_indices, trios, pos_meth_filter, pos_cna_filter,
                            pos_exp_filter, pca_E_final, pca_M_final,
                            gene_pos_indice, meth.posER.pc.sig,
                            pos_age, pos_race) {
  
  assign("isTrue", isTRUE, envir = .GlobalEnv)
  
  do.call(rbind, lapply(trio_indices, function(i) {
    df_trio <- get_trio_data(i, trios, pos_meth_filter, pos_cna_filter,
                             pos_exp_filter, pca_E_final, pca_M_final,
                             gene_pos_indice, meth.posER.pc.sig,
                             pos_age, pos_race)
    if (is.null(df_trio)) return(NULL)
    
    model <- as.data.frame(infer.trio(df_trio))
    
    cbind(data.frame(trio_row = i,
                     meth.row = trios[i,]$meth.row,
                     cna.row  = trios[i,]$cna.row,
                     gene.row = trios[i,]$gene.row),
          model)
  }))
}

# ==============================================
# Step 4: Model Application and Validation
# ==============================================

# Original case (Pos_M1.2 model)
corM1.2_models_df <- get_trio_models(cor_M1.2$pos_trios_row, trios, pos_meth_filter,
                                     pos_cna_filter, pos_exp_filter, pca_E_final,
                                     pca_M_final, gene_pos_indice, meth.posER.pc.sig,
                                     pos_age, pos_race)

## Note: by checking the result, some trios belong to M4(or other model) instead of M1.2

# Any other set of trios
corM1.1_models_df <- get_trio_models(cor_M1.1$pos_trios_row,  trios, pos_meth_filter,
                                     pos_cna_filter, pos_exp_filter, pca_E_final,
                                     pca_M_final, gene_pos_indice, meth.posER.pc.sig,
                                     pos_age, pos_race)

## Note: by checking the result, some trios belong to M4 instead of M1.1


### === check all the Pos_models of M1.2 ====
final_pos_M1.2_models <- get_trio_models(final_pos_M1.2$pos_trios_row,  trios, pos_meth_filter,
                                         pos_cna_filter, pos_exp_filter, pca_E_final,
                                         pca_M_final, gene_pos_indice, meth.posER.pc.sig,
                                         pos_age, pos_race)
