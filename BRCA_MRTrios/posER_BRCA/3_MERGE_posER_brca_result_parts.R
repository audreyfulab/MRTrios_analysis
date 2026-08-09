library(data.table)

## Ondemand path:
# input_dir <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_posER_BRCA_part"
# output_dir <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_posER_BRCA"

## Local path:
input_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA_part"
output_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA"

# ── Find all shard files (exclude merged ALL file) ─────────
shard_files <- sort(list.files(input_dir,
                               pattern = "trio_posER_BRCA_results_part.*\\.txt",
                               full.names = TRUE))

cat(sprintf("Found %d shard files. Merging...\n", length(shard_files)))
##  600 shard files

# ── Merge all shards ───────────────────────────────────────
all_results <- rbindlist(lapply(shard_files, fread), fill = TRUE)

cat(sprintf("Total rows: %d | Total cols: %d\n",
            nrow(all_results), ncol(all_results)))

## Total rows: 292525 | Total cols: 24

# ── Save merged file ───────────────────────────────────────
out_file <- file.path(output_dir, "trio_posER_BRCA_results_ALL.txt")
fwrite(all_results, file = out_file, sep = "\t")

cat(sprintf("Saved: %s\n", out_file))

## Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BRCA/Output_posER_BRCA/trio_posER_BRCA_results_ALL.txt

# ====================================================================================================
## Merge info:
# # 08/04/26 Fri - 8/9/26 Sun : Ondemand done with
# 1-10 array: 1-60 parts
## Merge above parts:
# 600 shard files:
## Total rows: 292525 | Total cols: 24

# colnames(all_results)
# [1] "trio_row"       "meth.row"       "cna.row"        "gene.row"       "Total.PC.Count"
# [6] "Confounders"    "b11"            "b12"            "b21"            "b22"           
#[11] "V1:T1"          "V1:T2"          "pb11"           "pb12"           "pb21"          
#[16] "pb22"           "pV1:T1"         "pV1:T2"         "Minor.freq"     "coef11"        
#[21] "coef12"         "coef21"         "coef22"         "Inferred.Model"

# table(all_results$Confounders)

# age|race|U_exp|U_meth       age|race|U_meth 
#               292372                   153 
# ====================================================================================================
