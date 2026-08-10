library(data.table)

## Ondemand path:
input_dir <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_LIHC_part"
output_dir <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_LIHC"

## Local path:
input_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_LIHC/Output_LIHC_part"
output_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_LIHC/Output_LIHC"

# ── Find all shard files (exclude merged ALL file) ─────────
shard_files <- sort(list.files(input_dir,
                               pattern = "trio_LIHC_results_part.*\\.txt",
                               full.names = TRUE))

cat(sprintf("Found %d shard files. Merging...\n", length(shard_files)))

# ── Merge all shards ───────────────────────────────────────
all_results <- rbindlist(lapply(shard_files, fread), fill = TRUE)

cat(sprintf("Total rows: %d | Total cols: %d\n",
            nrow(all_results), ncol(all_results)))

# ── Save merged file ───────────────────────────────────────
out_file <- file.path(output_dir, "trio_results_LIHC_ALL.txt")
fwrite(all_results, file = out_file, sep = "\t")

cat(sprintf("Saved: %s\n", out_file))


# ── OUTPUT INFO: ──────────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────────────────────────────────────────────
#Found 620 shard files. Merging...
#Total rows: 308028 | Total cols: 24
#Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_LIHC/Output_LIHC/trio_results_LIHC_ALL.txt

# table(all_results$Confounders)
 age|sex|race        age|sex|race|U_exp    age|sex|race|U_exp|U_meth       age|sex|race|U_meth
   13                      2416                    302983                     2616

# ────────────────────────────────────────────────────────────────────────────────────────────────────────



