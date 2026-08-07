library(data.table)

input_dir <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_posER_BRCA_part"
output_dir <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_posER_BRCA"

# ── Find all shard files (exclude merged ALL file) ─────────
shard_files <- sort(list.files(input_dir,
                               pattern = "trio_posER_BRCA_results_part.*\\.txt",
                               full.names = TRUE))

cat(sprintf("Found %d shard files. Merging...\n", length(shard_files)))

# ── Merge all shards ───────────────────────────────────────
all_results <- rbindlist(lapply(shard_files, fread), fill = TRUE)

cat(sprintf("Total rows: %d | Total cols: %d\n",
            nrow(all_results), ncol(all_results)))

# ── Save merged file ───────────────────────────────────────
out_file <- file.path(output_dir, "trio_posER_BRCA_results_ALL.txt")
fwrite(all_results, file = out_file, sep = "\t")

cat(sprintf("Saved: %s\n", out_file))


# ==========================================================
## Merge info:
# # 08/07/26 Fri: Ondemand done with
# array 10: 1-58
# array 1/2/6/7:1-56
# array 5/8: 1-54
# array 4: 1-52
# array 3 : 1-55
# ==== Merge above parts:
## 496 shard files:
## Total rows: 245899 | Total cols: 24
# ==========================================================
