library(data.table)

input_dir <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_LIHC_part"
output_dir <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_LIHC"

# ── Find all shard files (exclude merged ALL file) ─────────
shard_files <- sort(list.files(input_dir,
                               pattern = "trio_results_part.*\\.txt",
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
