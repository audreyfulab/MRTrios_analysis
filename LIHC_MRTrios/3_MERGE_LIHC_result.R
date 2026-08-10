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

# head(all_results)

#   trio_row meth.row cna.row gene.row Total.PC.Count               Confounders   b11   b12   b21
#      <int>    <int>   <int>    <int>          <int>                    <char> <int> <int> <int>
#1:        1        1   18353    14709              6 age|sex|race|U_exp|U_meth     1     0     0
#2:        2      238   18353    14709              8 age|sex|race|U_exp|U_meth     1     0     0
#3:        3    18644   18353    14709              8 age|sex|race|U_exp|U_meth     1     0     0
#4:        4    31783   18353    14709              7 age|sex|race|U_exp|U_meth     1     0     1
#5:        5    36140   18353    14709              5 age|sex|race|U_exp|U_meth     1     0     1
#6:        6    37919   18353    14709              8 age|sex|race|U_exp|U_meth     1     0     0

#     b22 V1:T1 V1:T2         pb11      pb12         pb21      pb22       pV1:T1       pV1:T2
#   <int> <int> <int>        <num>     <num>        <num>     <num>        <num>        <num>
#1:     0     1     0 5.636774e-16 0.1908092 3.946543e-02 0.1978022 6.108465e-17 6.691556e-01
#2:     0     1     0 3.764238e-21 0.5434565 1.322606e-01 0.5934066 6.108465e-17 7.740530e-01
#3:     0     1     0 7.466844e-19 0.3786214 1.353608e-01 0.3746254 6.108465e-17 6.883170e-01
#4:     0     1     0 5.058641e-12 0.1438561 1.886645e-03 0.1588412 6.108465e-17 2.481372e-02
#5:     0     1     1 3.118516e-14 0.1348651 9.369886e-10 0.1378621 6.108465e-17 4.396289e-07
#6:     0     1     0 1.817511e-20 0.6443556 8.843450e-01 0.6373626 6.108465e-17 5.805017e-01

#    Minor.freq    coef11     coef12     coef21     coef22 Inferred.Model
#         <num>     <num>      <num>      <num>      <num>         <char>
#1: 0.005540166  8.599503 -1.3338825 -2.0689119 -1.3338825           M0.1
#2: 0.005540166 10.115437  0.5586886  1.5089468  0.5586886           M0.1
#3: 0.005540166  9.423941  0.8779123 -1.4969108  0.8779123           M0.1
#4: 0.005540166  7.282706  1.4384839  3.1435829  1.4384839             M3
#5: 0.005540166  7.957595 -1.5199040 -6.3096812 -1.5199040             M3
#6: 0.005540166  9.910802 -0.4560321 -0.1455746 -0.4560321           M0.1

# table(all_results$Confounders)
 age|sex|race        age|sex|race|U_exp    age|sex|race|U_exp|U_meth       age|sex|race|U_meth
   13                      2416                    302983                     2616

# ────────────────────────────────────────────────────────────────────────────────────────────────────────



