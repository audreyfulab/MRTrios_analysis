library(data.table)

## Ondemand path:
input_dir <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_BLCA_part"
output_dir <- "/wsu/home/hb/hb68/hb6890/fulab/MRTrios/Output_BLCA"

## Local path:
input_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BLCA/Output_BLCA_part"
output_dir <- "/Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BLCA/Output_BLCA"

# ── Find all shard files (exclude merged ALL file) ─────────
shard_files <- sort(list.files(input_dir,
                               pattern = "trio_BLCA_results_part.*\\.txt",
                               full.names = TRUE))

cat(sprintf("Found %d shard files. Merging...\n", length(shard_files)))

# ── Merge all shards ───────────────────────────────────────
all_results <- rbindlist(lapply(shard_files, fread), fill = TRUE)

cat(sprintf("Total rows: %d | Total cols: %d\n",
            nrow(all_results), ncol(all_results)))

# ── Save merged file ───────────────────────────────────────
out_file <- file.path(output_dir, "trio_results_BLCA_ALL.txt")
fwrite(all_results, file = out_file, sep = "\t")

cat(sprintf("Saved: %s\n", out_file))


# ── OUPUT Info ───────────────────────────────────────

## Found 620 shard files. Merging...
## Total rows: 308405 | Total cols: 24
# head(all_results)
#   trio_row meth.row cna.row gene.row Total.PC.Count               Confounders   b11   b12   b21
#      <int>    <int>   <int>    <int>          <int>                    <char> <int> <int> <int>
#1:        1        1   18353    14709              9 age|sex|race|U_exp|U_meth     1     0     1
#2:        2      238   18353    14709             12 age|sex|race|U_exp|U_meth     1     0     0
#3:        3    18643   18353    14709             14 age|sex|race|U_exp|U_meth     1     0     0
#4:        4    31781   18353    14709              9 age|sex|race|U_exp|U_meth     1     0     1
#5:        5    36138   18353    14709             10 age|sex|race|U_exp|U_meth     1     0     1
#6:        6    37917   18353    14709              9 age|sex|race|U_exp|U_meth     1     0     0
 
#  b22 V1:T1 V1:T2         pb11       pb12         pb21      pb22       pV1:T1       pV1:T2
#   <int> <int> <int>        <num>      <num>        <num>     <num>        <num>        <num>
#1:     0     1     1 6.186955e-19 0.96303696 1.546396e-10 0.9570430 1.398269e-11 7.992355e-04
#2:     0     1     0 4.990942e-25 0.58241758 2.857416e-01 0.5594406 1.398269e-11 1.474907e-02
#3:     0     1     0 4.096671e-24 0.26473526 2.852126e-01 0.2407592 1.398269e-11 5.157219e-01
#4:     0     1     1 3.071985e-18 0.09090909 1.926809e-04 0.1238761 1.398269e-11 2.561813e-06
#5:     0     1     1 5.648493e-18 0.26673327 1.459399e-08 0.3026973 1.398269e-11 6.947230e-11
#6:     0     1     0 9.337462e-23 0.92707293 2.783355e-02 0.9030969 1.398269e-11 2.070866e-02

#    Minor.freq    coef11      coef12    coef21      coef22 Inferred.Model
#         <num>     <num>       <num>     <num>       <num>         <char>
#1: 0.004950495  9.397612  0.05182861 -6.588807  0.05182861             M3
#2: 0.004950495 11.133867 -0.60259205  1.069056 -0.60259205           M0.1
#3: 0.004950495 10.886368  1.10727358 -1.070240  1.10727358           M0.1
#4: 0.004950495  9.241631  1.69307280  3.770942  1.69307280             M3
#5: 0.004950495  9.099411 -0.99153860 -5.795145 -0.99153860             M3
#6: 0.004950495 10.496697  0.10492688  2.208264  0.10492688           M0.1

## Saved: /Users/lianzuo/LZ/ResearchProject/Fulab/MRTrios_BLCA/Output_BLCA/trio_results_BLCA_ALL.txt

