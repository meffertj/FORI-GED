for g in 666
do
for h in 282 290 455 717 1246 1340
do
./gurobi_imdb -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/IMDB/IMDB_20
./gurobi_imdb -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/IMDB/IMDB_20
./gurobi_imdb -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/IMDB/IMDB_20
done
done
