for g in 666
do
for h in 83 87 474 511 525 854
do
./gurobi_imdb -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/IMDB/IMDB_10
./gurobi_imdb -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/IMDB/IMDB_10
./gurobi_imdb -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/IMDB/IMDB_10
done
done
