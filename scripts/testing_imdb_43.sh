for g in 666
do
for h in 378 442 685 905 1037 1270
do
./gurobi_imdb -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/IMDB/IMDB_43
./gurobi_imdb -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/IMDB/IMDB_43
./gurobi_imdb -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/IMDB/IMDB_43
done
done
