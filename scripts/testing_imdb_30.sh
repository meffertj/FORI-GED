for g in 666
do
for h in 272 775 806 900 1290 1426
do
./gurobi_imdb -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/IMDB/IMDB_30
./gurobi_imdb -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/IMDB/IMDB_30
./gurobi_imdb -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/IMDB/IMDB_30
done
done
