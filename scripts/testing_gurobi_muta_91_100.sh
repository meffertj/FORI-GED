gi=0
for g in 2490 643 267 292 800 4084 3970 4204 3691 3755
do
hi=0
for h in 2490 643 267 292 800 4084 3970 4204 3691 3755
do
if [ $hi -gt $gi ]
then
./gurobi_muta -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Muta/Muta_91_100
./gurobi_muta -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Muta/Muta_91_100
./gurobi_muta -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Muta/Muta_91_100
./gurobi_muta -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Muta/Muta_91_100
fi
((hi++))
done
((gi++))
done
