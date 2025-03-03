gi=0
for g in 20646 1215 2948 42157 1555 18844 357 13023 18791 42153
do
hi=0
for h in 20646 1215 2948 42157 1555 18844 357 13023 18791 42153
do
if [ $hi -gt $gi ]
then
./gurobi_aids -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/AIDS/AIDS_41_50
./gurobi_aids -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/AIDS/AIDS_41_50
./gurobi_aids -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/AIDS/AIDS_41_50
./gurobi_aids -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/AIDS/AIDS_41_50
fi
((hi++))
done
((gi++))
done
