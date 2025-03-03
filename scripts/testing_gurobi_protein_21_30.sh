gi=0
for g in 416 358 446 193 586 451 248 202 375 51
do
hi=0
for h in 416 358 446 193 586 451 248 202 375 51
do
if [ $hi -gt $gi ]
then
./gurobi_protein -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Protein/Protein_21_30
./gurobi_protein -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Protein/Protein_21_30
./gurobi_protein -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Protein/Protein_21_30
./gurobi_protein -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Protein/Protein_21_30
fi
((hi++))
done
((gi++))
done
