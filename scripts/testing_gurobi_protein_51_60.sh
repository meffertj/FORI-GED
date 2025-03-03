gi=0
for g in 578 591 310 112 292 465 594 16 597 306
do
hi=0
for h in 578 591 310 112 292 465 594 16 597 306
do
if [ $hi -gt $gi ]
then
./gurobi_protein -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Protein/Protein_51_60
./gurobi_protein -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Protein/Protein_51_60
./gurobi_protein -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Protein/Protein_51_60
./gurobi_protein -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Protein/Protein_51_60
fi
((hi++))
done
((gi++))
done
