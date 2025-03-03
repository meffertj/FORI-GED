gi=0
for g in 4216 1209 373 3437 2938 922 3018 1297 2729 3932
do
hi=0
for h in 4216 1209 373 3437 2938 922 3018 1297 2729 3932
do
if [ $hi -gt $gi ]
then
./gurobi_muta -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Muta/Muta_31_40
./gurobi_muta -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Muta/Muta_31_40
./gurobi_muta -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Muta/Muta_31_40
./gurobi_muta -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Muta/Muta_31_40
fi
((hi++))
done
((gi++))
done
