gi=0
for g in 669 186 1068 212 30 3477 1109 2197 3450 662
do
hi=0
for h in 669 186 1068 212 30 3477 1109 2197 3450 662
do
if [ $hi -gt $gi ]
then
./gurobi_muta -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Muta/Muta_61_70
./gurobi_muta -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Muta/Muta_61_70
./gurobi_muta -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Muta/Muta_61_70
./gurobi_muta -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Muta/Muta_61_70
fi
((hi++))
done
((gi++))
done
