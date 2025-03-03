gi=0
for g in 2933 3903 104 2756 182 3137 1863 2892 3660 3050
do
hi=0
for h in 2933 3903 104 2756 182 3137 1863 2892 3660 3050
do
if [ $hi -gt $gi ]
then
./gurobi_muta -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Muta/Muta_41_50
./gurobi_muta -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Muta/Muta_41_50
./gurobi_muta -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Muta/Muta_41_50
./gurobi_muta -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Muta/Muta_41_50
fi
((hi++))
done
((gi++))
done
