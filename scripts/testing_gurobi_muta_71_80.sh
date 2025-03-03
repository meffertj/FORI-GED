gi=0
for g in 3998 2908 2273 3485 3329 4266 2366 520 3957 4075
do
hi=0
for h in 3998 2908 2273 3485 3329 4266 2366 520 3957 4075
do
if [ $hi -gt $gi ]
then
./gurobi_muta -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Muta/Muta_71_80
./gurobi_muta -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Muta/Muta_71_80
./gurobi_muta -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Muta/Muta_71_80
./gurobi_muta -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Muta/Muta_71_80
fi
((hi++))
done
((gi++))
done
