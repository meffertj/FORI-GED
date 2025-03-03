gi=0
for g in 849 4089 3911 1927 4007 294 4096 3943 1563 4334
do
hi=0
for h in 849 4089 3911 1927 4007 294 4096 3943 1563 4334
do
if [ $hi -gt $gi ]
then
./gurobi_muta -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Muta/Muta_81_90
./gurobi_muta -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Muta/Muta_81_90
./gurobi_muta -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Muta/Muta_81_90
./gurobi_muta -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Muta/Muta_81_90
fi
((hi++))
done
((gi++))
done
