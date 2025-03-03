gi=0
for g in 3265 1220 2988 1580 13 653 3621 3026 1270 3590
do
hi=0
for h in 3265 1220 2988 1580 13 653 3621 3026 1270 3590
do
if [ $hi -gt $gi ]
then
./gurobi_muta -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Muta/Muta_21_30
./gurobi_muta -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Muta/Muta_21_30
./gurobi_muta -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Muta/Muta_21_30
./gurobi_muta -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Muta/Muta_21_30
fi
((hi++))
done
((gi++))
done
