gi=0
for g in 3529 1751 1727 779 4328 596 2428 4124 2704 2914
do
hi=0
for h in 3529 1751 1727 779 4328 596 2428 4124 2704 2914
do
if [ $hi -gt $gi ]
then
./gurobi_muta -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/Muta/Muta_51_60
./gurobi_muta -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/Muta/Muta_51_60
./gurobi_muta -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/Muta/Muta_51_60
./gurobi_muta -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/Muta/Muta_51_60
fi
((hi++))
done
((gi++))
done
