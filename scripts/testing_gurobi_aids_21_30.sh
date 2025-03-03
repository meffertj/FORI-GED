gi=0
for g in 11059 11108 18434 18350 33065 22710 21992 41803 20076 41804
do
hi=0
for h in 11059 11108 18434 18350 33065 22710 21992 41803 20076 41804
do
if [ $hi -gt $gi ]
then
./gurobi_aids -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/AIDS/AIDS_21_30
./gurobi_aids -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/AIDS/AIDS_21_30
./gurobi_aids -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/AIDS/AIDS_21_30
./gurobi_aids -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/AIDS/AIDS_21_30
fi
((hi++))
done
((gi++))
done
