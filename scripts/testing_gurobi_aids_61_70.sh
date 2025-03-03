gi=0
for g in 25279 2009 27771 17832 15757 16981 32612 17831 16979 15752
do
hi=0
for h in 25279 2009 27771 17832 15757 16981 32612 17831 16979 15752
do
if [ $hi -gt $gi ]
then
./gurobi_aids -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/AIDS/AIDS_61_70
./gurobi_aids -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/AIDS/AIDS_61_70
./gurobi_aids -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/AIDS/AIDS_61_70
./gurobi_aids -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/AIDS/AIDS_61_70
fi
((hi++))
done
((gi++))
done
