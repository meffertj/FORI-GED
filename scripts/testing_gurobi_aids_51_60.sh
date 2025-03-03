gi=0
for g in 15750 2007 42363 1610 41885 1605 26764 16977 24733 2343
do
hi=0
for h in 15750 2007 42363 1610 41885 1605 26764 16977 24733 2343
do
if [ $hi -gt $gi ]
then
./gurobi_aids -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/AIDS/AIDS_51_60
./gurobi_aids -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/AIDS/AIDS_51_60
./gurobi_aids -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/AIDS/AIDS_51_60
./gurobi_aids -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/AIDS/AIDS_51_60
fi
((hi++))
done
((gi++))
done
