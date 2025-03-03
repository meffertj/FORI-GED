gi=0
for g in 17834 27733 30418 21686 30417 27674 11100 17835 23752 16179
do
hi=0
for h in 17834 27733 30418 21686 30417 27674 11100 17835 23752 16179
do
if [ $hi -gt $gi ]
then
./gurobi_aids -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/AIDS/AIDS_71_80
./gurobi_aids -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/AIDS/AIDS_71_80
./gurobi_aids -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/AIDS/AIDS_71_80
./gurobi_aids -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/AIDS/AIDS_71_80
fi
((hi++))
done
((gi++))
done
