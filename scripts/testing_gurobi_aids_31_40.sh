gi=0
for g in 27783 37270 28302 28795 27416 27414 27934 27304 28863 27782
do
hi=0
for h in 27783 37270 28302 28795 27416 27414 27934 27304 28863 27782
do
if [ $hi -gt $gi ]
then
./gurobi_aids -f FORI -g ${g} -h ${h}  -t 8 -l 600 -w FORI/AIDS/AIDS_31_40
./gurobi_aids -f F1+ -g ${g} -h ${h}  -t 8 -l 600 -w F1Plus/AIDS/AIDS_31_40
./gurobi_aids -f F2 -g ${g} -h ${h}  -t 8 -l 600 -w F2/AIDS/AIDS_31_40
./gurobi_aids -f F2+ -g ${g} -h ${h}  -t 8 -l 600 -w F2Plus/AIDS/AIDS_31_40
fi
((hi++))
done
((gi++))
done
