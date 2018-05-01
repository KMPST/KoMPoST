#bin/trento Pb Pb 10 -o out --grid-max 20.48 --grid-step 0.08 -n 2000 --b-max 1 -p 1
EventID="0"; Ns="512";
 tail -n+9 ${EventID}.dat | perl -pe 's/\s+/\n/g' | awk -v Ns=${Ns} 'BEGIN {for(x=0;x<Ns;x++){for(y=0;y<Ns;y++){T00[x+Ns*y]=0.0;}} Index=0;} {T00[Index]=$1; Index++;} END {for(y=0;y<Ns;y++){for(x=0;x<Ns;x++){print x,y,T00[x+Ns*y],0.5*T00[x+Ns*y],0.5*T00[x+Ns*y],0.0,0.0,0.0,0.0,0.0,0.0,0.0;} printf("\n"); }}' > EnergyMomentumTensorID${ID}.txt
#awk '{sum+=$3};END{print sum}' EnergyMomentumTensorID.txt
