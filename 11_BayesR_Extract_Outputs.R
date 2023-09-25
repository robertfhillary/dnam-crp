## Run in Terminal 


## Set Directory 
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Outputs/
  ## Set up loop for any output file with a *.csv extension 
  for i in *.csv
do
## Extract columns with sigma name - sigmaG (variance explained) and sigmaE (error term - unexplained proportion of variance)  
sigma1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "sigma" | cut -f 1 |  sed 's/:/\n/g' | awk 'NR==1')
sigma2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "sigma" | cut -f 1 |  sed 's/:/\n/g' | awk 'END{print $NF}')
## Create variable name 
A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d. -f1)
cat $i | cut -d ',' -f $sigma1-$sigma2 > ../Sigma/$A
## Extract columns with beta name - goes from first to last column 
beta1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "beta" | cut -f 1 | sed 's/:/\n/g' | awk 'NR==1')
beta2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "beta" | cut -f 1 | sed 's/:/\n/g' | awk 'END{print $NF}')
A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d_ -f1)
cat $i | cut -d ',' -f $beta1-$beta2 > ../Beta/$A
## Extract columns with components name - goes from first to last column; within each iteration, tells you which group marker was assigned to e.g. small, medium, large effect 
comp1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "comp" | cut -f 1 | sed 's/:/\n/g' | awk 'NR==1')
comp2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "comp" | cut -f 1 | sed 's/:/\n/g' | awk 'END{print $NF}')
A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d_ -f1)
cat $i | cut -d ',' -f $comp1-$comp2 > ../Comp/$A
done
