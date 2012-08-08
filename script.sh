#! /bin/bash

#what are the parameters to be estimated 
# alpha and beta
MESSAGE="Genetic algorithm is running"
echo $MESSAGE

echo Please enter values for the parameters alpha and epsilon
read alp eps 

ALPHA=$alp
EPS=$eps
#create initial population of size 
#how the hell do i encode the parameters

N=200
M=8
last=50
sel=0.5
M2=2*$(awk 'BEGIN{print int(0.5*(8/2)+0.5)}')
mutrate=0.01
#nmuts=
a=12.0
b=3.0
c=$(float_eval "$a / $b")

#run the model to generate an initial population (generate a set of solutions)
#pass the initial guesses for the parameters as command-line arguments to my simulator
./bin/Hello-World -p $ALPHA $EPS > o 
matlab -nodesktop -nosplash -r "farzan_zap_fft"

#run for MAX_ITR
for i in `seq 1 $200`; do
    if [ $i = 5 ]; then
	echo "This is 5"
    else
	echo $i
    fi
done

