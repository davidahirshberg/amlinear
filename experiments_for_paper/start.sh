!# /bin/bash

setupvals=(1 2 3 4)
nvals=(500 1000)
pvals=(6 12)
sigmavals=(0.5 1)
sparvals=(3 4)

reps=3

for ((i1=0; i1<${#setupvals[@]} ;i1++))
do
for ((i3=0; i3<${#nvals[@]} ;i3++))
do
for ((i4=0; i4<${#pvals[@]} ;i4++))
do
for ((i5=0; i5<${#sigmavals[@]} ;i5++))
do
for ((i6=0; i6<${#sparvals[@]} ;i6++))
do
for ((i7=0; i7<${#etavals[@]} ;i7++))
do
    setup=${setupvals[$i1]}
    n=${nvals[$i3]}
    p=${pvals[$i4]}
    sigma=${sigmavals[$i5]}
    spar=${sparvals[$i6]}
    eta=${etavals[$i7]}

    fnm="logging/progress-$setup-$n-$p-$sigma-$spar-$eta-$reps.out"
    echo $fnm

    nohup R CMD BATCH --no-save --no-restore "--args $setup $n $p $sigma $spar $reps" run_simu.R $fnm &
done
done
done
done
done
done
