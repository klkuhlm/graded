#!/bin/bash

IFS=$'\n' :; PARS=($(cat parameters.in))

NT=$(wc -l tD.in | awk '{print $1}')

#echo "eta ${PARS[0]}" "kappa ${PARS[1]}" "m ${PARS[2]}" "r ${PARS[3]}" "nt ${NT}"

sed -e "s/%%%eta%%%/${PARS[0]}/" \
    -e "s/%%%kappa%%%/${PARS[1]}/" \
    -e "s/%%%m%%%/${PARS[2]}/" \
    -e "s/%%%r%%%/${PARS[3]}/" \
    -e "s/%%%nt%%%/${NT}/" < powerlaw-template.in > powerlaw.in

./powerlaw
