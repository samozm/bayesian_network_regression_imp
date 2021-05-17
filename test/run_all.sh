#!/usr/bin/env bash

usage()
{
  echo "Usage: run_all.sh [-h | --help][ -b | --nburn <int>] [ -s | --nsamp <int>]"
}

help_output()
{
  usage
  echo "-b or --nburn allows the input of the number of Gibbs samples to use for burn-in."
  echo "-s or --nsamp allows the input of the number of Gibbs samples to then retain for sampling."
}

nburn=30000
nsamp=20000

nburn_in=`echo "$@" | sed -n -E -e 's/.*-nburn[[:space:]]([0-9]+).*/\1/p' -e 's/.*-b[[:space:]]([0-9]+).*/\1/p'`
nsamp_in=`echo "$@" | sed -n -E -e 's/.*-nsamp[[:space:]]([0-9]+).*/\1/p' -e 's/.*-s[[:space:]]([0-9]+).*/\1/p'`

if [ -z $nburn_in ] || ! [[ $nburn_in =~ ^[0-9]+$ ]]
then
  echo "Invalid or no nburn value given, using 30000"
else
  nburn=$nburn_in
fi

if [ -z $nsamp_in ] || ! [[ $nsamp_in =~ ^[0-9]+$ ]]
then
  echo "Invalid or no nsamp value given, using 20000"
else
  nsamp=$nsamp_in
fi

julia test/simulation1.jl --nburn $nburn --nsamp $nsamp

julia test/simulation2.jl -c 1 --nburn $nburn --nsamp $nsamp

julia test/simulation2.jl -c 2 --nburn $nburn --nsamp $nsamp

open plots/test/gamma_cis_sim1.png
open plots/test/gamma_cis_sim2_case1.png
open plots/test/gamma_cis_sim2_case2.png
