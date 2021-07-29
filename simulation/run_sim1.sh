#!/bin/bash

echo "unrealistic simulations (sim 1)"

bash simulation/run_unrealistic.sh 5 10

sleep 4s

bash simulation/run_unrealistic.sh 9 10

sleep 4s

bash simulation/run_unrealistic.sh 10 15

sleep 4s

bash simulation/run_unrealistic.sh 15 20