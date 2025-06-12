#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

for i in $(seq 0 12);
do
    ${SCRIPT_DIR}/build/ntt_bench_bar "p" $i > "${SCRIPT_DIR}/ntt_q${i}.coe"
    ${SCRIPT_DIR}/build/ntt_bench_bar "n" $i > "${SCRIPT_DIR}/intt_q${i}.coe"
done

