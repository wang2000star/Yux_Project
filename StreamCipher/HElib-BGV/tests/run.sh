#!/bin/bash

# 获取脚本所在的目录
SCRIPT_DIR=$(cd $(dirname $0) && pwd)

# 切换到脚本所在的目录
cd $SCRIPT_DIR

# 定义参数数对
param_pairs=(
  "4 0"
  "4 1"
  "4 6"
  "4 7"
  "5 0"
  "5 1"
  "5 2"
  "5 3"
  "6 0"
  "6 1"
  "6 2"
  "7 0"
  "8 0"
)

for param_pair in "${param_pairs[@]}"; do
  set -- $param_pair
  Nr=$1
  idx=$2
  echo "Nr=$Nr, idx=$idx"
  for i in {1..10}; do
    echo "第${i}次测试："
    ./../build/tests/test_Yus_p_C64_ClientAndServer6_${Nr}_${idx} $Nr $idx
    if [ $i -lt 10 ]; then
      sleep 10
    fi
  done
done

python3 calculate.py

## chmod +x run_test_Yus_p_C64_ClientAndServer4.sh
## ./run_test_Yus_p_C64_ClientAndServer4.sh
