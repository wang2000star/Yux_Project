#!/bin/bash

# 获取脚本所在的目录
SCRIPT_DIR=$(cd $(dirname $0) && pwd)

# 切换到脚本所在的目录
cd $SCRIPT_DIR

# echo "test_Yus_4开始!"
# ./../build/tests/test_Yus_4
# sleep 10
# echo "test_Yus_6开始!"
# ./../build/tests/test_Yus_6
# sleep 10
# echo "test_HERA_4开始!"
# ./../build/tests/test_HERA_4
# sleep 10
# echo "test_HERA_5开始!"
# ./../build/tests/test_HERA_5
# sleep 10
# echo "test_pasta_4开始!"
# ./../build/tests/test_pasta_4
# sleep 10
# echo "test_pasta2_4开始!"
# ./../build/tests/test_pasta2_4
# sleep 10
# echo "test_pasta_3开始!"
# ./../build/tests/test_pasta_3
# sleep 10
# echo "test_pasta2_3开始!"
# ./../build/tests/test_pasta2_3
# sleep 10

# echo "test_Yus_4_33开始!"
# ./../build/tests/test_Yus_4_33
# sleep 10
# echo "test_Yus_6_33开始!"
# ./../build/tests/test_Yus_6_33
# sleep 10
# echo "test_HERA_4_33开始!"
# ./../build/tests/test_HERA_4_33
# sleep 10
# echo "test_HERA_5_33开始!"
# ./../build/tests/test_HERA_5_33
# sleep 10
# echo "test_pasta_4_33开始!"
# ./../build/tests/test_pasta_4_33
# sleep 10
# echo "test_pasta2_4_33开始!"
# ./../build/tests/test_pasta2_4_33
# sleep 10
# echo "test_pasta_3_33开始!"
# ./../build/tests/test_pasta_3_33
# sleep 10
# echo "test_pasta2_3_33开始!"
# ./../build/tests/test_pasta2_3_33
echo "test_Yus_5开始!"
./../build/tests/test_Yus_5
sleep 10
echo "test_Yus_4_33开始!"
./../build/tests/test_Yus_4_33
sleep 10
echo "test_Yus_5_33开始!"
./../build/tests/test_Yus_5_33
sleep 10
echo "test_Yus_6_33开始!"
./../build/tests/test_Yus_6_33

# echo "test_Yux_14开始!"
# ./../build/tests/test_Yux_14
echo "测试完成!"

# # 定义参数数对
# param_pairs=(
#   "4 0"
#   "4 1"
#   "4 6"
#   "4 7"
#   "5 0"
#   "5 1"
#   "5 2"
#   "5 3"
#   "6 0"
#   "6 1"
#   "6 2"
#   "7 0"
#   "8 0"
# )

# for param_pair in "${param_pairs[@]}"; do
#   set -- $param_pair
#   Nr=$1
#   idx=$2
#   echo "Nr=$Nr, idx=$idx"
#   for i in {1..10}; do
#     echo "第${i}次测试："
#     ./../build/tests/test_Yus_p_C64_ClientAndServer6_${Nr}_${idx} $Nr $idx
#     if [ $i -lt 10 ]; then
#       sleep 10
#     fi
#   done
# done

# python3 calculate.py

## chmod +x run_test_Yus_p_C64_ClientAndServer4.sh
## ./run_test_Yus_p_C64_ClientAndServer4.sh
