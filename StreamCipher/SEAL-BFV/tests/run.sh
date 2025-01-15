#!/bin/bash

# 获取脚本所在的目录
SCRIPT_DIR=$(cd $(dirname $0) && pwd)

# 切换到脚本所在的目录
cd $SCRIPT_DIR

echo "--------------- SEAL-BFV测试开始!!! ---------------"
echo "--------------- test_Yus_4开始! ---------------"
./../build/tests/test_Yus_4
echo "--------------- test_Yus_4完成! ---------------"
sleep 10
echo "--------------- test_Yus_5开始! ---------------"
./../build/tests/test_Yus_5
echo "--------------- test_Yus_5完成! ---------------"
sleep 10
echo "--------------- test_Yus_6开始! ---------------"
./../build/tests/test_Yus_6
echo "--------------- test_Yus_6完成! ---------------"
sleep 10
echo "--------------- test_HERA_4开始! ---------------"
./../build/tests/test_HERA_4
echo "--------------- test_HERA_4完成! ---------------"
sleep 10
echo "--------------- test_HERA_5开始! ---------------"
./../build/tests/test_HERA_5
echo "--------------- test_HERA_5完成! ---------------"
sleep 10
echo "--------------- test_pasta_4开始! ---------------"
./../build/tests/test_pasta_4
echo "--------------- test_pasta_4完成! ---------------"
sleep 10
echo "--------------- test_pasta2_4开始! ---------------"
./../build/tests/test_pasta2_4
echo "--------------- test_pasta2_4完成! ---------------"
# sleep 10
# echo "--------------- test_pasta_3开始! ---------------"
# ./../build/tests/test_pasta_3
# echo "--------------- test_pasta_3完成! ---------------"
# sleep 10
# echo "--------------- test_pasta2_3开始! ---------------"
# ./../build/tests/test_pasta2_3
# echo "--------------- test_pasta2_3完成! ---------------"

echo "--------------- SEAL-BFV测试完成!!! ---------------"

## chmod +x run.sh
## ./run.sh
