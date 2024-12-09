#!/bin/bash

# 获取脚本所在的目录
SCRIPT_DIR=$(cd $(dirname $0) && pwd)

# 切换到脚本所在的目录
cd $SCRIPT_DIR

for i in {1..10}
do
    echo "第${i}次测试："
    ./../build/tests/test_Yus_p_C64_ClientAndServer4
    if [ $i -lt 10 ]; then
        sleep 5
    fi
done

python3 calculate.py

## chmod +x run_test_Yus_p_C64_ClientAndServer4.sh
## ./run_test_Yus_p_C64_ClientAndServer4.sh
