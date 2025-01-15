#!/bin/bash

# 获取脚本所在的目录
SCRIPT_DIR=$(cd $(dirname $0) && pwd)

# 切换到脚本所在的目录
cd $SCRIPT_DIR

./HElib-BGV/tests/run.sh

sleep 10

./SEAL-BFV/tests/run.sh

echo "--------------- 所有测试完成!!! ---------------"

## chmod +x run.sh
## ./run.sh
