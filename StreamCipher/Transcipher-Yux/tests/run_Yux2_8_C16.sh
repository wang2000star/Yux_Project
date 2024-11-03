#!/bin/bash
./test_Yux2_8_C16_Clientoff
echo "test_Yux2_8_C16_Clientoff has finished, starting test_Yux2_8_C16_Serveroff..."
./test_Yux2_8_C16_Serveroff
echo "test_Yux2_8_C16_Serveroff has finished, starting test_Yux2_8_C16_Clienton..."
./test_Yux2_8_C16_Clienton
echo "test_Yux2_8_C16_Clienton has finished, starting test_Yux2_8_C16_Serveron..."
./test_Yux2_8_C16_Serveron

##chmod +x run_Yux2_8_C16.sh
