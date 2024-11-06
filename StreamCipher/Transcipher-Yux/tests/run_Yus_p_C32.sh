#!/bin/bash
echo "Starting test_Yus_p_C32_Clientoff..."
./test_Yus_p_C32_Clientoff
echo "test_Yus_p_C32_Clientoff has finished, starting test_Yus_p_C32_Serveroff..."
./test_Yus_p_C32_Serveroff
#echo "test_Yux_p_C16_Serveroff has finished, starting test_Yux_p_C16_Clienton..."
#./test_Yux_p_C16_Clienton
#echo "test_Yux_p_C16_Clienton has finished, starting test_Yux_p_C16_Serveron..."
#./test_Yux_p_C16_Serveron

## chmod +x run_Yus_p_C32.sh
## ./run_Yus_p_C32.sh
