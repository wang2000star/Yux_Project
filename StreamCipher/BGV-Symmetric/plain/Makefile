CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2 -march=native

# 目标文件名
TARGET = main

# 头文件和源文件
HEADERS = random_bit.hpp tool.hpp Yux2_8.hpp
SOURCES = main.cpp tool.cpp Yux2_8.cpp

# 生成的可执行文件
$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES)

# 清除生成的文件
.PHONY: clean
clean:
	rm -f $(TARGET)