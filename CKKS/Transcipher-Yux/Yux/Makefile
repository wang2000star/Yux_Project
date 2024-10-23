CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2 -march=native

# 目标文件名
TARGET = main

# 头文件和源文件
HEADERS = tool.hpp Yux2_8.hpp Client.hpp Server.hpp params.hpp random_bit.hpp FHEtool.hpp
SOURCES = main.cpp tool.cpp Yux2_8.cpp Client.cpp Server.cpp

# 库路径和库文件
LIB_PATHS = -L./thirdparty/HElib_v2.2.2 -L./thirdparty/ntl-11.5.1
LIBS = -lhelib -lntl

# 包含路径
INCLUDE_PATHS = -I./thirdparty/HElib_v2.2.2/include -I./thirdparty/ntl-11.5.1/include

# 生成的可执行文件
$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATHS) -o $(TARGET) $(SOURCES) $(LIB_PATHS) $(LIBS)

# 清除生成的文件
.PHONY: clean
clean:
	rm -f $(TARGET)