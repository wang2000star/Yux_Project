# Hybrid Sysmmetric Scheme

## Build and Run

### 前期准备

* Ubuntu-22.04
* 依赖项

  ```bash
  sudo apt update
  sudo apt install build-essential
  ```
* GNU make >=3.82

  版本查看：

  ```bash
  make --version
  ```
* **Pthreads**

  什么是Pthreads？[Pthreads概述](https://www.cnblogs.com/blueclue/archive/2010/07/16/1779006.html)版本查看

  ```bash
  getconf GNU_LIBPTHREAD_VERSION
  ```
* git >= 1.83

  版本查看：

```bash
git --version
```

* g++ >= 7.3.1

  ```bash
  # 假设你已经编译并安装了 g++，并且它位于 /usr/local/bin/g++

  # 创建符号链接
  sudo ln -sf /usr/local/bin/g++ /usr/bin/g++

  # 使用 update-alternatives 管理 g++ 版本
  sudo update-alternatives --install /usr/bin/g++ g++ /usr/local/bin/g++ 100
  sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++ 50

  # 选择默认的 g++ 版本
  sudo update-alternatives --config g++

  # 验证设置
  g++ --version
  ```

  版本查看：

  ```bash
  g++ --version
  ```
* cmake >= 3.30

  版本查看：

  ```bash
  cmake --version
  ```

Our implementation relies on [HElib v2.2.1](https://github.com/homenc/HElib/tree/v2.2.1), which is a famous homomomorphic encryption scheme library and `v2.2.1` is now (July, 2022) the newest release version. Therefore, you firstly need to install HElib. The HElib iteratively relies on [NTL](https://libntl.org/), [GMP](https://gmplib.org) and M4 libraries. We suggest that install them to your default user's program path `/usr/local`, and the suggeted version of them are showed in following table. What's more, using `./configure SHARED=on` to get a shared library when compile NTL is also suggested.

| Libraries | Version | Website                                                           |
| --------- | ------- | ----------------------------------------------------------------- |
| M4        | v1.4.18 | [http://mirrors.kernel.org/gnu/m4](http://mirrors.kernel.org/gnu/m4) |
| GMP       | v6.2.1  | [https://gmplib.org](https://gmplib.org)                             |
| NTL       | v11.4.3 | [https://libntl.org](https://libntl.org)                             |
| HElib     | v2.2.1  | [https://github.com/homenc/HElib](https://github.com/homenc/HElib)   |

m4下载安装：

```bash
wget https://mirrors.ustc.edu.cn/gnu/m4/m4-latest.tar.xz
tar -xf m4-latest.tar.xz
rm m4-latest.tar.xz
#根据自己实际文件名填
cd m4-1.4.19
./configure
make
sudo make install
```

m4版本查看：

```bash
hash -r
m4 --version
```

* patchelf：

```bash
wget https://github.com/NixOS/patchelf/releases/download/0.18.0/patchelf-0.18.0.tar.gz
tar -axf patchelf-0.18.0.tar.gz
rm patchelf-0.18.0.tar.gz
cd patchelf-0.18.0
./configure
make
sudo make install
# 验证安装
hash -r
patchelf --version
```

* GMP安装教程

  ```bash
  wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
  tar -xf gmp-6.3.0.tar.xz
  rm gmp-6.3.0.tar.xz
  cd gmp-6.3.0
  ./configure
  make或make -j$(nproc)或make -j16
  make check
  sudo make install
  ```

注意：需要保证gmp.h版本和libgmp版本一致

gmp.h版本查看方法：

```bash
grep "define __GNU_MP_VERSION" /usr/local/include/gmp.h
```

libgmp版本查看方法：

创建check_gmp_version.c文件，编辑以下内容

```c
#include <stdio.h>
#include <gmp.h>

int main() {
    printf("GMP version: %s\n", gmp_version);
    return 0;
}
```

运行

```bash
gcc check_gmp_version.c -lgmp -o check_gmp_version
./check_gmp_version
```

* NTL 安装教程

```bash
  wget https://libntl.org/ntl-11.5.1.tar.gz
  tar -xf ntl-11.5.1.tar.gz
  rm ntl-11.5.1.tar.gz
  cd ntl-11.5.1/src
  ./configure
  make
  make check
  sudo make install
```

NTL版本查看：

```bash
grep "define NTL_VERSION" /usr/local/include/NTL/version.h
```

We suggest that follow the command below to install HElib. You can also refer to the website [HElib GitHub repository](https://github.com/homenc/HElib).

```bash
wget https://github.com/homenc/HElib/archive/refs/tags/v2.3.0.tar.gz
tar -xf v2.3.0.tar.gz
rm v2.3.0.tar.gz
cd HElib-2.3.0
mkdir build
cd build

（方法1）
cmake -DPACKAGE_BUILD=ON -DENABLE_TEST=ON ..
make
ctest
sudo make install

（方法2）
cmake -DPACKAGE_BUILD=ON -DENABLE_TEST=ON -DCMAKE_INSTALL_PREFIX=/home/Alice/helib_install ..
make
make -j$(nproc)
ctest
make install

# 更新动态链接库缓存
sudo ldconfig
export PATH=/home/Alice/helib_install/bin:$PATH
export LD_LIBRARY_PATH=/home/Alice/helib_install/lib:$LD_LIBRARY_PATH

```

错误排查：HElib-2.3.0文件夹里面的VERSION需要编辑为2.3.0

helib_pack/share/cmake/helib里面的：

（1）helibConfig.cmake文件修改为

```cmake
include(CMakeFindDependencyMacro)

# Set the version number
set(helib_VERSION 2.3.0)

# Other configuration settings
# ...

# Include the version file
include("${CMAKE_CURRENT_LIST_DIR}/helibConfigVersion.cmake")
```

```bash
CMake Error at CMakeLists.txt:39 (find_package):
  Could not find a configuration file for package "helib" that exactly
  matches requested version "2.2.0".

  The following configuration files were considered but not accepted:

    /home/wfz/helib_install/helib_pack/share/cmake/helib/helibConfig.cmake, version: unknown
```

（2）helibConfigVersion.cmake文件修改为：

```cmake
# helibConfigVersion.cmake
set(PACKAGE_VERSION "2.3.0")

if(PACKAGE_VERSION VERSION_LESS PACKAGE_FIND_VERSION)
  set(PACKAGE_VERSION_COMPATIBLE FALSE)
else()
  set(PACKAGE_VERSION_COMPATIBLE TRUE)
  if(PACKAGE_FIND_VERSION STREQUAL PACKAGE_VERSION)
    set(PACKAGE_VERSION_EXACT TRUE)
  endif()
endif()
```

本人强烈推荐的一个安装教程：[安装HElib并运行示例程序](https://blog.csdn.net/baishuiniyaonulia/article/details/122737035)

测试：

```bash
cd HElib-2.3.0/examples
mkdir build
cd build
cmake -Dhelib_DIR=/home/Alice/helib_install/helib_pack/share/cmake/helib ..
make
cd ./bin
```

After install [HElib](https://github.com/homenc/HElib), you can run our tests by

```bash
mkdir build
cd build
cmake ..
make
./test/test-trans-Yu2x-8-C16

```

Symmetic test:

```bash
g++ test-Yu2x-8.cpp -Iinclude ../Yux/Yu2x-8.cpp -o test-Yu2x-8
```
