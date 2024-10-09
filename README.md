# Yux

* 编译流程

```cmd
cd CBSmode
mkdir build
cd build
cmake -DENABLE_TEST=ON ..
make -j$(nproc)
cd Yux
./homoYux2
```
