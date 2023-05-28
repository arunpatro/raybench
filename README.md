# raybench
benchmarking cpp-vs-rust for graphics tasks

## build and run
For cpp projects, we use CMake as a build tool. For the `bvh` task:
```sh
cd cpp/bvh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make; ./bvh
```
For rust projects, look at the main() in main.rs, and choose which task to do
```sh
cd rust
cargo run --release
```


