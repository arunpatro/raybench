## raybench
benchmaring rust-vs-cpp (also mojo and python) for ray tracing

Mentions:
- https://docs.google.com/spreadsheets/d/1OQ_U8fY8DCzwz8CoYbobhHeJ4oKW49x4uOhf4Yy9mwU/edit#gid=0
- https://github.com/niofis/raybench
- https://xania.org/201505/on-rust-performance


### Introduction
The graphics class at NYU was a great motivation for me to dive deep into compiled languages like `C++` and `Rust`. I implement simple raytracing scenes in cpp (class default) and later in rust by trying to replicate similar code structure albeit a few obv differences like `OOP` vs `Struct-Impl-Trait`. This benchmark is purely on CPUs, as GPUs involve Shading Languages and more thoughtful logic. It would be interesting to study building modern day graphics based on `webgpu` , `vulcan`, `metal-3` and and also target the `wasm` architecture. All of them are really optimized! 

We start with a simple explanation of the task in english, followed by python to express it more formally. We show how one might implement this in C++, which is what I had to do in the course. And finally, compare it with a similar rust implementation. I comment on my observations in this section: 

I profile both code to find out the pros/cons of coding in either language. Cpp is the baseline.


### Tasks
We test the languages with 3 tasks:
1. Ray Tracing (3 images)
2. Ray Tracing with BVH accelaration (2 objects)
3. CPU Rasterization (2 videos)


### 1> Ray Tracing
This is a simple explanation of ray tracing in python, that I borrow from https://gist.github.com/omaraflak/08eed161f5390c27fc8fed136f2ff53d. 

Explanation
```
TODO 
```




### Observations
We have some observations:
1. cpp compiles faster and produces a smaller binary and is not significantly faster than rust.
2. rust is easier to debug
3. floating point math is inaccurate and gives different intermediate results for both
4. link-time-optimization (LTO) in rust is faster to compile than rust without LTO ??
5. Rust with LTO is slower than without LTO ??
6. How can thin LTO be faster than LTO ??

### Profiling C++ vs Rust
- I use xcode Instruments to profile both C++ and Rust code. Rust executables have a ordered way of creating debug symbols and this helps in finding the order of the function call inside another. This is not the case with cpp, which just bundles all function calls to one. 
- This can be due to recursion of reflected rays, I can try to profile without reflection. 


### [Image 01]: Reflection with 7 Spheres and 1 Plane
<img src="./scene.png">

```
Objects: 7 Spheres and 1 Plane
Material: Reflective and Specular
Camera: fov: 45, focal length = 5, position = (0, 1, 10), Perspective
```


Benchmarking results:
<table>
  <thead>
    <tr>
      <th></th>
      <th>compile time</th>
      <th>runtime 1600x800</th>
      <th>runtime 2400x1200</th>
      <th>binary size</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>C++</td>
      <td>6.077s</td>
      <td>3.110s</td>
      <td>6.70 +- 0.46s</td>
      <td>109KB</td>
    </tr>
    <tr>
      <td>C++ (G++ LTO)</td>
      <td>-</td>
      <td>-</td>
      <td>4.60 +- 0.49s</td>
      <td>70KB</td>
    </tr>
    <tr>
      <td>Rust</td>
      <td>1m 23s</td>
      <td>1.567s</td>
      <td>2.778s</td>
      <td>1.5MB</td>
    </tr>
    <tr>
      <td>Rust (LTO)</td>
      <td>54.227s</td>
      <td>3.101s</td>
      <td>6.914s</td>
      <td>760KB</td>
    </tr>
    <tr>
      <td>Rust (Thin LTO)</td>
      <td>1m 15s</td>
      <td>2.056s</td>
      <td>4.170s</td>
      <td>1.4MB</td>
    </tr>
  </tbody>
</table>


### Notes
- Rust uses `nlagebra`, CPP uses `Eigen`
- 
