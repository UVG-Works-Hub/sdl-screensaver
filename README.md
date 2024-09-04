# Julia-Mandelbrot Screensaver

This project implements a real-time Julia-Mandelbrot fractal screensaver using SDL2 for graphics. It consists of two main versions:

1. Sequential: A basic implementation without parallelization.
2. Parallel: An optimized version using OpenMP for parallel computation and SIMD instructions for improved performance.

The parallel version aims to significantly improve rendering speed, especially on multi-core systems and modern hardware with AVX2 support.

## Features

- Real-time rendering of Julia-Mandelbrot fractals
- Interactive controls for zooming, panning, and adjusting parameters
- Adjustable iteration count, cycle length, and time speed
- FPS display and performance metrics
- Fullscreen and borderless window modes

## Requirements

- C++17 compatible compiler (e.g., GCC, Clang)
- SDL2 library
- OpenMP support (for parallel version)
- AVX2 instruction set support (for SIMD optimizations in parallel version)

## Building

To compile the project, use the provided Makefile:

```bash
make all
```

This will create four executables in the `bin` directory:

- `sequential`: The basic sequential version
- `parallel`: The optimized parallel version
- `sequential_tests`: Tests for the sequential version
- `parallel_tests`: Tests for the parallel version

You can also build individual executables:

```bash
make sequential
make parallel
make sequential_tests
make parallel_tests
```

## Running

To run the sequential version:

```bash
./bin/sequential [options]
```

To run the parallel version:

```bash
./bin/parallel [options]
```

Available options:

- `-w, --width`: Set window width (default: 800)
- `-h, --height`: Set window height (default: 600)
- `-x, --centerX`: Set initial X center (default: -0.3)
- `-y, --centerY`: Set initial Y center (default: 0.0)
- `-z, --zoom`: Set initial zoom level (default: 2.0)
- `-zs, --zoomSpeed`: Set zoom speed (default: 1.01)
- `-ms, --moveSpeed`: Set movement speed (default: 0.1)
- `-mi, --maxIterations`: Set maximum iterations (default: 1000)
- `-ts, --timeSpeed`: Set time speed (default: 0.003)
- `-cl, --cycleLength`: Set cycle length (default: 30.0)
- `-t, --threads`: Set number of threads (parallel version only, default: max available)
- `-fs, --fullscreen`: Enable fullscreen mode (0 or 1)
- `-bw, --borderless`: Enable borderless window mode (0 or 1)

## Performance Comparison

To run performance comparisons between the sequential and parallel versions, use the provided R script:

```bash
Rscript compare.r
```

This will generate comparison graphs and other metrics showing the performance differences between the two implementations.

## Controls

- Arrow keys: Move the fractal center
- Z/X: Zoom in/out
- I/K: Increase/Decrease max iterations
- C/V: Increase/Decrease cycle length
- Space: Pause/unpause animation
- Number keys (1-8): Set number of threads (parallel version only)
- ,/. : Decrease/Increase time speed
- F11: Toggle fullscreen mode
- F12: Toggle borderless window mode

## Gallery

![image](https://github.com/user-attachments/assets/d3538d50-331e-45ee-9593-46ae451dd665)
![image](https://github.com/user-attachments/assets/ec462740-4ddb-4a0c-9533-ac12a76c8c86)
![image](https://github.com/user-attachments/assets/e5d6b396-0695-445e-a9f8-aaafef559531)
![image](https://github.com/user-attachments/assets/084cd6a2-9dda-434c-9189-2ab6e871e9ca)

## License

This project is licensed under the MIT License.
