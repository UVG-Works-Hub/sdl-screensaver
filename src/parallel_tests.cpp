/*
    Julia-Mandelbrot Screensaver Program with Test Mode
    File Name: parallel_tests.cpp
    Authors: Jose Gómez, Adrian Rodriguez, Samuel Chamalé
    Date: September 3, 2024

    Description:
    This program renders a Julia-Mandelbrot fractal in real-time using SDL2 for graphics and
    OpenMP for parallel computation. The user can zoom in and out, move around the fractal,
    and adjust several parameters such as maximum iterations, speed, and cycle length.
    It includes an optimized SIMD-based Mandelbrot set calculation for improved performance
    on modern hardware.

    In addition to the real-time rendering, it includes a test mode that runs a series of
    performance tests at different zoom levels and iteration counts, outputting results to a CSV file.

    Features:
    - Zoom in/out and move using keyboard controls.
    - Adjustable maximum iterations and threading for better performance.
    - Real-time rendering with FPS display.
    - Pausing and unpausing the animation.
    - Test mode for benchmarking performance, results saved to CSV file.

    Build Dependencies:
    - SDL2 library
    - OpenMP library
    - AVX instructions for SIMD optimization (x86-64 architecture)

    Controls:
    - Arrow keys: Move the fractal center
    - Z/X: Zoom in/out
    - I/K: Increase/Decrease max iterations
    - C/V: Increase/Decrease cycle length
    - Space: Pause/unpause animation
    - Number keys (1-8): Set number of threads
    - ,/. : Decrease/Increase time speed

    Test Mode:
    - Run with the `-t` or `--test` flag to activate performance tests.
    - Results saved to "parallel_test_results.csv".

    Compilation:
    make

    License:
    MIT License
*/

#include <SDL2/SDL.h>
#include <cmath>
#include <string>
#include <chrono>
#include <omp.h>
#include <memory>
#include <immintrin.h>
#include <random>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

// Structure to hold color information
struct Color {
    Uint8 r, g, b, a;
};

// Structure to hold pixel information
struct Pixel {
    Uint8 r, g, b, a;
};

// Class to manage the image buffer
class ImageBuffer {
private:
    std::unique_ptr<Pixel[]> buffer;
    int width, height;

public:
    ImageBuffer(int w, int h) : width(w), height(h) {
        buffer = std::make_unique<Pixel[]>(w * h);
    }

    Pixel& at(int x, int y) {
        return buffer[y * width + x];
    }

    const Pixel& at(int x, int y) const {
        return buffer[y * width + x];
    }

    int getWidth() const { return width; }
    int getHeight() const { return height; }

    const void* getBufferData() const {
        return buffer.get();
    }
};

// Function to map a value from one range to another using double precision
double mapValue(double value, double fromLow, double fromHigh, double toLow, double toHigh) {
    return (value - fromLow) * (toHigh - toLow) / (fromHigh - fromLow) + toLow;
}

// SIMD optimized function for Julia-Mandelbrot set calculation
inline __m256d julia_mandelbrot_simd(__m256d x0, __m256d y0, __m256d cr, __m256d ci, int maxIterations) {
    __m256d x = x0;
    __m256d y = y0;
    __m256d x2 = _mm256_mul_pd(x, x);
    __m256d y2 = _mm256_mul_pd(y, y);
    __m256d iterations = _mm256_setzero_pd();
    __m256d four = _mm256_set1_pd(4.0);
    __m256d one = _mm256_set1_pd(1.0);
    __m256d mask = _mm256_setzero_pd();

    for (int i = 0; i < maxIterations; ++i) {
        mask = _mm256_cmp_pd(_mm256_add_pd(x2, y2), four, _CMP_LE_OS);
        if (_mm256_movemask_pd(mask) == 0)
            break;

        y = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(2.0), x), y), ci);
        x = _mm256_add_pd(_mm256_sub_pd(x2, y2), cr);
        x2 = _mm256_mul_pd(x, x);
        y2 = _mm256_mul_pd(y, y);
        iterations = _mm256_add_pd(iterations, _mm256_and_pd(one, mask));
    }

    return iterations;
}

// Function to calculate complex constant for Julia set
std::pair<double, double> calculateComplexConstant(double t, double cycleLength) {
    t = fmod(t, cycleLength);
    double angle = (t / cycleLength) * 2 * M_PI;
    double cr = 0.7885 * std::cos(angle);
    double ci = 0.7885 * std::sin(angle);
    return {cr, ci};
}

// Function to get color based on iteration count
SDL_Color getColor(double iteration, int max_iteration)
{
    if (iteration >= max_iteration)
        return {0, 0, 0, 255}; // Black for points inside the set

    // Smooth coloring
    double smoothed = iteration + 1 - log(log(2.0)) / log(2.0);
    smoothed = smoothed / max_iteration;

    // Create an electric blue to white gradient
    double t = pow(smoothed, 0.5);
    Uint8 r = static_cast<Uint8>(9 * (1 - t) * t * t * t * 255);
    Uint8 g = static_cast<Uint8>(15 * (1 - t) * (1 - t) * t * t * 255);
    Uint8 b = static_cast<Uint8>(8.5 * (1 - t) * (1 - t) * (1 - t) * t * 255);

    return {r, g, b, 255};
}

// Function to calculate elements generated
int calculateElementsGenerated(int screenWidth, int screenHeight, double zoom) {
    return static_cast<int>(screenWidth * screenHeight * zoom * zoom / 16.0);
}

// Main rendering function
void renderJuliaMandelbrot(SDL_Renderer *renderer, ImageBuffer &imageBuffer,
                           double centerX, double centerY, double zoom, int maxIterations,
                           int screenWidth, int screenHeight, double time, double cycleLength) {
    double aspectRatio = static_cast<double>(screenWidth) / screenHeight;
    double rangeY = 4.0 / zoom;
    double rangeX = rangeY * aspectRatio;

    double minReal = centerX - rangeX / 2;
    double maxReal = centerX + rangeX / 2;
    double minImaginary = centerY - rangeY / 2;
    double maxImaginary = centerY + rangeY / 2;

    auto [cr, ci] = calculateComplexConstant(time, cycleLength);
    __m256d cr_vec = _mm256_set1_pd(cr);
    __m256d ci_vec = _mm256_set1_pd(ci);

    #pragma omp parallel for collapse(2) schedule(dynamic, 1)
    for (int py = 0; py < screenHeight; py++) {
        for (int px = 0; px < screenWidth; px += 4) {
            double y0 = mapValue(py, 0, screenHeight - 1, minImaginary, maxImaginary);

            __m256d x0 = _mm256_setr_pd(
                mapValue(px, 0, screenWidth - 1, minReal, maxReal),
                mapValue(px + 1, 0, screenWidth - 1, minReal, maxReal),
                mapValue(px + 2, 0, screenWidth - 1, minReal, maxReal),
                mapValue(px + 3, 0, screenWidth - 1, minReal, maxReal)
            );
            __m256d y0_vec = _mm256_set1_pd(y0);

            __m256d iterations = julia_mandelbrot_simd(x0, y0_vec, cr_vec, ci_vec, maxIterations);

            double iter_array[4];
            _mm256_storeu_pd(iter_array, iterations);

            #pragma omp simd
            for (int i = 0; i < 4; ++i) {
                double iteration = iter_array[i];
                if (iteration < maxIterations) {
                    iteration = iteration + 1 - log(log(2.0)) / log(2.0);
                }

                SDL_Color color = getColor(iteration, maxIterations);

                if ((px + i) < screenWidth) {
                    imageBuffer.at(px + i, py) = {color.r, color.g, color.b, color.a};
                }
            }
        }
    }

    SDL_Texture *texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, imageBuffer.getWidth(), imageBuffer.getHeight());
    SDL_UpdateTexture(texture, NULL, imageBuffer.getBufferData(), imageBuffer.getWidth() * sizeof(Pixel));
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_DestroyTexture(texture);
}

// New structure to hold test results
struct TestResult {
    double zoom;
    int maxIterations;
    double averageFPS;
    double averageRenderTime;
    int elementsGenerated;
};

// Function to run tests and save results
void runTests(SDL_Renderer *renderer, ImageBuffer &imageBuffer, int screenWidth, int screenHeight,
              double centerX, double centerY, double cycleLength, int numThreads) {
    std::vector<TestResult> results;
    std::vector<double> zoomLevels = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0};
    std::vector<int> iterationLevels = {100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200};

    for (double zoom : zoomLevels) {
        for (int maxIterations : iterationLevels) {
            TestResult result;
            result.zoom = zoom;
            result.maxIterations = maxIterations;

            double totalFPS = 0.0;
            double totalRenderTime = 0.0;
            int measurementCount = 7;

            while (measurementCount < 17) {
                auto frameStart = std::chrono::high_resolution_clock::now();

                double time = static_cast<double>(measurementCount);
                renderJuliaMandelbrot(renderer, imageBuffer, centerX, centerY, zoom, maxIterations,
                                      screenWidth, screenHeight, time, cycleLength);
                SDL_RenderPresent(renderer);

                auto frameEnd = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> frameDuration = frameEnd - frameStart;

                double fps = 1.0 / frameDuration.count();
                totalFPS += fps;
                totalRenderTime += frameDuration.count();
                measurementCount++;

                if (fps < 30.0 && measurementCount >= 3) {
                    break;
                }
            }

            result.averageFPS = totalFPS / measurementCount;
            result.averageRenderTime = totalRenderTime / measurementCount;
            result.elementsGenerated = calculateElementsGenerated(screenWidth, screenHeight, zoom);

            results.push_back(result);

            std::cout << "Completed test: Zoom " << zoom << ", Max Iterations " << maxIterations << std::endl;
        }
    }

    std::ofstream outFile("parallel_test_results.csv");
    outFile << "Zoom,MaxIterations,AverageFPS,AverageRenderTime,ElementsGenerated,ScreenWidth,ScreenHeight,Threads\n";
    for (const auto &result : results) {
        outFile << std::fixed << std::setprecision(2)
                << result.zoom << ","
                << result.maxIterations << ","
                << result.averageFPS << ","
                << result.averageRenderTime << ","
                << result.elementsGenerated << ","
                << screenWidth << ","
                << screenHeight << ","
                << numThreads << "\n";
    }
    outFile.close();

    std::cout << "Test results saved to parallel_test_results.csv" << std::endl;
}

int main() {
    // Default parameters
    int screenWidth = 800;
    int screenHeight = 600;
    bool fullscreen = false;
    bool borderless = false;
    double centerX = -0.3;
    double centerY = 0.0;
    double cycleLength = 30.0;
    int numThreads = omp_get_max_threads();
    bool runTestMode = true;

    // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "SDL could not initialize! SDL_Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    Uint32 windowFlags = SDL_WINDOW_RESIZABLE | SDL_WINDOW_SHOWN;
    if (fullscreen) {
        windowFlags |= SDL_WINDOW_FULLSCREEN;
    } else if (borderless) {
        windowFlags |= SDL_WINDOW_BORDERLESS;
    }

    SDL_Window* window = SDL_CreateWindow("Julia-Mandelbrot Screensaver", SDL_WINDOWPOS_UNDEFINED,
                                          SDL_WINDOWPOS_UNDEFINED, screenWidth,
                                          screenHeight, windowFlags);
    if (window == nullptr) {
        std::cerr << "Window could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == nullptr) {
        std::cerr << "Renderer could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    omp_set_num_threads(numThreads);
    ImageBuffer imageBuffer(screenWidth, screenHeight);

    if (runTestMode) {
        runTests(renderer, imageBuffer, screenWidth, screenHeight, centerX, centerY, cycleLength, numThreads);
    } else {
        // Normal rendering and control flow (your original loop code would go here)
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
