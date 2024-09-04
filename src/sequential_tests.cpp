/*
    Julia-Mandelbrot Screensaver Program with Test Mode (Sequential)
    File Name: sequential_tests.cpp
    Authors: Jose Gómez, Adrian Rodriguez, Samuel Chamalé
    Date: September 3, 2024

    Description:
    This program renders a Julia-Mandelbrot fractal in real-time using SDL2 for graphics.
    The user can zoom in and out, move around the fractal,
    and adjust several parameters such as maximum iterations, speed, and cycle length.
    It renders the fractal sequentially without any parallel computation or SIMD optimizations.

    In addition to the real-time rendering, it includes a test mode that runs a series of
    performance tests at different zoom levels and iteration counts, outputting results to a CSV file.

    Features:
    - Zoom in/out and move using keyboard controls.
    - Adjustable maximum iterations.
    - Real-time rendering with FPS display.
    - Pausing and unpausing the animation.
    - Test mode for benchmarking performance, results saved to CSV file.

    Build Dependencies:
    - SDL2 library

    Test Mode:
    - Run with the `-t` or `--test` flag to activate performance tests.
    - Results saved to "sequential_test_results.csv".

    License:
    MIT License
*/

#include <SDL2/SDL.h>
#include <cmath>
#include <string>
#include <chrono>
#include <memory>
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

// Sequential function for Julia-Mandelbrot set calculation
double julia_mandelbrot(double x0, double y0, double cr, double ci, int maxIterations) {
    double x = x0;
    double y = y0;
    int iteration = 0;

    while (x * x + y * y <= 4.0 && iteration < maxIterations) {
        double xtemp = x * x - y * y + cr;
        y = 2.0 * x * y + ci;
        x = xtemp;
        iteration++;
    }

    return iteration;
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

// Main rendering function (Sequential)
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

    for (int py = 0; py < screenHeight; py++) {
        for (int px = 0; px < screenWidth; px++) {
            double y0 = mapValue(py, 0, screenHeight - 1, minImaginary, maxImaginary);
            double x0 = mapValue(px, 0, screenWidth - 1, minReal, maxReal);

            double iteration = julia_mandelbrot(x0, y0, cr, ci, maxIterations);

            SDL_Color color = getColor(iteration, maxIterations);
            imageBuffer.at(px, py) = {color.r, color.g, color.b, color.a};
        }
    }

    SDL_Texture *texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, imageBuffer.getWidth(), imageBuffer.getHeight());
    SDL_UpdateTexture(texture, NULL, imageBuffer.getBufferData(), imageBuffer.getWidth() * sizeof(Pixel));
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_DestroyTexture(texture);
}

// Structure to hold test results
struct TestResult {
    double zoom;
    int maxIterations;
    double averageFPS;
    double averageRenderTime;
    int elementsGenerated;
};

// Function to run tests and save results (Sequential)
void runTests(SDL_Renderer *renderer, ImageBuffer &imageBuffer, int screenWidth, int screenHeight,
              double centerX, double centerY, double cycleLength) {
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

    std::ofstream outFile("sequential_test_results.csv");
    outFile << "Zoom,MaxIterations,AverageFPS,AverageRenderTime,ElementsGenerated,ScreenWidth,ScreenHeight\n";
    for (const auto &result : results) {
        outFile << std::fixed << std::setprecision(2)
                << result.zoom << ","
                << result.maxIterations << ","
                << result.averageFPS << ","
                << result.averageRenderTime << ","
                << result.elementsGenerated << ","
                << screenWidth << ","
                << screenHeight << "\n";
    }
    outFile.close();

    std::cout << "Test results saved to sequential_test_results.csv" << std::endl;
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

    ImageBuffer imageBuffer(screenWidth, screenHeight);

    if (runTestMode) {
        runTests(renderer, imageBuffer, screenWidth, screenHeight, centerX, centerY, cycleLength);
    } else {
        // Normal rendering and control flow (your original loop code would go here)
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
