/*
    Julia-Mandelbrot Screensaver Program
    File Name: julia_mandelbrot_screensaver.cpp
    Authors: Jose Daniel Gómez Cabrera,
    Date: September 3, 2024

    Description:
    This program renders a Julia-Mandelbrot fractal in real-time using SDL2 for graphics and
    OpenMP for parallel computation. The user can zoom in and out, move around the fractal,
    and adjust several parameters such as maximum iterations, speed, and cycle length.
    It includes an optimized SIMD-based Mandelbrot set calculation for improved performance
    on modern hardware.

    Features:
    - Zoom in/out and move using keyboard controls.
    - Adjustable maximum iterations and threading for better performance.
    - Real-time rendering with FPS display.
    - Pausing and unpausing the animation.

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
            // Calculate y0 for the current row
            double y0 = mapValue(py, 0, screenHeight - 1, minImaginary, maxImaginary);

            // Prepare x0 values using SIMD for four pixels at a time
            __m256d x0 = _mm256_setr_pd(
                mapValue(px, 0, screenWidth - 1, minReal, maxReal),
                mapValue(px + 1, 0, screenWidth - 1, minReal, maxReal),
                mapValue(px + 2, 0, screenWidth - 1, minReal, maxReal),
                mapValue(px + 3, 0, screenWidth - 1, minReal, maxReal)
            );
            __m256d y0_vec = _mm256_set1_pd(y0);

            // Perform fractal iterations using SIMD
            __m256d iterations = julia_mandelbrot_simd(x0, y0_vec, cr_vec, ci_vec, maxIterations);

            double iter_array[4];
            _mm256_storeu_pd(iter_array, iterations);

            // Process the first 4 pixels in parallel using SIMD
            #pragma omp simd
            for (int i = 0; i < 4; ++i) {
                double iteration = iter_array[i];
                if (iteration < maxIterations) {
                    iteration = iteration + 1 - log(log(2.0)) / log(2.0);
                }

                // Calculate color based on iteration count
                SDL_Color color = getColor(iteration, maxIterations);

                // Apply the gradient effect and store the color in the buffer
                if ((px + i) < screenWidth) { // Handle screen boundary here
                    imageBuffer.at(px + i, py) = {color.r, color.g, color.b, color.a};
                }
            }
        }
    }


    // Create texture and update it with the image buffer data
    SDL_Texture *texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, imageBuffer.getWidth(), imageBuffer.getHeight());
    SDL_UpdateTexture(texture, NULL, imageBuffer.getBufferData(), imageBuffer.getWidth() * sizeof(Pixel));
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_DestroyTexture(texture);
}

int main(int argc, char* argv[]) {
    // Default parameters
    int screenWidth = 800;
    int screenHeight = 600;
    bool fullscreen = false;
    bool borderless = false;
    double centerX = -0.3;
    double centerY = 0.0;
    double zoom = 2.0;
    double zoomSpeed = 1.01;
    double moveSpeed = 0.1;
    int maxIterations = 1000;
    double timeSpeed = 0.003;
    double cycleLength = 30.0;
    int numThreads = omp_get_max_threads();

    // Parse command line arguments
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 < argc) {
            try {
            std::string arg = argv[i];
            if (arg == "-w" || arg == "--width") screenWidth = std::stoi(argv[i+1]);
            else if (arg == "-h" || arg == "--height") screenHeight = std::stoi(argv[i+1]);
            else if (arg == "-x" || arg == "--centerX") centerX = std::stod(argv[i+1]);
            else if (arg == "-y" || arg == "--centerY") centerY = std::stod(argv[i+1]);
            else if (arg == "-z" || arg == "--zoom") zoom = std::stod(argv[i+1]);
            else if (arg == "-zs" || arg == "--zoomSpeed") zoomSpeed = std::stod(argv[i+1]);
            else if (arg == "-ms" || arg == "--moveSpeed") moveSpeed = std::stod(argv[i+1]);
            else if (arg == "-mi" || arg == "--maxIterations") maxIterations = std::stoi(argv[i+1]);
            else if (arg == "-ts" || arg == "--timeSpeed") timeSpeed = std::stod(argv[i+1]);
            else if (arg == "-cl" || arg == "--cycleLength") cycleLength = std::stod(argv[i+1]);
            else if (arg == "-t" || arg == "--threads") numThreads = std::stoi(argv[i+1]);
            else if (arg == "-fs" || arg == "--fullscreen") fullscreen = std::stoi(argv[i+1]) != 0;
            else if (arg == "-bw" || arg == "--borderless") borderless = std::stoi(argv[i+1]) != 0;
            } catch (const std::exception& e) {
                std::cerr << "Error parsing argument: " << argv[i] << std::endl;
                return 1;
            }
        }
    }

    // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "SDL could not initialize! SDL_Error: " << SDL_GetError() << std::endl;
        return 1;
    }


    // Create window
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

    // Create renderer
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == nullptr) {
        std::cerr << "Renderer could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    bool quit = false;
    SDL_Event e;

    auto lastTime = std::chrono::high_resolution_clock::now();
    int frameCount = 0;

    omp_set_num_threads(numThreads);

    ImageBuffer imageBuffer(screenWidth, screenHeight);

    double time = 0.0;
    bool paused = false;

    // Adjust moveSpeed based on zoom, setting a maximum zoom for movement calculation
    double maxMovementZoom = 1e12;  // Define a maximum zoom for the movement scaling
    double effectiveZoom = (zoom > maxMovementZoom) ? maxMovementZoom : zoom;

    while (!quit) {
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) {
                quit = true;
            } else if (e.type == SDL_KEYDOWN) {
                switch (e.key.keysym.sym) {
                    case SDLK_UP:
                        centerY -= moveSpeed / effectiveZoom;
                        break;
                    case SDLK_DOWN:
                        centerY += moveSpeed / effectiveZoom;
                        break;
                    case SDLK_LEFT:
                        centerX -= moveSpeed / effectiveZoom;
                        break;
                    case SDLK_RIGHT:
                        centerX += moveSpeed / effectiveZoom;
                        break;
                    case SDLK_z:
                        zoom *= zoomSpeed;
                        effectiveZoom = (zoom > maxMovementZoom) ? maxMovementZoom : zoom;
                        break;
                    case SDLK_x:
                        zoom /= zoomSpeed;
                        effectiveZoom = (zoom > maxMovementZoom) ? maxMovementZoom : zoom;
                        break;
                    case SDLK_i:
                        maxIterations += 100;
                        break;
                    case SDLK_k:
                        maxIterations = std::max(100, maxIterations - 100);
                        break;
                    case SDLK_c:
                        cycleLength += 1.0;
                        break;
                    case SDLK_v:
                        cycleLength = std::max(1.0, cycleLength - 1.0);
                        break;
                    case SDLK_SPACE:
                        paused = !paused;
                        break;
                    case SDLK_F11:  // Toggle fullscreen
                        fullscreen = !fullscreen;
                        if (fullscreen) {
                            SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);
                        } else {
                            SDL_SetWindowFullscreen(window, 0);
                        }
                        break;
                    case SDLK_F12:  // Toggle borderless windowed mode (windowless)
                        borderless = !borderless;
                        SDL_SetWindowBordered(window, borderless ? SDL_FALSE : SDL_TRUE);
                        break;
                }
            } else if (e.type == SDL_WINDOWEVENT && e.window.event == SDL_WINDOWEVENT_RESIZED) {
                // Handle window resizing
                screenWidth = e.window.data1;
                screenHeight = e.window.data2;

                // Recreate the image buffer with the new size
                imageBuffer = ImageBuffer(screenWidth, screenHeight);
            }
        }

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        auto start = std::chrono::high_resolution_clock::now();
        renderJuliaMandelbrot(renderer, imageBuffer, centerX, centerY, zoom, maxIterations, screenWidth, screenHeight, time, cycleLength);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        SDL_RenderPresent(renderer);

        // Update time if not paused
        if (!paused) {
            time += timeSpeed;
        }

        // Calculate the current position in the cycle (0.0 to 1.0)
        double cyclePosition = fmod(time, cycleLength) / cycleLength;

        // Calculate FPS
        frameCount++;
        auto currentTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedTime = currentTime - lastTime;

        if (elapsedTime.count() >= 1.0) {
            double fps = frameCount / elapsedTime.count();
            frameCount = 0;
            lastTime = currentTime;

            std::string windowTitle = "fps: " + std::to_string(fps) +
                                      " - z: " + std::to_string(zoom) +
                                      " - y: " + std::to_string(centerY) +
                                      " - x: " + std::to_string(centerX) +
                                      " - mi: " + std::to_string(maxIterations) +
                                      " - threads: " + std::to_string(numThreads) +
                                      " - render time: " + std::to_string(elapsed.count()) + "s" +
                                      " - cycle: " + std::to_string(cycleLength) + "s" +
                                      " - pos: " + std::to_string(cyclePosition) +
                                      (paused ? " - PAUSED" : "");
            SDL_SetWindowTitle(window, windowTitle.c_str());
        }

        // Handle thread count changes
        const Uint8 *state = SDL_GetKeyboardState(NULL);
        for (int i = SDL_SCANCODE_1; i <= SDL_SCANCODE_8; ++i) {
            if (state[i]) {
                numThreads = i - SDL_SCANCODE_1 + 1;
                omp_set_num_threads(numThreads);
                break;
            }
        }
        if (state[SDL_SCANCODE_COMMA])
            timeSpeed *= 0.99; // Slow down
        if (state[SDL_SCANCODE_PERIOD])
            timeSpeed *= 1.01; // Speed up
    }

    // Clean up
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}

/*
    References:
        - https://www3.nd.edu/~dthain/courses/cse30341/spring2020/project3/
        - https://scala-blitz.github.io/home/documentation/examples//mandelbrot.html
        - https://en.wikipedia.org/wiki/Julia_set
        - https://youtu.be/uc2yok_pLV4
        - https://complex-analysis.com/contenido/conjunto_de_julia.html.
        - https://www3.gobiernodecanarias.org/medusa/ecoblog/mrodperv/fractales/conjuntos-de-julia-y-mandelbrot/
*/
