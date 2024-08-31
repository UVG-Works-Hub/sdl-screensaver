#include <SDL2/SDL.h>
#include <cmath>
#include <string>
#include <chrono>
#include <omp.h>
#include <vector>
#include <memory>
#include <immintrin.h>

SDL_Color getColor(float iteration, int max_iteration) {
    if (iteration >= max_iteration) return {0, 0, 0, 255};  // Black for points inside the set

    // Smooth coloring
    float smoothed = iteration + 1 - log(log(2.0f)) / log(2.0f);
    smoothed = smoothed / max_iteration;

    // Create an electric blue to white gradient
    float t = pow(smoothed, 0.5);
    Uint8 r = static_cast<Uint8>(9 * (1-t) * t * t * t * 255);
    Uint8 g = static_cast<Uint8>(15 * (1-t) * (1-t) * t * t * 255);
    Uint8 b = static_cast<Uint8>(8.5 * (1-t) * (1-t) * (1-t) * t * 255);

    return {r, g, b, 255};
}

struct Pixel {
    Uint8 r, g, b, a;
};

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

    // Add this new method to get a pointer to the buffer data
    const void* getBufferData() const {
        return buffer.get();
    }
};


float mapValue(float value, float fromLow, float fromHigh, float toLow, float toHigh) {
    return (value - fromLow) * (toHigh - toLow) / (fromHigh - fromLow) + toLow;
}

inline __m256 mandelbrot_simd(__m256 x0, __m256 y0, int maxIterations) {
    __m256 x = _mm256_setzero_ps();
    __m256 y = _mm256_setzero_ps();
    __m256 x2 = _mm256_setzero_ps();
    __m256 y2 = _mm256_setzero_ps();
    __m256 iterations = _mm256_setzero_ps();
    __m256 four = _mm256_set1_ps(4.0f);
    __m256 two = _mm256_set1_ps(2.0f);  // Optimization: Pre-load 2.0f to reduce memory access
    __m256 one = _mm256_set1_ps(1.0f);
    // Optimization: Initialize mask with all bits set to 1
    __m256 mask = _mm256_castsi256_ps(_mm256_set1_epi32(-1));

    for (int i = 0; i < maxIterations; ++i) {
        // Optimization: Simplified calculations
        __m256 xy = _mm256_mul_ps(x, y);
        __m256 x_new = _mm256_add_ps(_mm256_sub_ps(x2, y2), x0);
        y = _mm256_add_ps(_mm256_mul_ps(two, xy), y0);
        x = x_new;

        x2 = _mm256_mul_ps(x, x);
        y2 = _mm256_mul_ps(y, y);

        // Optimization: Moved mask update after x2 and y2 calculations
        __m256 cmp = _mm256_cmp_ps(_mm256_add_ps(x2, y2), four, _CMP_LE_OS);
        mask = _mm256_and_ps(mask, cmp);
        iterations = _mm256_add_ps(iterations, _mm256_and_ps(one, mask));

        // Optimization: Early exit check moved to end of loop
        if (_mm256_movemask_ps(mask) == 0) break;
    }

    return iterations;
}

void renderMandelbrot(SDL_Renderer* renderer, ImageBuffer& imageBuffer,
                      float centerX, float centerY, float zoom, int maxIterations, int SCREEN_WIDTH, int SCREEN_HEIGHT) {
    // Pre-calculate constants to reduce computations in the inner loops
    float aspectRatio = static_cast<float>(SCREEN_WIDTH) / SCREEN_HEIGHT;
    float rangeY = 4.0f / zoom;
    float rangeX = rangeY * aspectRatio;
    float minReal = centerX - rangeX / 2;
    float maxReal = centerX + rangeX / 2;
    float minImaginary = centerY - rangeY / 2;
    float maxImaginary = centerY + rangeY / 2;

    // SIMD vector constants
    __m256 v_minReal = _mm256_set1_ps(minReal);
    __m256 v_factor = _mm256_set1_ps((maxReal - minReal) / (SCREEN_WIDTH - 1));
    __m256 v_centerX = _mm256_set1_ps(centerX);
    __m256 v_centerY = _mm256_set1_ps(centerY);
    __m256 v_rangeX_half = _mm256_set1_ps(rangeX / 2);
    __m256 v_four = _mm256_set1_ps(4.0f);
    __m256 v_maxIterations = _mm256_set1_ps(static_cast<float>(maxIterations));

    // Optimization: Use a larger chunk size for better cache utilization
    const int CHUNK_SIZE = 32;

    // Parallel processing of rows
    #pragma omp parallel for schedule(dynamic, 1)
    for (int py = 0; py < SCREEN_HEIGHT; py++) {
        float y0 = mapValue(py, 0, SCREEN_HEIGHT - 1, minImaginary, maxImaginary);
        __m256 v_y0 = _mm256_set1_ps(y0);

        // Process pixels in chunks for better vectorization
        for (int px = 0; px < SCREEN_WIDTH; px += CHUNK_SIZE) {
            __m256 iterations[CHUNK_SIZE / 8];
            __m256 x[CHUNK_SIZE / 8];

            // Calculate x0 for the chunk
            for (int i = 0; i < CHUNK_SIZE / 8; ++i) {
                __m256 v_px = _mm256_set_ps(px + i*8 + 7, px + i*8 + 6, px + i*8 + 5, px + i*8 + 4,
                                            px + i*8 + 3, px + i*8 + 2, px + i*8 + 1, px + i*8);
                x[i] = _mm256_add_ps(v_minReal, _mm256_mul_ps(v_px, v_factor));
            }

            // Compute Mandelbrot iterations for the chunk
            for (int i = 0; i < CHUNK_SIZE / 8; ++i) {
                iterations[i] = mandelbrot_simd(x[i], v_y0, maxIterations);
            }

            // Process the results for the chunk
            for (int i = 0; i < CHUNK_SIZE && px + i < SCREEN_WIDTH; ++i) {
                float iteration = reinterpret_cast<float*>(&iterations[i / 8])[i % 8];

                // Smooth coloring
                if (iteration < maxIterations) {
                    iteration = iteration + 1 - log2f(log2f(x[i / 8][i % 8] * x[i / 8][i % 8] + y0 * y0));
                }

                SDL_Color color = getColor(iteration, maxIterations);

                // Optimized gradient calculation
                float dx = x[i / 8][i % 8] - centerX;
                float dy = y0 - centerY;
                float distanceSquared = dx * dx + dy * dy;
                float gradientFactor = 1.0f - std::min(sqrtf(distanceSquared) / (rangeX / 2), 1.0f);

                // Apply gradient effect
                color.r = static_cast<Uint8>(color.r * gradientFactor);
                color.g = static_cast<Uint8>(color.g * gradientFactor);
                color.b = static_cast<Uint8>(std::min(255.0f, color.b * gradientFactor + 40));

                imageBuffer.at(px + i, py) = {color.r, color.g, color.b, color.a};
            }
        }
    }

    // Optimization: Use streaming texture update for faster rendering
    SDL_Texture* texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, SCREEN_WIDTH, SCREEN_HEIGHT);
    SDL_UpdateTexture(texture, NULL, imageBuffer.getBufferData(), SCREEN_WIDTH * sizeof(Pixel));
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_DestroyTexture(texture);
}

int main() {
    // Initialize SDL
    SDL_Init(SDL_INIT_VIDEO);

    const int SCREEN_WIDTH = 800;
    const int SCREEN_HEIGHT = 600;

    SDL_Window* window = SDL_CreateWindow("Mandelbrot Zoom", SDL_WINDOWPOS_UNDEFINED,
                                          SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH,
                                          SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    // Mandelbrot set parameters
    float centerX = -0.67f;
    float centerY = 0.0f;
    float zoom = 1.74f;
    float zoomSpeed = 1.01f;
    float moveSpeed = 0.1f;
    int maxIterations = 1000;

    // OpenMP setup
    int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    // Create multiple image buffers for double buffering
    const int NUM_BUFFERS = 2;
    std::vector<ImageBuffer> imageBuffers;
    for (int i = 0; i < NUM_BUFFERS; ++i) {
        imageBuffers.emplace_back(SCREEN_WIDTH, SCREEN_HEIGHT);
    }

    // Variables for FPS calculation
    auto lastTime = std::chrono::high_resolution_clock::now();
    int frameCount = 0;
    float fps = 0.0f;

    // Main loop variables
    bool quit = false;
    SDL_Event e;
    int currentBuffer = 0;
    bool bufferReady = false;

    // Optimization: Pre-allocate textures
    std::vector<SDL_Texture*> textures;
    for (int i = 0; i < NUM_BUFFERS; ++i) {
        textures.push_back(SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, SCREEN_WIDTH, SCREEN_HEIGHT));
    }

    // Optimization: Use a separate thread for rendering
    std::thread renderThread([&]() {
        while (!quit) {
            int bufferToRender = (currentBuffer + 1) % NUM_BUFFERS;
            renderMandelbrot(renderer, imageBufffers[bufferToRender], centerX, centerY, zoom, maxIterations, SCREEN_WIDTH, SCREEN_HEIGHT);
            
            #pragma omp atomic write
            bufferReady = true;
        }
    });

    while (!quit) {
        // Handle events
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) {
                quit = true;
            } else if (e.type == SDL_KEYDOWN) {
                // Optimization: Use OpenMP tasks for input handling
                #pragma omp task
                {
                    switch (e.key.keysym.sym) {
                        case SDLK_UP:    centerY -= moveSpeed / zoom; break;
                        case SDLK_DOWN:  centerY += moveSpeed / zoom; break;
                        case SDLK_LEFT:  centerX -= moveSpeed / zoom; break;
                        case SDLK_RIGHT: centerX += moveSpeed / zoom; break;
                        case SDLK_z:     zoom *= zoomSpeed; break;
                        case SDLK_x:     zoom /= zoomSpeed; break;
                        case SDLK_i:     maxIterations += 100; break;
                        case SDLK_k:     maxIterations = std::max(100, maxIterations - 100); break;
                    }
                }
            }
        }

        // Check if a new buffer is ready
        bool isReady = false;
        #pragma omp atomic read
        isReady = bufferReady;

        if (isReady) {
            // Update texture with the new buffer
            SDL_UpdateTexture(textures[currentBuffer], NULL, imageBuffers[currentBuffer].getBufferData(), SCREEN_WIDTH * sizeof(Pixel));

            // Clear screen and render
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderClear(renderer);
            SDL_RenderCopy(renderer, textures[currentBuffer], NULL, NULL);
            SDL_RenderPresent(renderer);

            // Switch buffers
            currentBuffer = (currentBuffer + 1) % NUM_BUFFERS;
            
            #pragma omp atomic write
            bufferReady = false;

            // Calculate FPS
            frameCount++;
            auto currentTime = std::chrono::high_resolution_clock::now();
            std::chrono::duration<float> elapsedTime = currentTime - lastTime;

            if (elapsedTime.count() >= 1.0f) {
                fps = frameCount / elapsedTime.count();
                frameCount = 0;
                lastTime = currentTime;

                // Update window title with stats
                std::string windowTitle = "fps: " + std::to_string(fps) +
                                          " - z: " + std::to_string(zoom) +
                                          " - y: " + std::to_string(centerY) +
                                          " - x: " + std::to_string(centerX) +
                                          " - mi: " + std::to_string(maxIterations) +
                                          " - threads: " + std::to_string(numThreads);
                SDL_SetWindowTitle(window, windowTitle.c_str());
            }
        }

        // Optimization: Yield to other threads if no new frame is ready
        if (!isReady) {
            std::this_thread::yield();
        }
    }

    // Clean up
    renderThread.join();
    for (auto texture : textures) {
        SDL_DestroyTexture(texture);
    }
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}