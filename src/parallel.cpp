#include <SDL2/SDL.h>
#include <cmath>
#include <string>
#include <chrono>
#include <omp.h>
#include <vector>
#include <memory>
#include <immintrin.h>

SDL_Color getColor(float iteration, int max_iteration)
{
    if (iteration >= max_iteration)
        return {0, 0, 0, 255}; // Black for points inside the set

    // Smooth coloring
    float smoothed = iteration + 1 - log(log(2.0f)) / log(2.0f);
    smoothed = smoothed / max_iteration;

    // Create an electric blue to white gradient
    float t = pow(smoothed, 0.5);
    Uint8 r = static_cast<Uint8>(9 * (1 - t) * t * t * t * 255);
    Uint8 g = static_cast<Uint8>(15 * (1 - t) * (1 - t) * t * t * 255);
    Uint8 b = static_cast<Uint8>(8.5 * (1 - t) * (1 - t) * (1 - t) * t * 255);

    return {r, g, b, 255};
}

struct Pixel
{
    Uint8 r, g, b, a;
};

class ImageBuffer
{
private:
    std::unique_ptr<Pixel[]> buffer;
    int width, height;

public:
    ImageBuffer(int w, int h) : width(w), height(h)
    {
        buffer = std::make_unique<Pixel[]>(w * h);
    }

    Pixel &at(int x, int y)
    {
        return buffer[y * width + x];
    }

    const Pixel &at(int x, int y) const
    {
        return buffer[y * width + x];
    }

    int getWidth() const { return width; }
    int getHeight() const { return height; }

    // Add this new method to get a pointer to the buffer data
    const void *getBufferData() const
    {
        return buffer.get();
    }
};

float mapValue(float value, float fromLow, float fromHigh, float toLow, float toHigh)
{
    return (value - fromLow) * (toHigh - toLow) / (fromHigh - fromLow) + toLow;
}

inline __m256 julia_mandelbrot_simd(__m256 x0, __m256 y0, __m256 cr, __m256 ci, int maxIterations)
{
    __m256 x = x0;
    __m256 y = y0;
    __m256 x2 = _mm256_mul_ps(x, x);
    __m256 y2 = _mm256_mul_ps(y, y);
    __m256 iterations = _mm256_setzero_ps();
    __m256 four = _mm256_set1_ps(4.0f);
    __m256 one = _mm256_set1_ps(1.0f);
    __m256 mask = _mm256_setzero_ps();

    for (int i = 0; i < maxIterations; ++i)
    {
        mask = _mm256_cmp_ps(_mm256_add_ps(x2, y2), four, _CMP_LE_OS);
        if (_mm256_movemask_ps(mask) == 0)
            break;

        y = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(2.0f), x), y), ci);
        x = _mm256_add_ps(_mm256_sub_ps(x2, y2), cr);
        x2 = _mm256_mul_ps(x, x);
        y2 = _mm256_mul_ps(y, y);
        iterations = _mm256_add_ps(iterations, _mm256_and_ps(one, mask));
    }

    return iterations;
}

std::pair<float, float> calculateComplexConstant(float t, float cycleLength)
{
    // Ensure t is within [0, cycleLength]
    t = fmodf(t, cycleLength);

    // Map t to [0, 2π]
    float angle = (t / cycleLength) * 2 * M_PI;

    float cr = 0.7885f * std::cos(angle );
    float ci = 0.7885f * std::sin(angle ); 
    return {cr, ci};
}

void renderJuliaMandelbrot(SDL_Renderer *renderer, ImageBuffer &imageBuffer,
                           float centerX, float centerY, float zoom, int maxIterations,
                           int SCREEN_WIDTH, int SCREEN_HEIGHT, float time, float cycleLength)
{
    float aspectRatio = static_cast<float>(SCREEN_WIDTH) / SCREEN_HEIGHT;
    float rangeY = 4.0f / zoom;
    float rangeX = rangeY * aspectRatio;

    float minReal = centerX - rangeX / 2;
    float maxReal = centerX + rangeX / 2;
    float minImaginary = centerY - rangeY / 2;
    float maxImaginary = centerY + rangeY / 2;

    auto [cr, ci] = calculateComplexConstant(time, cycleLength);
    __m256 cr_vec = _mm256_set1_ps(cr);
    __m256 ci_vec = _mm256_set1_ps(ci);

#pragma omp parallel for schedule(dynamic, 1)
    for (int py = 0; py < SCREEN_HEIGHT; py++)
    {
        float y0 = mapValue(py, 0, SCREEN_HEIGHT - 1, minImaginary, maxImaginary);

        for (int px = 0; px < SCREEN_WIDTH; px += 8)
        {
            __m256 x0 = _mm256_setr_ps(
                mapValue(px, 0, SCREEN_WIDTH - 1, minReal, maxReal),
                mapValue(px + 1, 0, SCREEN_WIDTH - 1, minReal, maxReal),
                mapValue(px + 2, 0, SCREEN_WIDTH - 1, minReal, maxReal),
                mapValue(px + 3, 0, SCREEN_WIDTH - 1, minReal, maxReal),
                mapValue(px + 4, 0, SCREEN_WIDTH - 1, minReal, maxReal),
                mapValue(px + 5, 0, SCREEN_WIDTH - 1, minReal, maxReal),
                mapValue(px + 6, 0, SCREEN_WIDTH - 1, minReal, maxReal),
                mapValue(px + 7, 0, SCREEN_WIDTH - 1, minReal, maxReal));
            __m256 y0_vec = _mm256_set1_ps(y0);

            __m256 iterations = julia_mandelbrot_simd(x0, y0_vec, cr_vec, ci_vec, maxIterations);

            float iter_array[8];
            _mm256_storeu_ps(iter_array, iterations);

            for (int i = 0; i < 8 && px + i < SCREEN_WIDTH; ++i)
            {
                float iteration = iter_array[i];

                if (iteration < maxIterations)
                {
                    iteration = iteration + 1 - log(log(2.0f)) / log(2.0f);
                }

                SDL_Color color = getColor(iteration, maxIterations);

                float x0_scalar = mapValue(px + i, 0, SCREEN_WIDTH - 1, minReal, maxReal);
                float distanceFromCenter = sqrt((x0_scalar - centerX) * (x0_scalar - centerX) + (y0 - centerY) * (y0 - centerY));
                float gradientFactor = 1.0f - std::min(distanceFromCenter / (rangeX / 2), 1.0f);
                color.r = static_cast<Uint8>(color.r * gradientFactor);
                color.g = static_cast<Uint8>(color.g * gradientFactor);
                color.b = static_cast<Uint8>(std::min(255.0f, color.b * gradientFactor + 40));

                imageBuffer.at(px + i, py) = {color.r, color.g, color.b, color.a};
            }
        }
    }

    // Create texture and update it with the image buffer data
    SDL_Texture *texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, imageBuffer.getWidth(), imageBuffer.getHeight());
    SDL_UpdateTexture(texture, NULL, imageBuffer.getBufferData(), imageBuffer.getWidth() * sizeof(Pixel));
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_DestroyTexture(texture);
}

int main()
{
    SDL_Init(SDL_INIT_VIDEO);

    const int SCREEN_WIDTH = 800;
    const int SCREEN_HEIGHT = 600;

    SDL_Window *window = SDL_CreateWindow("Mandelbrot Zoom", SDL_WINDOWPOS_UNDEFINED,
                                          SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH,
                                          SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    float centerX = -0.3f;
    float centerY = 0.0f;
    float zoom = 2;
    float zoomSpeed = 1.01f;
    float moveSpeed = 0.1f;
    int maxIterations = 1000;

    bool quit = false;
    SDL_Event e;

    auto lastTime = std::chrono::high_resolution_clock::now();
    int frameCount = 0;

    int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    ImageBuffer imageBuffer(SCREEN_WIDTH, SCREEN_HEIGHT);

    float time = 0.0f;
    float timeSpeed = 0.003f;
    float cycleLength = 30.0f; // Length of one complete cycle in seconds

    while (!quit)
    {
        while (SDL_PollEvent(&e) != 0)
        {
            if (e.type == SDL_QUIT)
            {
                quit = true;
            }
            else if (e.type == SDL_KEYDOWN)
            {
                switch (e.key.keysym.sym)
                {
                case SDLK_UP: // Move up
                    centerY -= moveSpeed / zoom;
                    break;
                case SDLK_DOWN: // Move down
                    centerY += moveSpeed / zoom;
                    break;
                case SDLK_LEFT: // Move left
                    centerX -= moveSpeed / zoom;
                    break;
                case SDLK_RIGHT: // Move right
                    centerX += moveSpeed / zoom;
                    break;
                case SDLK_z: // Zoom in
                    zoom *= zoomSpeed;
                    break;
                case SDLK_x: // Zoom out
                    zoom /= zoomSpeed;
                    break;
                case SDLK_i: // Increase max iterations, means more detail
                    maxIterations += 100;
                    break;
                case SDLK_k: // Decrease max iterations, means less detail
                    maxIterations = std::max(100, maxIterations - 100);
                    break;
                case SDLK_c: // Increase cycle length
                    cycleLength += 1.0f;
                    break;
                case SDLK_v: // Decrease cycle length
                    cycleLength = std::max(1.0f, cycleLength - 1.0f);
                    break;
                }
            }
        }
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        auto start = std::chrono::high_resolution_clock::now();
        renderJuliaMandelbrot(renderer, imageBuffer, centerX, centerY, zoom, maxIterations, SCREEN_WIDTH, SCREEN_HEIGHT, time, cycleLength);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        SDL_RenderPresent(renderer);

        // Update time
        time += timeSpeed;

        // Calculate the current position in the cycle (0.0 to 1.0)
        float cyclePosition = fmodf(time, cycleLength) / cycleLength;
        // Calculate FPS
        frameCount++;
        auto currentTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float> elapsedTime = currentTime - lastTime;

        if (elapsedTime.count() >= 1.0f)
        {
            float fps = frameCount / elapsedTime.count();
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
                                      " - pos: " + std::to_string(cyclePosition);
            SDL_SetWindowTitle(window, windowTitle.c_str());
        }

        // Handle thread count changes
        const Uint8 *state = SDL_GetKeyboardState(NULL);
        for (int i = SDL_SCANCODE_1; i <= SDL_SCANCODE_8; ++i)
        {
            if (state[i])
            {
                numThreads = i - SDL_SCANCODE_1 + 1;
                omp_set_num_threads(numThreads);
                break;
            }
        }
        if (state[SDL_SCANCODE_COMMA])
            timeSpeed *= 0.99f; // Slow down
        if (state[SDL_SCANCODE_PERIOD])
            timeSpeed *= 1.01f; // Speed up
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}