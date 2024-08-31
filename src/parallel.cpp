#include <SDL2/SDL.h>
#include <cmath>
#include <string>
#include <chrono>
#include <omp.h>
#include <vector>

float mapValue(float value, float fromLow, float fromHigh, float toLow, float toHigh) {
    return (value - fromLow) * (toHigh - toLow) / (fromHigh - fromLow) + toLow;
}

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

void renderMandelbrot(SDL_Renderer* renderer, int SCREEN_WIDTH, int SCREEN_HEIGHT,
                      float centerX, float centerY, float zoom, int maxIterations) {
    float aspectRatio = static_cast<float>(SCREEN_WIDTH) / SCREEN_HEIGHT;
    float rangeY = 4.0f / zoom;
    float rangeX = rangeY * aspectRatio;

    float minReal = centerX - rangeX / 2;
    float maxReal = centerX + rangeX / 2;
    float minImaginary = centerY - rangeY / 2;
    float maxImaginary = centerY + rangeY / 2;

    #pragma omp parallel
    {
        std::vector<SDL_Color> lineColors(SCREEN_WIDTH);

        #pragma omp for schedule(dynamic, 1)
        for (int py = 0; py < SCREEN_HEIGHT; py++) {
            for (int px = 0; px < SCREEN_WIDTH; px++) {
                float x0 = mapValue(px, 0, SCREEN_WIDTH - 1, minReal, maxReal);
                float y0 = mapValue(py, 0, SCREEN_HEIGHT - 1, minImaginary, maxImaginary);
                float x = 0.0f;
                float y = 0.0f;
                float iteration = 0.0f;
                float zn2 = 0.0f;

                while (zn2 <= 4 && iteration < maxIterations) {
                    float xtemp = x*x - y*y + x0;
                    y = 2*x*y + y0;
                    x = xtemp;
                    zn2 = x*x + y*y;
                    iteration += 1.0f;
                }

                if (iteration < maxIterations) {
                    iteration = iteration + 1 - log(log(sqrt(zn2))) / log(2.0);
                }

                SDL_Color color = getColor(iteration, maxIterations);

                float distanceFromCenter = sqrt((x0 - centerX) * (x0 - centerX) + (y0 - centerY) * (y0 - centerY));
                float gradientFactor = 1.0f - std::min(distanceFromCenter / (rangeX / 2), 1.0f);
                color.r = static_cast<Uint8>(color.r * gradientFactor);
                color.g = static_cast<Uint8>(color.g * gradientFactor);
                color.b = static_cast<Uint8>(std::min(255.0f, color.b * gradientFactor + 40));

                lineColors[px] = color;
            }

            #pragma omp critical
            {
                for (int px = 0; px < SCREEN_WIDTH; px++) {
                    SDL_SetRenderDrawColor(renderer, lineColors[px].r, lineColors[px].g, lineColors[px].b, lineColors[px].a);
                    SDL_RenderDrawPoint(renderer, px, py);
                }
            }
        }
    }
}

int main() {
    SDL_Init(SDL_INIT_VIDEO);

    const int SCREEN_WIDTH = 800;
    const int SCREEN_HEIGHT = 600;

    SDL_Window* window = SDL_CreateWindow("Mandelbrot Zoom", SDL_WINDOWPOS_UNDEFINED,
                                          SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH,
                                          SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    float centerX = -0.67f;
    float centerY = 0.0f;
    float zoom = 1.74f;
    float zoomSpeed = 1.01f;
    float moveSpeed = 0.1f;
    int maxIterations = 1000;

    bool quit = false;
    SDL_Event e;

    auto lastTime = std::chrono::high_resolution_clock::now();
    int frameCount = 0;

    int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    while (!quit) {
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) {
                quit = true;
            } else if (e.type == SDL_KEYDOWN) {
                switch (e.key.keysym.sym) {
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
                }
            }
        }

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        auto start = std::chrono::high_resolution_clock::now();
        renderMandelbrot(renderer, SCREEN_WIDTH, SCREEN_HEIGHT, centerX, centerY, zoom, maxIterations);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        SDL_RenderPresent(renderer);

        // Calculate FPS
        frameCount++;
        auto currentTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float> elapsedTime = currentTime - lastTime;

        if (elapsedTime.count() >= 1.0f) {
            float fps = frameCount / elapsedTime.count();
            frameCount = 0;
            lastTime = currentTime;

            std::string windowTitle = "fps: " + std::to_string(fps) +
                                      " - z: " + std::to_string(zoom) +
                                      " - y: " + std::to_string(centerY) +
                                      " - x: " + std::to_string(centerX) +
                                      " - mi: " + std::to_string(maxIterations) +
                                      " - threads: " + std::to_string(numThreads) +
                                      " - render time: " + std::to_string(elapsed.count()) + "s";
            SDL_SetWindowTitle(window, windowTitle.c_str());
        }

        // Handle thread count changes
        const Uint8* state = SDL_GetKeyboardState(NULL);
        for (int i = SDL_SCANCODE_1; i <= SDL_SCANCODE_8; ++i) {
            if (state[i]) {
                numThreads = i - SDL_SCANCODE_1 + 1;
                omp_set_num_threads(numThreads);
                break;
            }
        }
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
