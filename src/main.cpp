#include <SDL2/SDL.h>
#include <cmath>

float mapValue(float value, float fromLow, float fromHigh, float toLow, float toHigh) {
    return (value - fromLow) * (toHigh - toLow) / (fromHigh - fromLow) + toLow;
}

SDL_Color getColor(int iteration, int max_iteration) {
    if (iteration == max_iteration) return {0, 0, 0, 255};

    float hue = 360.0f * iteration / max_iteration;
    float saturation = 1.0f;
    float value = 1.0f;

    float c = value * saturation;
    float x = c * (1 - std::abs(std::fmod(hue / 60.0f, 2.0f) - 1));
    float m = value - c;

    float r, g, b;
    if (hue < 60) { r = c; g = x; b = 0; }
    else if (hue < 120) { r = x; g = c; b = 0; }
    else if (hue < 180) { r = 0; g = c; b = x; }
    else if (hue < 240) { r = 0; g = x; b = c; }
    else if (hue < 300) { r = x; g = 0; b = c; }
    else { r = c; g = 0; b = x; }

    return {
        static_cast<Uint8>((r + m) * 255),
        static_cast<Uint8>((g + m) * 255),
        static_cast<Uint8>((b + m) * 255),
        255
    };
}

// Based on pseudocode from Wikipedia.
// Reference: https://en.wikipedia.org/wiki/Mandelbrot_set#:~:text=The%20Mandelbrot%20set%20(%2F%CB%88m,aesthetic%20appeal%20and%20fractal%20structures.
void renderMandelbrot(SDL_Renderer* renderer, int SCREEN_WIDTH, int SCREEN_HEIGHT, 
                      float centerX, float centerY, float zoom) {
    float aspectRatio = static_cast<float>(SCREEN_WIDTH) / SCREEN_HEIGHT;
    float rangeY = 4.0f / zoom;
    float rangeX = rangeY * aspectRatio;
    
    float minReal = centerX - rangeX / 2;
    float maxReal = centerX + rangeX / 2;
    float minImaginary = centerY - rangeY / 2;
    float maxImaginary = centerY + rangeY / 2;
    
    int maxIterations = 1000;

    for (int px = 0; px < SCREEN_WIDTH; px++) {
        for (int py = 0; py < SCREEN_HEIGHT; py++) {
            // We map each pixel value to the Mandelbrot bounds
            float x0 = mapValue(px, 0, SCREEN_WIDTH - 1, minReal, maxReal);
            float y0 = mapValue(py, 0, SCREEN_HEIGHT - 1, minImaginary, maxImaginary);
            float x = 0.0f;
            float y = 0.0f;
            int iteration = 0;

            // We see its behaviour after many iterations
            while (x*x + y*y <= 4 && iteration < maxIterations) {
                float xtemp = x*x - y*y + x0;
                y = 2*x*y + y0;
                x = xtemp;
                iteration++;
            }

            // Based on the results, we set a color
            SDL_Color color = getColor(iteration, maxIterations);
            SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
            SDL_RenderDrawPoint(renderer, px, py);
        }
    }
}

int main(int argc, char* argv[]) {
    SDL_Init(SDL_INIT_VIDEO);
    
    const int SCREEN_WIDTH = 800;
    const int SCREEN_HEIGHT = 600;
    
    SDL_Window* window = SDL_CreateWindow("Mandelbrot Zoom", SDL_WINDOWPOS_UNDEFINED, 
                                          SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, 
                                          SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    
    float centerX = -0.5f;
    float centerY = 0.0f;
    float zoom = 1.0f;
    float zoomSpeed = 1.01f;
    
    bool quit = false;
    SDL_Event e;
    
    while (!quit) {
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) {
                quit = true;
            }
        }
        
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        
        renderMandelbrot(renderer, SCREEN_WIDTH, SCREEN_HEIGHT, centerX, centerY, zoom);
        
        SDL_RenderPresent(renderer);
        
        zoom *= zoomSpeed;
        SDL_Delay(16);  // Cap at roughly 60 FPS
    }
    
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return 0;
}