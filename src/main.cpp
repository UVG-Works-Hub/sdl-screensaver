#include <SDL2/SDL.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>

const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;
const int NUM_CIRCLES = 20;

struct Circle {
    float x, y;
    float dx, dy;
    SDL_Color color;
};

int main(int argc, char* args[]) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cout << "SDL could not initialize! SDL_Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    SDL_Window* window = SDL_CreateWindow("Simple Screensaver", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    if (window == nullptr) {
        std::cout << "Window could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == nullptr) {
        std::cout << "Renderer could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    std::vector<Circle> circles(NUM_CIRCLES);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (auto& circle : circles) {
        circle.x = dis(gen) * SCREEN_WIDTH;
        circle.y = dis(gen) * SCREEN_HEIGHT;
        circle.dx = (dis(gen) - 0.5) * 5;
        circle.dy = (dis(gen) - 0.5) * 5;
        circle.color = {static_cast<Uint8>(dis(gen) * 255), 
                        static_cast<Uint8>(dis(gen) * 255), 
                        static_cast<Uint8>(dis(gen) * 255), 255};
    }

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

        for (auto& circle : circles) {
            circle.x += circle.dx;
            circle.y += circle.dy;

            if (circle.x < 0 || circle.x > SCREEN_WIDTH) circle.dx = -circle.dx;
            if (circle.y < 0 || circle.y > SCREEN_HEIGHT) circle.dy = -circle.dy;

            SDL_SetRenderDrawColor(renderer, circle.color.r, circle.color.g, circle.color.b, circle.color.a);
            for (int w = 0; w < 8; w++) {
                for (int h = 0; h < 8; h++) {
                    SDL_RenderDrawPoint(renderer, static_cast<int>(circle.x) + w, static_cast<int>(circle.y) + h);
                }
            }
        }

        SDL_RenderPresent(renderer);
        SDL_Delay(16);  // Cap at roughly 60 fps
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}